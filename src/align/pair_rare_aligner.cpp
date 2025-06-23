#include "rare_aligner.h"
#include "anchor.h"

PairRareAligner::PairRareAligner(const FilePath work_dir,
	const uint_t thread_num,
	uint_t chunk_size,
	uint_t overlap_size,
	uint_t min_anchor_length,
	uint_t max_anchor_frequency)
	: work_dir(work_dir)
	, index_dir(work_dir / INDEX_DIR)
	, thread_num(thread_num)
	, chunk_size(chunk_size)
	, overlap_size(overlap_size)
	, min_anchor_length(min_anchor_length)
	, max_anchor_frequency(max_anchor_frequency)
{
	if (!std::filesystem::exists(index_dir)) {
		std::filesystem::create_directories(index_dir);
	}


	this->group_id = 0;
	this->round_id = 0;

}

MatchVec3DPtr PairRareAligner::alignPairGenome(
	SpeciesName query_name,
	SeqPro::ManagerVariant& query_fasta_manager,
	SearchMode         search_mode,
	bool allow_MEM,
	sdsl::int_vector<0>& ref_global_cache,
	SeqPro::Length sampling_interval) {

	ThreadPool shared_pool(thread_num);

	MatchVec3DPtr anchors = findQueryFileAnchor(
		query_name, query_fasta_manager, search_mode, allow_MEM, shared_pool, ref_global_cache, sampling_interval);

	return anchors;

}

FilePath PairRareAligner::buildIndex(const std::string prefix, SeqPro::ManagerVariant& ref_fasta_manager_, bool fast_build) {

	FilePath ref_index_path = index_dir / prefix;

	if (!std::filesystem::exists(ref_index_path)) {
		std::filesystem::create_directories(ref_index_path);
	}

	ref_seqpro_manager = &ref_fasta_manager_;
	// 使用 std::visit 来获取 fasta_path
	std::string fasta_path_str;
	std::visit([&fasta_path_str](auto&& manager_ptr) {
		using PtrType = std::decay_t<decltype(manager_ptr)>;
		if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
			// 如果是 SequenceManager，直接调用 getFastaPath()
			fasta_path_str = manager_ptr->getFastaPath();
		} else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
			// 如果是 MaskedSequenceManager，先 getOriginalManager() 再 getFastaPath()
			fasta_path_str = manager_ptr->getOriginalManager().getFastaPath();
		} else {
			throw std::runtime_error("Unhandled manager type in variant.");
		}
	}, ref_fasta_manager_);

	// index_dir的路径加上prefix的前缀加上fasta_manager.fasta_path_的扩展名
	FilePath output_path = ref_index_path / (prefix + std::filesystem::path(fasta_path_str).extension().string());

	ref_index.emplace(prefix, ref_fasta_manager_);

	FilePath idx_file_path = ref_index_path / (prefix + "." + FMINDEX_EXTESION);

	spdlog::info("Indexing with prefix: {}, index path: {}", prefix, ref_index_path.string());

	if (!std::filesystem::exists(idx_file_path)) {
		ref_index->buildIndex(output_path, fast_build, thread_num);
		ref_index->saveToFile(idx_file_path.string());
	}
	else {
		ref_index->loadFromFile(idx_file_path.string());
	}
	spdlog::info("Indexing finished, index path: {}", ref_index_path.string());

	return ref_index_path;
}

MatchVec3DPtr PairRareAligner::findQueryFileAnchor(
	const std::string prefix,
	SeqPro::ManagerVariant& query_fasta_manager,
	SearchMode         search_mode,
	bool allow_MEM,
	ThreadPool& pool,
	sdsl::int_vector<0>& ref_global_cache,
	SeqPro::Length sampling_interval)
{
	/* ---------- 1. 结果文件路径，与多基因组保持同一目录 ---------- */
	FilePath result_dir = work_dir / RESULT_DIR
		/ ("group_" + std::to_string(group_id))
		/ ("round_" + std::to_string(round_id));
	std::filesystem::create_directories(result_dir);

	spdlog::info("[findQueryFileAnchor] begin to algin {}", prefix);

	FilePath anchor_file = result_dir /
		(prefix + "_" + SearchModeToString(search_mode) + "." + ANCHOR_EXTENSION);

	/* ---------- 若已存在结果文件，直接加载 ---------- */
	if (std::filesystem::exists(anchor_file)) {
		MatchVec3DPtr result = std::make_shared<MatchVec3D>();
		loadMatchVec3D(anchor_file, result);
		return result;
	}

	/* ---------- 读取 FASTA 并分片 ---------- */
	RegionVec chunks = preAllocateChunks(query_fasta_manager, chunk_size, overlap_size);
	// 智能分块策略：自动根据序列数量和长度选择最优的分块方式

	/* ---------- ① 计时：搜索 Anchor ---------- */
	auto t_search0 = std::chrono::steady_clock::now();

	// ThreadPool pool(thread_num);
	std::vector<std::future<MatchVec2DPtr>> futures;
	// 根据线程数量和chunk数量决定每个线程处理的chunk数量
	// 确保每个线程至少处理一个chunk，同时避免线程过多
	size_t num_chunks_per_thread = chunks.size() / thread_num;
	if (chunks.size() % thread_num != 0) {
		num_chunks_per_thread++; // 如果不能整除，则向上取整，确保所有chunk都被处理
	}
	if (num_chunks_per_thread == 0 && !chunks.empty()) { // 至少处理一个chunk
	    num_chunks_per_thread = 1;
	}

	futures.reserve(chunks.size() / num_chunks_per_thread + (chunks.size() % num_chunks_per_thread != 0 ? 1 : 0));

	for (size_t i = 0; i < chunks.size(); i += num_chunks_per_thread) {
		std::vector<Region> chunk_group;
		for (size_t j = i; j < std::min(i + num_chunks_per_thread, chunks.size()); ++j) {
			chunk_group.push_back(chunks[j]);
		}

		futures.emplace_back(
			pool.enqueue(
				[this, chunk_group, &query_fasta_manager, search_mode, allow_MEM, &ref_global_cache, sampling_interval]() -> MatchVec2DPtr {
					MatchVec2DPtr group_matches = std::make_shared<MatchVec2D>();
					for (const auto& ck : chunk_group) {
						std::string seq = std::visit([&ck](auto&& manager_ptr) -> std::string {
							using PtrType = std::decay_t<decltype(manager_ptr)>;
							if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
								return manager_ptr->getSubSequence(ck.chr_name, ck.start, ck.length);
							} else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
								return manager_ptr->getSubSequence(ck.chr_name, ck.start, ck.length);
							} else {
								throw std::runtime_error("Unhandled manager type in variant.");
							}
						}, query_fasta_manager);
						if (seq.length() <ck.length) continue;
						MatchVec2DPtr forwoard_matches = ref_index->findAnchors(
							ck.chr_name, seq, search_mode,
							Strand::FORWARD,
							allow_MEM,
							ck.start,
							0,
							max_anchor_frequency,
							ref_global_cache,
							sampling_interval);

						// 合并当前chunk的matches到group_matches
						for (const auto& match_list : *forwoard_matches) {
							group_matches->push_back(match_list);
						}	
					}
					return group_matches;
				}));
		futures.emplace_back(
			pool.enqueue(
				[this, chunk_group, &query_fasta_manager, search_mode, allow_MEM, &ref_global_cache, sampling_interval]() -> MatchVec2DPtr {
					MatchVec2DPtr group_matches = std::make_shared<MatchVec2D>();
					for (const auto& ck : chunk_group) {
						std::string seq = std::visit([&ck](auto&& manager_ptr) -> std::string {
							using PtrType = std::decay_t<decltype(manager_ptr)>;
							if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
								return manager_ptr->getSubSequence(ck.chr_name, ck.start, ck.length);
							} else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
								return manager_ptr->getSubSequence(ck.chr_name, ck.start, ck.length);
							} else {
								throw std::runtime_error("Unhandled manager type in variant.");
							}
						}, query_fasta_manager);
						if (seq.length() <ck.length) continue;
						MatchVec2DPtr reverse_matches = ref_index->findAnchors(
							ck.chr_name, seq, MIDDLE_SEARCH,
							Strand::REVERSE,
							allow_MEM,
							ck.start,
							min_anchor_length,
							max_anchor_frequency,
							ref_global_cache,
							sampling_interval);
						// 合并当前chunk的matches到group_matches
						for (const auto& match_list : *reverse_matches) {
							group_matches->push_back(match_list);
						}
					}
					return group_matches;
				}));

	}

	MatchVec3DPtr result = std::make_shared<MatchVec3D>();

	result->reserve(futures.size());

	for (auto& fut : futures) {
		MatchVec2DPtr part = fut.get();
		result->emplace_back(std::move(*part));
	}

	auto t_search1 = std::chrono::steady_clock::now();
	double search_ms = std::chrono::duration<double, std::milli>(t_search1 - t_search0).count();
	spdlog::info("[findQueryFileAnchor] search  = {:.3f} ms", search_ms);
	/* ---------- ③ 计时：保存 ---------- */
	auto t_save0 = std::chrono::steady_clock::now();
	// saveMatchVec3D(anchor_file, result);
	auto t_save1 = std::chrono::steady_clock::now();
	double save_ms = std::chrono::duration<double, std::milli>(t_save1 - t_save0).count();

	/* ---------- 打印统计 ---------- */
	spdlog::info("[findQueryFileAnchor] save    = {:.3f} ms", save_ms);

	return result;          // NRVO / move-elided
}


/* ============================================================= *
 *  把三维 clusters  ->  按 ref 的一维 clusters
 *  并行对每个 ref 走 keepWithSplitGreedy，再写回三维结构
 * ============================================================= */
 /**
  * @brief  全局 MatchClusterVec 上执行“两级 map + 最大堆贪婪拆分”过滤
  *
  * @param cluster_vec_ptr  所有 clusters 的 shared_ptr
  * @param pool             ThreadPool（本实现单线程，参数仅留作占位）
  * @param min_span         最小跨度阈值
  */
void PairRareAligner::Clusters2GraphByGreedy(SpeciesName query_name, MatchClusterVecPtr cluster_vec_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span)
{
	ThreadPool pool(thread_num);

	if (!cluster_vec_ptr || cluster_vec_ptr->empty()) return;

	/* ---------- 1. 取出所有 cluster 建最大堆 ---------- */
	struct Node {
		MatchCluster cl;
		int_t        span;
	};
	auto cmp = [](const Node& a, const Node& b) { return a.span < b.span; };

	std::vector<Node> heap;
	heap.reserve(cluster_vec_ptr->size());

	for (auto& cl : *cluster_vec_ptr) {
		if (cl.empty()) continue;
		int_t sc = clusterSpan(cl);
		if (sc >= min_span)
			heap.push_back({ std::move(cl), sc });
	}
	cluster_vec_ptr->clear();                 // 原向量腾空
	std::make_heap(heap.begin(), heap.end(), cmp);

	/* ---------- 2. 两级 interval map（按 ChrName） ---------- */
	std::unordered_map<ChrName, IntervalMap> rMaps;   // 参考端
	std::unordered_map<ChrName, IntervalMap> qMaps;   // 查询端

	MatchClusterVec kept; kept.reserve(heap.size());

	/* ---------- 3. 贪婪循环 ---------- */
	while (!heap.empty())
	{
		std::pop_heap(heap.begin(), heap.end(), cmp);
		Node cur = std::move(heap.back());
		heap.pop_back();

		if (cur.span < min_span || cur.cl.empty()) continue;

		const ChrName& refChr = cur.cl.front().ref_region.chr_name;
		const ChrName& qChr = cur.cl.front().query_region.chr_name;

		uint_t rb = start1(cur.cl.front());
		uint_t re = start1(cur.cl.back()) + len1(cur.cl.back());
		uint_t qb = start2(cur.cl.front());
		uint_t qe = start2(cur.cl.back()) + len2(cur.cl.back());

		int_t RL = 0, RR = 0, QL = 0, QR = 0;
		bool ref_hit = overlap1D(rMaps[refChr], rb, re, RL, RR);
		bool query_hit = overlap1D(qMaps[qChr], qb, qe, QL, QR);

		if (!ref_hit && !query_hit) {
			/* 完全无冲突 —— 记录并接纳 */
			insertInterval(rMaps[refChr], rb, re);
			insertInterval(qMaps[qChr], qb, qe);
			kept.emplace_back(std::move(cur.cl));
			continue;
		}

		/* 有冲突 —— 拆分为 ≤ 2 段，仍按跨度丢回堆 */
		for (auto& part : splitCluster(cur.cl, ref_hit, RL, RR,
			query_hit, QL, QR)) {
			int_t sc = clusterSpan(part);
			if (sc >= min_span) {
				heap.push_back({ std::move(part), sc });
				std::push_heap(heap.begin(), heap.end(), cmp);
			}
		}
	}

	/* ---------- 4. 输出过滤结果 ---------- */
	*cluster_vec_ptr = std::move(kept);      // 顺序无关
}


MatchClusterVecPtr PairRareAligner::filterPairSpeciesAnchors(SpeciesName query_name, MatchVec3DPtr& anchors, SeqPro::ManagerVariant& query_fasta_manager, RaMesh::RaMeshMultiGenomeGraph& graph)
{
	ThreadPool shared_pool(thread_num);
	MatchByStrandByQueryRefPtr unique_anchors = std::make_shared<MatchByStrandByQueryRef>();;
	MatchByStrandByQueryRefPtr repeat_anchors = std::make_shared<MatchByStrandByQueryRef>();;

	groupMatchByQueryRef(anchors, unique_anchors, repeat_anchors,
		*ref_seqpro_manager, query_fasta_manager, shared_pool);

	shared_pool.waitAllTasksDone();
	//anchors->clear();
	//anchors->shrink_to_fit();
	anchors.reset();
	spdlog::info("groupMatchByQueryRef done");
	sortMatchByQueryStart(unique_anchors, shared_pool);
	sortMatchByQueryStart(repeat_anchors, shared_pool);

	shared_pool.waitAllTasksDone();
	spdlog::info("sortMatchByQueryStart done");

	ClusterVecPtrByStrandByQueryRefPtr cluster_ptr = clusterAllChrMatch(
		unique_anchors,
		repeat_anchors,
		shared_pool);

	shared_pool.waitAllTasksDone();

	spdlog::info("clusterAllChrMatch done");

	MatchClusterVecPtr cluster_vec_ptr = groupClustersToVec(cluster_ptr, shared_pool, thread_num);
	
	shared_pool.waitAllTasksDone();

	return cluster_vec_ptr;

}


