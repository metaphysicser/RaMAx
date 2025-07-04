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

	ref_name = prefix;
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
void PairRareAligner::constructGraphByGreedy(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, ClusterVecPtrByStrandByQueryRefPtr cluster_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, uint_t min_span)
{
#ifdef _DEBUG_
	using Clock = std::chrono::steady_clock;
	using nsec_t = uint64_t;
	static std::atomic<nsec_t> ns_producer{ 0 }, ns_extend{ 0 }, ns_insert{ 0 };
	static std::atomic<uint64_t> cnt_submit{ 0 }, cnt_extend{ 0 }, cnt_insert{ 0 }, cnt_finish{ 0 };

	static std::atomic<nsec_t> ns_check_overlap{ 0 }, ns_insert_interval{ 0 },
		ns_clone_cluster{ 0 }, ns_enqueue_task{ 0 }, ns_store_kept{ 0 };

	auto ns2ms = [](nsec_t ns) { return static_cast<double>(ns) / 1'000'000.0; };
#endif

	ThreadPool pool(thread_num);

	MatchClusterVecPtr cluster_vec_ptr = groupClustersToVec(cluster_ptr, pool, thread_num);

	pool.waitAllTasksDone();

	if (!cluster_vec_ptr || cluster_vec_ptr->empty()) return;

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
	cluster_vec_ptr->clear();
	std::make_heap(heap.begin(), heap.end(), cmp);

	std::unordered_map<ChrName, IntervalMap> rMaps;
	std::unordered_map<ChrName, IntervalMap> qMaps;

	MatchClusterVec kept;
	kept.reserve(heap.size());

	while (!heap.empty()) {
#ifdef _DEBUG_
		auto t0 = Clock::now();
#endif

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

#ifdef _DEBUG_
		auto t_chk = Clock::now();
#endif
		bool ref_hit = overlap1D(rMaps[refChr], rb, re, RL, RR);
		bool query_hit = overlap1D(qMaps[qChr], qb, qe, QL, QR);
#ifdef _DEBUG_
		ns_check_overlap += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_chk).count();
#endif

		if (!ref_hit && !query_hit) {
#ifdef _DEBUG_
			auto t_ins = Clock::now();
#endif
			insertInterval(rMaps[refChr], rb, re);
			insertInterval(qMaps[qChr], qb, qe);
#ifdef _DEBUG_
			ns_insert_interval += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_ins).count();
#endif

#ifdef _DEBUG_
			auto t_clone = Clock::now();
#endif
			auto task_cl = std::make_shared<MatchCluster>(cur.cl);
#ifdef _DEBUG_
			ns_clone_cluster += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_clone).count();
#endif

#ifdef _DEBUG_
			auto t_enq = Clock::now();
#endif
			pool.enqueue([this, &graph, &query_name, task_cl, &query_seqpro_manager] {
#ifdef _DEBUG_
				auto t1 = Clock::now();
#endif
				// AnchorVec anchor_vec = extendClusterToAnchor(*task_cl, *ref_seqpro_manager, query_seqpro_manager);
#ifdef _DEBUG_
				auto t2 = Clock::now();
				ns_extend += std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
				cnt_extend.fetch_add(1, std::memory_order_relaxed);
#endif
				// graph.insertAnchorIntoGraph(ref_name, query_name, anchor_vec);
				graph.insertClusterIntoGraph(ref_name, query_name, *task_cl);
#ifdef _DEBUG_
				auto t3 = Clock::now();
				ns_insert += std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count();
				cnt_insert.fetch_add(1, std::memory_order_relaxed);
				cnt_finish.fetch_add(1, std::memory_order_relaxed);
#endif
				});
#ifdef _DEBUG_
			ns_enqueue_task += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_enq).count();
#endif

#ifdef _DEBUG_
			auto t_store = Clock::now();
#endif
			kept.emplace_back(cur.cl);
#ifdef _DEBUG_
			ns_store_kept += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_store).count();
#endif

#ifdef _DEBUG_
			ns_producer += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t0).count();
			cnt_submit.fetch_add(1, std::memory_order_relaxed);
#endif

			continue;
		}

		for (auto& part : splitCluster(cur.cl, ref_hit, RL, RR, query_hit, QL, QR)) {
			int_t sc = clusterSpan(part);
			if (sc >= min_span) {
				heap.push_back({ part, sc });
				std::push_heap(heap.begin(), heap.end(), cmp);
			}
		}
	}

	pool.waitAllTasksDone();

#ifdef _DEBUG_
	std::cerr << "\n==========  DEBUG: constructGraphByGreedy  ==========\n"
		<< "Producer (split/enqueue) : " << ns2ms(ns_producer.load()) << " ms\n"
		<< "extendClusterToAnchor    : " << ns2ms(ns_extend.load()) << " ms  ("
		<< cnt_extend << " calls, avg "
		<< ns2ms(ns_extend.load()) / std::max<uint64_t>(1, cnt_extend) << " ms)\n"
		<< "insertAnchorIntoGraph    : " << ns2ms(ns_insert.load()) << " ms  ("
		<< cnt_insert << " calls, avg "
		<< ns2ms(ns_insert.load()) / std::max<uint64_t>(1, cnt_insert) << " ms)\n"
		<< "Tasks submitted / finished : "
		<< cnt_submit << " / " << cnt_finish << '\n'
		<< "----- Producer step breakdown (ms) -----\n"
		<< "check overlap         : " << ns2ms(ns_check_overlap) << '\n'
		<< "insertInterval        : " << ns2ms(ns_insert_interval) << '\n'
		<< "clone cluster (shared): " << ns2ms(ns_clone_cluster) << '\n'
		<< "enqueue task          : " << ns2ms(ns_enqueue_task) << '\n'
		<< "store kept            : " << ns2ms(ns_store_kept) << '\n'
		<< "=====================================================\n";
#endif

	// *cluster_vec_ptr = std::move(kept);
}



ClusterVecPtrByStrandByQueryRefPtr PairRareAligner::filterPairSpeciesAnchors(SpeciesName query_name, MatchVec3DPtr& anchors, SeqPro::ManagerVariant& query_fasta_manager, RaMesh::RaMeshMultiGenomeGraph& graph)
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

	return cluster_ptr;

}

void PairRareAligner::constructGraphByDP(
	SpeciesName                               query_name,
	SeqPro::ManagerVariant& query_seqpro_manager,
	ClusterVecPtrByStrandByQueryRefPtr        matrix_ptr,
	RaMesh::RaMeshMultiGenomeGraph& graph,
	uint_t                                    min_span)
{
#ifdef _DEBUG_
	using Clock = std::chrono::steady_clock;
	using nsec_t = uint64_t;
	static std::atomic<nsec_t> ns_scan{ 0 }, ns_rebal{ 0 }, ns_enq{ 0 }, ns_ext{ 0 }, ns_ins{ 0 };
	static std::atomic<uint64_t> cnt_clu{ 0 }, cnt_tasks{ 0 };
	auto ns2ms = [](nsec_t ns) {return static_cast<double>(ns) / 1'000'000.0;};
#endif
	if (!matrix_ptr || matrix_ptr->empty()) return;
	const size_t T = thread_num ? thread_num : 1;

	/* ---------- 0. 轻量哈希 ---------- */
	auto hash_qr = [](size_t q, size_t r) noexcept->size_t {
		size_t h = q + 0x9e3779b97f4a7c15ull;
		h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ull; h ^= (h >> 27);
		h = (h * 0x94d049bb133111ebull) ^ (h >> 31);
		return h ^ (r + (r << 6) + (r >> 2));
		};

	/* ---------- 1. 扫描+分桶 ---------- */
#ifdef _DEBUG_
	auto t0 = Clock::now();
#endif
	struct Bin {
		std::vector<MatchClusterVecPtr> vecs;
		size_t clusters = 0;
	};
	std::vector<Bin> bins(T);

	for (size_t s = 0;s < matrix_ptr->size();++s) {
		auto& qArr = (*matrix_ptr)[s];
		for (size_t q = 0;q < qArr.size();++q) {
			auto& rArr = qArr[q];
			for (size_t r = 0;r < rArr.size();++r) {
				auto& ptr = rArr[r];
				if (!ptr || ptr->empty())continue;
				size_t b = hash_qr(q, r) % T;
				bins[b].vecs.push_back(ptr);
				bins[b].clusters += ptr->size();
			}
		}
	}
#ifdef _DEBUG_
	ns_scan += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t0).count();
#endif

	/* ---------- 2. rebalance：最大桶按 cluster 数折半填空桶 ---------- */
/* ---------- 2. rebalance：最大桶深度折半填空桶 ---------- */
#ifdef _DEBUG_
	auto t_reb0 = Clock::now();
#endif
	{
		std::vector<size_t> empty;
		for (size_t i = 0; i < bins.size(); ++i)
			if (bins[i].clusters == 0) empty.push_back(i);

		// 最大堆 (clusters, idx)
		using Info = std::pair<size_t, size_t>;
		auto cmp = [](const Info& a, const Info& b) { return a.first < b.first; };

		auto build_heap = [&] {
			std::priority_queue<Info, std::vector<Info>, decltype(cmp)> pq(cmp);
			for (size_t i = 0; i < bins.size(); ++i)
				if (bins[i].clusters > 1) pq.emplace(bins[i].clusters, i);
			return pq;
			};

		auto pq = build_heap();

		while (!empty.empty() && !pq.empty())
		{
			auto [clu, idx] = pq.top(); pq.pop();
			if (clu <= 1) break;                         // 已无法再切
			Bin& big = bins[idx];

			// 期望搬走的一半 cluster
			size_t target = clu / 2;
			size_t moved = 0;
			Bin newBin;

			// ① 先按 pointer 搬
			while (!big.vecs.empty() && moved + big.vecs.back()->size() <= target) {
				auto ptr = std::move(big.vecs.back());
				big.vecs.pop_back();
				big.clusters -= ptr->size();
				newBin.clusters += ptr->size();
				newBin.vecs.push_back(std::move(ptr));
				moved += newBin.vecs.back()->size();
			}

			// ② 若还缺，且 big 里至少有 1 个 vecPtr，可深度拆 vec 本身
			if (moved < target && !big.vecs.empty()) {
				// 取最后一个 vecPtr（最大概率大）
				auto& src_ptr = big.vecs.back();
				if (src_ptr->size() > 1) {
					size_t need = target - moved;
					size_t half = std::min(need, src_ptr->size() / 2);
					auto split_it = src_ptr->end() - half;             // 从尾部切

					auto new_ptr = std::make_shared<MatchClusterVec>(split_it, src_ptr->end());
					src_ptr->erase(split_it, src_ptr->end());

					newBin.vecs.push_back(new_ptr);
					newBin.clusters += new_ptr->size();
					big.clusters -= new_ptr->size();
					moved += new_ptr->size();
				}
			}

			// 把 newBin 填入空槽
			if (newBin.clusters == 0) break;              // 防守：没切成功
			size_t dst = empty.back(); empty.pop_back();
			bins[dst] = std::move(newBin);

			// 把分裂后的桶重新入堆
			pq = build_heap();                            // 重新建堆保持最大性
		}
	}
#ifdef _DEBUG_
	ns_rebal += std::chrono::duration_cast<std::chrono::nanoseconds>(
		Clock::now() - t_reb0).count();
#endif


	/* ---------- 3. 并行处理 ---------- */
	ThreadPool pool(T);
	for (auto& bin : bins) {
		if (bin.clusters == 0) continue;
#ifdef _DEBUG_
		cnt_tasks.fetch_add(1, std::memory_order_relaxed);
		auto t_en = Clock::now();
#endif
		pool.enqueue([this, &graph, &query_seqpro_manager, &query_name,
			clusters = std::move(bin.vecs), min_span]() mutable
			{
				for (auto& vecPtr : clusters) {
					if (!vecPtr) continue;
					for (const auto& clu : *vecPtr) {
						if (clu.empty() || clusterSpan(clu) < static_cast<int_t>(min_span)) continue;
#ifdef _DEBUG_
						cnt_clu.fetch_add(1, std::memory_order_relaxed);
						auto tA = Clock::now();
#endif
						AnchorVec anchors = extendClusterToAnchor(clu, *ref_seqpro_manager, query_seqpro_manager);
#ifdef _DEBUG_
						auto tB = Clock::now();
						ns_ext += std::chrono::duration_cast<std::chrono::nanoseconds>(tB - tA).count();
#endif
						graph.insertAnchorIntoGraph(clu.front().ref_region.chr_name, query_name, anchors);
#ifdef _DEBUG_
						auto tC = Clock::now();
						ns_ins += std::chrono::duration_cast<std::chrono::nanoseconds>(tC - tB).count();
#endif
					}
				}
			});
#ifdef _DEBUG_
		ns_enq += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t_en).count();
#endif
	}
	pool.waitAllTasksDone();

#ifdef _DEBUG_
	spdlog::debug("\n========== DEBUG: constructGraphByDP ==========\n"
		"scan + binning        : {:.3f} ms\n"
		"rebalance split       : {:.3f} ms\n"
		"enqueue tasks ({})    : {:.3f} ms\n"
		"extendClusterToAnchor : {:.3f} ms ({} clusters, avg {:.3f} ms)\n"
		"insertAnchorIntoGraph : {:.3f} ms ({} clusters, avg {:.3f} ms)\n"
		"===========================================================\n",
		ns2ms(ns_scan), ns2ms(ns_rebal),
		cnt_tasks.load(), ns2ms(ns_enq),
		ns2ms(ns_ext), cnt_clu.load(),
		ns2ms(ns_ext) / std::max<uint64_t>(1, cnt_clu.load()),
		ns2ms(ns_ins), cnt_clu.load(),
		ns2ms(ns_ins) / std::max<uint64_t>(1, cnt_clu.load())
	);
#endif
}








