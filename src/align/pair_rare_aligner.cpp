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

	/* ---------- ③ 计时：保存 ---------- */
	auto t_save0 = std::chrono::steady_clock::now();
	// saveMatchVec3D(anchor_file, result);
	auto t_save1 = std::chrono::steady_clock::now();
	double save_ms = std::chrono::duration<double, std::milli>(t_save1 - t_save0).count();

	/* ---------- Performance Statistics ---------- */
	spdlog::info("");
	spdlog::info("┌─────────────────────────────────────────────────────────┐");
	spdlog::info("│               findQueryFileAnchor Performance           │");
	spdlog::info("├─────────────────────────────────────────────────────────┤");
	spdlog::info("│  Search phase    : {:>8.3f} ms                          │", search_ms);
	spdlog::info("│  Save phase      : {:>8.3f} ms                          │", save_ms);
	spdlog::info("│  Total time      : {:>8.3f} ms                          │", search_ms + save_ms);
	spdlog::info("└─────────────────────────────────────────────────────────┘");
	spdlog::info("");

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

		Strand strand = cur.cl.front().strand;
		uint_t rb = start1(cur.cl.front());
		uint_t re = start1(cur.cl.back()) + len1(cur.cl.back());
		uint_t qb = 0;
		uint_t qe = 0;
		if (strand == FORWARD) {
			qb = start2(cur.cl.front());
			qe = start2(cur.cl.back()) + len2(cur.cl.back());
		}
		else {
			qb = start2(cur.cl.back());
			qe = start2(cur.cl.front()) + len2(cur.cl.front());
		}

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
				AnchorVec anchor_vec = extendClusterToAnchor(*task_cl, *ref_seqpro_manager, query_seqpro_manager);
				graph.insertAnchorIntoGraph(ref_name, query_name, anchor_vec);
				// graph.insertClusterIntoGraph(ref_name, query_name, *task_cl);
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
	spdlog::debug("");
	spdlog::debug("================================================================");
	spdlog::debug("                  constructGraphByGreedy Performance            ");
	spdlog::debug("================================================================");
	spdlog::debug("Main Operations:");
	spdlog::debug("  Producer (split/enqueue)  : {:>10.3f} ms", ns2ms(ns_producer.load()));
	spdlog::debug("  extendClusterToAnchor     : {:>10.3f} ms ({} calls, avg {:.3f} ms)",
	              ns2ms(ns_extend.load()), cnt_extend.load(),
	              ns2ms(ns_extend.load()) / std::max<uint64_t>(1, cnt_extend));
	spdlog::debug("  insertAnchorIntoGraph     : {:>10.3f} ms ({} calls, avg {:.3f} ms)",
	              ns2ms(ns_insert.load()), cnt_insert.load(),
	              ns2ms(ns_insert.load()) / std::max<uint64_t>(1, cnt_insert));
	spdlog::debug("");
	spdlog::debug("Task Statistics:");
	spdlog::debug("  Tasks submitted           : {}", cnt_submit.load());
	spdlog::debug("  Tasks finished            : {}", cnt_finish.load());
	spdlog::debug("");
	spdlog::debug("Producer Step Breakdown:");
	spdlog::debug("  check overlap             : {:>10.3f} ms", ns2ms(ns_check_overlap));
	spdlog::debug("  insertInterval            : {:>10.3f} ms", ns2ms(ns_insert_interval));
	spdlog::debug("  clone cluster (shared)    : {:>10.3f} ms", ns2ms(ns_clone_cluster));
	spdlog::debug("  enqueue task              : {:>10.3f} ms", ns2ms(ns_enqueue_task));
	spdlog::debug("  store kept                : {:>10.3f} ms", ns2ms(ns_store_kept));
	spdlog::debug("================================================================");
	spdlog::debug("");
#endif

	// *cluster_vec_ptr = std::move(kept);
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
void PairRareAligner::constructGraphByGreedyByRef(SpeciesName query_name, SeqPro::ManagerVariant& query_seqpro_manager, MatchClusterVecPtr cluster_vec_ptr, RaMesh::RaMeshMultiGenomeGraph& graph, ThreadPool& pool, uint_t min_span)
{
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

		std::pop_heap(heap.begin(), heap.end(), cmp);
		Node cur = std::move(heap.back());
		heap.pop_back();

		if (cur.span < min_span || cur.cl.empty()) continue;

		const ChrName& refChr = cur.cl.front().ref_region.chr_name;
		const ChrName& qChr = cur.cl.front().query_region.chr_name;

		Strand strand = cur.cl.front().strand;
		uint_t rb = start1(cur.cl.front());
		uint_t re = start1(cur.cl.back()) + len1(cur.cl.back());
		uint_t qb = 0;
		uint_t qe = 0;
		if (strand == FORWARD) {
			qb = start2(cur.cl.front());
			qe = start2(cur.cl.back()) + len2(cur.cl.back());
		}
		else {
			qb = start2(cur.cl.back());
			qe = start2(cur.cl.front()) + len2(cur.cl.front());
		}

		int_t RL = 0, RR = 0, QL = 0, QR = 0;

		bool ref_hit = overlap1D(rMaps[refChr], rb, re, RL, RR);
		bool query_hit = overlap1D(qMaps[qChr], qb, qe, QL, QR);

		if (!ref_hit && !query_hit) {
			insertInterval(rMaps[refChr], rb, re);
			insertInterval(qMaps[qChr], qb, qe);
			auto task_cl = std::make_shared<MatchCluster>(cur.cl);

			pool.enqueue([this, &graph, &query_name, task_cl, &query_seqpro_manager] {
				AnchorVec anchor_vec = extendClusterToAnchor(*task_cl, *ref_seqpro_manager, query_seqpro_manager);
				graph.insertAnchorIntoGraph(ref_name, query_name, anchor_vec);
				// graph.insertClusterIntoGraph(ref_name, query_name, *task_cl);
				});

			kept.emplace_back(cur.cl);

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








