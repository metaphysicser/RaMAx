#include "rare_aligner.h"

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
	FastaManager& query_fasta_manager,
	SearchMode         search_mode,
	bool allow_MEM) {

	ThreadPool shared_pool(thread_num);

	MatchVec3DPtr anchors = findQueryFileAnchor(
		query_name, query_fasta_manager, search_mode, allow_MEM, shared_pool);

	return anchors;

}

FilePath PairRareAligner::buildIndex(const std::string prefix, FastaManager& ref_fasta_manager_, bool fast_build) {

	FilePath ref_index_path = index_dir / prefix;

	if (!std::filesystem::exists(ref_index_path)) {
		std::filesystem::create_directories(ref_index_path);
	}

	ref_fasta_manager_ptr = &ref_fasta_manager_;

	// index_dir的路径加上prefix的前缀加上fasta_manager.fasta_path_的扩展名
	FilePath output_path = ref_index_path / (prefix + std::filesystem::path(ref_fasta_manager_ptr->fasta_path_).extension().string());

	ref_index = FM_Index(prefix, ref_fasta_manager_ptr);

	FilePath idx_file_path = ref_index_path / (prefix + "." + FMINDEX_EXTESION);

	spdlog::info("Indexing with prefix: {}, index path: {}", prefix, ref_index_path.string());

	if (!std::filesystem::exists(idx_file_path)) {
		ref_index.buildIndex(output_path, fast_build, thread_num);
		ref_index.saveToFile(idx_file_path.string());
	}
	else {
		ref_index.loadFromFile(idx_file_path.string());
	}
	spdlog::info("Indexing finished, index path: {}", ref_index_path.string());

	return ref_index_path;
}

MatchVec3DPtr PairRareAligner::findQueryFileAnchor(
	const std::string prefix,
	FastaManager& query_fasta_manager,
	SearchMode         search_mode,
	bool allow_MEM,
	ThreadPool& pool)
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
	RegionVec chunks = query_fasta_manager.preAllocateChunks(chunk_size, overlap_size);

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
				[this, chunk_group, &query_fasta_manager, search_mode, allow_MEM]() -> MatchVec2DPtr {
					MatchVec2DPtr group_matches = std::make_shared<MatchVec2D>();
					for (const auto& ck : chunk_group) {
						std::string seq = query_fasta_manager.getSubSequence(ck.chr_name, ck.start, ck.length);
						MatchVec2DPtr forwoard_matches = ref_index.findAnchors(
							ck.chr_name, seq, search_mode,
							Strand::FORWARD,
							allow_MEM,
							ck.start,
							0,
							max_anchor_frequency);

						// 合并当前chunk的matches到group_matches
						for (const auto& match_list : *forwoard_matches) {
							group_matches->push_back(match_list);
						}	
					}
					return group_matches;
				}));
		futures.emplace_back(
			pool.enqueue(
				[this, chunk_group, &query_fasta_manager, search_mode, allow_MEM]() -> MatchVec2DPtr {
					MatchVec2DPtr group_matches = std::make_shared<MatchVec2D>();
					for (const auto& ck : chunk_group) {
						std::string seq = query_fasta_manager.getSubSequence(ck.chr_name, ck.start, ck.length);
						
						MatchVec2DPtr reverse_matches = ref_index.findAnchors(
							ck.chr_name, seq, MIDDLE_SEARCH,
							Strand::REVERSE,
							allow_MEM,
							ck.start,
							min_anchor_length,
							max_anchor_frequency);
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
	saveMatchVec3D(anchor_file, result);
	auto t_save1 = std::chrono::steady_clock::now();
	double save_ms = std::chrono::duration<double, std::milli>(t_save1 - t_save0).count();

	/* ---------- 打印统计 ---------- */
	spdlog::info("[findQueryFileAnchor] save    = {:.3f} ms", save_ms);

	return result;          // NRVO / move-elided
}



//--------------------------------------------------------------------
// 主函数：直接将 slice 写入 (queryIdx, refIdx) 桶
//--------------------------------------------------------------------
void PairRareAligner::filterPairSpeciesAnchors(MatchVec3DPtr& anchors,
	FastaManager& query_fasta_manager)
{
	//MatchByQueryRef unique_anchors;
	//MatchByQueryRef repeat_anchors;
	MatchByStrandByQueryRef unique_anchors;
	MatchByStrandByQueryRef repeat_anchors;

	groupMatchByQueryRef(anchors, unique_anchors, repeat_anchors,
		*ref_fasta_manager_ptr, query_fasta_manager, thread_num);

	//// 测试得到的anchor对应的子串是否相同
	//for (const auto& first_match : unique_anchors[0][0]) {
	//		std::string query_subseq = query_fasta_manager.getSubSequence(
	//			first_match.query_region.chr_name, first_match.query_region.start, first_match.query_region.length);
	//		std::string ref_subseq = ref_fasta_manager.getSubSequence(
	//			first_match.ref_region.chr_name, first_match.ref_region.start, first_match.ref_region.length);
	//		if (query_subseq != ref_subseq) {
	//			spdlog::warn("Mismatch found in anchor: query {} vs ref {}",
	//				query_subseq, ref_subseq);
	//
	//	}
	//}

	sortMatchByQueryStart(unique_anchors, thread_num);
	sortMatchByQueryStart(repeat_anchors, thread_num);
	ThreadPool pool(thread_num);

	//--------------------------------------------------------------------
// ⬇ 1. 额外的 3-D 结果桶，和 unique_anchors 同形
//--------------------------------------------------------------------
	ClusterVecPtrByStrandByQueryRef cluster_results;
	cluster_results.resize(2);                      // strand: 0 = FORWARD, 1 = REVERSE

	for (auto& query_layer : cluster_results) {
		query_layer.resize(unique_anchors.size());          // 每条 strand 下的所有 query-chr
		for (auto& ref_row : query_layer)
			ref_row.resize(unique_anchors.front().size());            // 每个 (strand, query) 下的所有 ref-chr
	}

	//--------------------------------------------------------------------
	// 2. 主循环 —— 为每个 (i,j) 任务预留位置并提交线程池
	//--------------------------------------------------------------------
	for (uint_t k = 0; k < 2; k++) {
		for (uint_t i = 0; i < unique_anchors.size(); ++i) {
			for (uint_t j = 0; j < unique_anchors[i].size(); ++j) {
					MatchVec& unique_vec = unique_anchors[k][i][j];
					MatchVec& repeat_vec = repeat_anchors[k][i][j];

					// 先把下标复制到局部，避免 lambda 捕获引用后被后续循环修改
					auto kk = k;
					auto ii = i;
					auto jj = j;

					pool.enqueue([&cluster_results, kk, ii, jj,
						&unique_vec, &repeat_vec]()
						{
							// ⬇ 收集返回值
							cluster_results[kk][ii][jj] = clusterChrMatch(
								unique_vec,
								repeat_vec);   // 已排序
						});
			}
		}
	}

	pool.waitAllTasksDone();
//
//	//----------------------------------------------------------------
//	// 3) 如有后续操作，可直接用 cluster_results
//	//----------------------------------------------------------------
//	// anchors->clear(); …                // 原来的清理逻辑保持不变
//
//
//	//----------------------------------------------------------------
//	// 3) 释放原始 3D 数据以节省内存
//	//----------------------------------------------------------------
//	anchors->clear();
//	anchors->shrink_to_fit();
	return;

}


