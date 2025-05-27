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

}


std::vector<uint64_t> read_sa(const std::string& sa_file)
{
	std::ifstream file(sa_file, std::ios::binary | std::ios::ate);
	if (!file) throw std::runtime_error("Cannot open SA file");
	std::streamsize size = file.tellg();
	file.seekg(0);

	if (size % 5 != 0) {
		throw std::runtime_error("SA file size not multiple of 5 bytes");
	}

	size_t count = size / 5;
	std::vector<uint64_t> sa(count);
	std::vector<char> buffer(size);
	file.read(buffer.data(), size);

	const uint8_t* ptr = reinterpret_cast<const uint8_t*>(buffer.data());
	for (size_t i = 0; i < count; i++, ptr += 5)
	{
		// Little-Endian interpretation:
		uint64_t val =
			(uint64_t(ptr[0])) |
			(uint64_t(ptr[1]) << 8) |
			(uint64_t(ptr[2]) << 16) |
			(uint64_t(ptr[3]) << 24) |
			(uint64_t(ptr[4]) << 32);
		sa[i] = val;
	}
	return sa;
}

FilePath PairRareAligner::buildIndex(const std::string prefix, const FilePath fasta_path, bool fast_build) {

	FilePath ref_index_path = index_dir / prefix;

	if (!std::filesystem::exists(ref_index_path)) {
		std::filesystem::create_directories(ref_index_path);
	}


	ref_fasta_manager = FastaManager(fasta_path, getFaiIndexPath(fasta_path));

	// index_dir的路径加上prefix的前缀加上fasta_manager.fasta_path_的扩展名
	FilePath output_path = ref_index_path / (prefix + std::filesystem::path(fasta_path).extension().string());

	ref_index = FM_Index(prefix, &ref_fasta_manager);

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
	bool allow_MEM)
{
	namespace ch = std::chrono;
	FilePath result_dir_path = work_dir / RESULT_DIR;
	std::filesystem::create_directories(result_dir_path);

	FilePath anchor_file = result_dir_path /
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
	auto t_search0 = ch::steady_clock::now();

	ThreadPool pool(thread_num);
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
						MatchVec2DPtr matches = ref_index.findAnchors(
							ck.chr_name, seq, search_mode,
							Strand::FORWARD,
							allow_MEM,
							ck.start,
							min_anchor_length,
							max_anchor_frequency);
						// 合并当前chunk的matches到group_matches
						for(const auto& match_list : *matches){
						    group_matches->push_back(match_list);
						}
					}
					return group_matches;
				}));
	}

	pool.waitAllTasksDone();

	auto t_search1 = ch::steady_clock::now();
	double search_ms = ch::duration<double, std::milli>(t_search1 - t_search0).count();
	spdlog::info("[findQueryFileAnchor] search  = {:.3f} ms", search_ms);
	/* ---------- ② 计时：合并 ---------- */
	auto t_merge0 = ch::steady_clock::now();

	MatchVec3DPtr result = std::make_shared<MatchVec3D>();

	result->reserve(futures.size());
	size_t total_lists = 0;
	for (auto& fut : futures) {
		MatchVec2DPtr part = fut.get();
		total_lists += part->size();
		result->emplace_back(std::move(*part));
	}

	auto t_merge1 = ch::steady_clock::now();
	double merge_ms = ch::duration<double, std::milli>(t_merge1 - t_merge0).count();
	spdlog::info("[findQueryFileAnchor] merge   = {:.3f} ms", merge_ms);
	/* ---------- ③ 计时：保存 ---------- */
	auto t_save0 = ch::steady_clock::now();
	saveMatchVec3D(anchor_file, result);
	auto t_save1 = ch::steady_clock::now();
	double save_ms = ch::duration<double, std::milli>(t_save1 - t_save0).count();

	/* ---------- 打印统计 ---------- */
	spdlog::info("[findQueryFileAnchor] save    = {:.3f} ms", save_ms);

	return result;          // NRVO / move-elided
}



//--------------------------------------------------------------------
// 主函数：直接将 slice 写入 (queryIdx, refIdx) 桶
//--------------------------------------------------------------------
void PairRareAligner::clusterPairSpeciesAnchors(MatchVec3DPtr& anchors,
	FastaManager& query_fasta_manager)
{
	MatchByQueryRef unique_anchors;
	MatchByQueryRef repeat_anchors;

	groupMatchByQueryRef(anchors, unique_anchors, repeat_anchors,
		ref_fasta_manager, query_fasta_manager, thread_num);

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
	ClusterVecPtrByRefByQuery cluster_results(unique_anchors.size());
	for (auto& row : cluster_results)
		row.resize(/*与这一行对应的 ref-bucket 个数*/ 0);   // 稍后再 resize

	//--------------------------------------------------------------------
	// 2. 主循环 —— 为每个 (i,j) 任务预留位置并提交线程池
	//--------------------------------------------------------------------
	for (uint_t i = 0; i < unique_anchors.size(); ++i) {
		cluster_results[i].resize(unique_anchors[i].size());    // ← 现在知道列数

		for (uint_t j = 0; j < unique_anchors[i].size(); ++j) {
			MatchVec& unique_vec = unique_anchors[i][j];
			MatchVec& repeat_vec = repeat_anchors[i][j];

			// 先把下标复制到局部，避免 lambda 捕获引用后被后续循环修改
			auto ii = i;
			auto jj = j;

			pool.enqueue([&cluster_results, ii, jj,
				&unique_vec, &repeat_vec]()
				{
					// ⬇ 收集返回值
					cluster_results[ii][jj] = clusterChrMatch(
						unique_vec,
						repeat_vec);   // 已排序
				});
		}
	}

	pool.waitAllTasksDone();

	//----------------------------------------------------------------
	// 3) 如有后续操作，可直接用 cluster_results
	//----------------------------------------------------------------
	// anchors->clear(); …                // 原来的清理逻辑保持不变


	//----------------------------------------------------------------
	// 3) 释放原始 3D 数据以节省内存
	//----------------------------------------------------------------
	anchors->clear();
	anchors->shrink_to_fit();

}


