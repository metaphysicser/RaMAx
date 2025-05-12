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
void PairRareAligner::alignQueryFile(const std::string prefix,
	const FilePath fasta_path, SearchMode search_mode)
{
	// 1) 打开 FASTA 并分片
	FastaManager query_fasta_manager(fasta_path,
		getFaiIndexPath(fasta_path));
	RegionVec chunks = query_fasta_manager.preAllocateChunks(chunk_size,
		overlap_size);

	// 2) 启动线程池并提交任务
	ThreadPool pool(thread_num);
	std::vector<std::future<AnchorPtrListVec>> futures;
	futures.reserve(chunks.size());

	for (auto& chunk : chunks) {
		std::string seq = query_fasta_manager.getSubSequence(
			chunk.chr_name, chunk.start, chunk.length);

		futures.emplace_back(
			pool.enqueue(
				[this, chunk, seq, search_mode]() -> AnchorPtrListVec {
					return ref_index.findAnchors(
						chunk.chr_name,      // query_chr
						seq,                 // query 序列
						search_mode,
						Strand::FORWARD,     // strand
						chunk.start,         // query_offset
						min_anchor_length,   // 最小锚点长度
						max_anchor_frequency // 最大频次过滤
					);
				}
			)
		);
	}

	// 3) 等待所有任务完成
	pool.waitAllTasksDone();

	// 4) 合并所有 future 的结果到一个 AnchorPtrListVec
	AnchorPtrListVec all_results;
	// （可视情况预估总组数并 reserve，这里简单用 chunks.size() 预留最外层容量）
	all_results.reserve(chunks.size());

	for (auto& fut : futures) {
		AnchorPtrListVec part = fut.get();
		all_results.insert(
			all_results.end(),
			std::make_move_iterator(part.begin()),
			std::make_move_iterator(part.end())
		);
	}

	return;

	// 5) 后续处理——例如写出文件或进一步分析
	// writeAnchors(prefix + ".anchors", all_results);
}

