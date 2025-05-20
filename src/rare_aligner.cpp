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
#include <chrono>

AnchorPtrListVec PairRareAligner::findQueryFileAnchor(
	const std::string prefix,
	const FilePath fasta_path,
	SearchMode         search_mode)
{
	namespace ch = std::chrono;
	FilePath result_dir_path = work_dir / RESULT_DIR;
	std::filesystem::create_directories(result_dir_path);

	FilePath anchor_file = result_dir_path /
		(prefix + "_" + SearchModeToString(search_mode) + "." + ANCHOR_EXTENSION);

	/* ---------- 若已存在结果文件，直接加载 ---------- */
	if (std::filesystem::exists(anchor_file)) {
		AnchorPtrListVec result;
		loadAnchors(anchor_file.string(), result);
		return result;
	}

	/* ---------- 读取 FASTA 并分片 ---------- */
	FastaManager query_fasta_manager(fasta_path, getFaiIndexPath(fasta_path));
	RegionVec chunks = query_fasta_manager.preAllocateChunks(chunk_size, overlap_size);

	/* ---------- ① 计时：搜索 Anchor ---------- */
	auto t_search0 = ch::steady_clock::now();

	ThreadPool pool(thread_num);
	std::vector<std::future<AnchorPtrListVec>> futures;
	futures.reserve(chunks.size());

	for (const auto& ck : chunks) {
		std::string seq = query_fasta_manager.getSubSequence(ck.chr_name, ck.start, ck.length);

		futures.emplace_back(
			pool.enqueue(
				[this, ck, seq, search_mode]() -> AnchorPtrListVec {
					return ref_index.findAnchors(
						ck.chr_name, seq, search_mode,
						Strand::FORWARD,
						ck.start,
						min_anchor_length,
						max_anchor_frequency);
				}));
	}

	pool.waitAllTasksDone();

	auto t_search1 = ch::steady_clock::now();
	double search_ms = ch::duration<double, std::milli>(t_search1 - t_search0).count();
	spdlog::info("[findQueryFileAnchor] search  = {:.3f} ms", search_ms);
	/* ---------- ② 计时：合并 ---------- */
	auto t_merge0 = ch::steady_clock::now();
	std::vector<AnchorPtrListVec> parts;
	parts.reserve(futures.size());
	size_t total_lists = 0;
	for (auto& fut : futures) {
		AnchorPtrListVec part = fut.get();
		total_lists += part.size();
		parts.emplace_back(std::move(part));
	}

	AnchorPtrListVec all_results;
	all_results.reserve(total_lists);
	for (auto& part : parts)
		for (auto& lst : part)
			all_results.emplace_back(std::move(lst));

	auto t_merge1 = ch::steady_clock::now();
	double merge_ms = ch::duration<double, std::milli>(t_merge1 - t_merge0).count();
	spdlog::info("[findQueryFileAnchor] merge   = {:.3f} ms", merge_ms);
	/* ---------- ③ 计时：保存 ---------- */
	auto t_save0 = ch::steady_clock::now();
	// saveAnchors(anchor_file.string(), all_results);
	auto t_save1 = ch::steady_clock::now();
	double save_ms = ch::duration<double, std::milli>(t_save1 - t_save0).count();

	/* ---------- 打印统计 ---------- */
	spdlog::info("[findQueryFileAnchor] save    = {:.3f} ms", save_ms);

	return all_results;          // NRVO / move-elided
}

// TODO 这个函数是用于筛选初步得到的锚点
// 现在已有的代码是把MUM和MEM分开
// 优先处理MUM，MEM是对MUM比对结果的补充
// 最后返回到结果应该是AnchorPtrListVec。(我们现在做的内容是比对结果是一个ref区域对应一个query区域，AnchorPtr够用了，但后续可能会扩展到多对多，所以预先留好这个接口)
void PairRareAligner::filterAnchors(AnchorPtrListVec& anchors)
{
	// 清空上次的结果
	unique_anchors.clear();
	repeat_anchors.clear();

	// 遍历所有 AnchorPtrList
	for (auto& alist : anchors) {
		if (alist.size() == 1) {
			// 只有一个 Anchor，那就是 MUM（唯一匹配）
			unique_anchors.push_back(std::move(alist));
		}
		else {
			// 多于一个的，都归为重复匹配
			repeat_anchors.push_back(std::move(alist));
		}
	}
	return;
}
