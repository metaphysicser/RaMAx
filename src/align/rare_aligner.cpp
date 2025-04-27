#include "align.h"


PairRareAligner::PairRareAligner(const FilePath work_dir, const uint_t thread_num) {
	this->work_dir = work_dir;
	this->index_dir = work_dir / INDEX_DIR;
	this->thread_num = thread_num;

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

FilePath PairRareAligner::buildIndex(const std::string prefix, const FilePath fasta_path) {
	FilePath ref_index_path = index_dir / prefix;

	if (!std::filesystem::exists(ref_index_path)) {
		std::filesystem::create_directories(ref_index_path);
	}

	ref_fasta_manager = FastaManager(fasta_path, getFaiIndexPath(fasta_path));

	// index_dir的路径加上prefix的前缀加上fasta_manager.fasta_path_的扩展名
	FilePath output_path = ref_index_path / (prefix + std::filesystem::path(fasta_path).extension().string());

	FM_Index fm_index(prefix, &ref_fasta_manager);

	FilePath idx_file_path = ref_index_path / (prefix + "." + FMINDEX_EXTESION);

	spdlog::info("Indexing with prefix: {}, index path: {}", prefix, ref_index_path.string());

	

	if (!std::filesystem::exists(idx_file_path)) {
		fm_index.buildIndex(output_path, false, thread_num);
		fm_index.saveToFile(idx_file_path.string());
	}
	else {
		fm_index.loadFromFile(idx_file_path.string());
	}
	std::string query = "NCACAGCGAGCTATCGATCGTAGCTAGCTAGCTAGCTCGTAGCTAACACTGTGTGTACTACGACTAGCTACAACACAGCGAGCTATCGATCGTAGCTAGCTAGCTAG";
	AnchorPtrVec p = fm_index.findAnchors("qeury", query, FORWARD, 0, 20, 50);
	spdlog::info("Indexing finished, index path: {}", ref_index_path.string());

	return ref_index_path;
}