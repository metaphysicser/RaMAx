
#include "index.h"


IndexManager::IndexManager(const FilePath work_dir, const uint_t thread_num) {
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

FilePath IndexManager::buildIndex(const std::string prefix, FastaManager& fasta_manager) {
	FilePath index_path = index_dir / prefix;

	if (!std::filesystem::exists(index_path)) {
		std::filesystem::create_directories(index_path);
	}

	// index_dir的路径加上prefix的前缀加上fasta_manager.fasta_path_的扩展名
	FilePath output_path = index_path / (prefix + std::filesystem::path(fasta_manager.fasta_path_).extension().string());

	FM_Index fm_index(&fasta_manager);
	spdlog::info("Indexing with prefix: {}, index path: {}", prefix, index_path.string());
	
	fm_index.buildIndex(fasta_manager, output_path, false, thread_num);

	spdlog::info("Indexing finished, index path: {}", index_path.string());


	std::cout << fasta_manager.getSubConcatSequence(0, 800000) << std::endl;
	//FilePath sa_path = output_path;
	//sa_path += ".sa";
	//std::vector<uint64_t> sa = read_sa(sa_path);
	//size_t sa_size = sa.size();


	//FilePath bwt_path = output_path;
	//bwt_path += ".bwt";
	//std::ifstream bwt_file(bwt_path, std::ios::binary | std::ios::ate);


	//std::streamsize size = bwt_file.tellg();   // 获取文件大小
	//bwt_file.seekg(0, std::ios::beg);          // 回到文件开头

	//std::string bwt(size, '\0');        // 分配足够的空间
	//bwt_file.read(&bwt[0], size);           // 一次性读入
	//bwt.erase(0, 1);  // 从索引0开始，删除1个字符

	//// 打印bwt中除了ATCG的值
	//for (size_t i = 0; i < bwt.size(); ++i) {
	//	if (bwt[i] != 'A' && bwt[i] != 'G' && bwt[i] != 'C' && bwt[i] != 'T') {
	//		std::cout << "bwt[" << i << "] = " << int(bwt[i]) << " " << bwt[i] << std::endl;
	//	}
	//}



	//std::string text = fasta_manager.concatRecords();

	//// 打印text中除了ATCG的值
	////for (size_t i = 0; i < text.size(); ++i) {
	////	if (text[i] != 'A' && text[i] != 'C' && text[i] != 'G' && text[i] != 'T') {
	////		std::cout << "text[" << i << "] = " << text[i] << std::endl;
	////	}
	////}

	//bool passed = true;
	//size_t n = text.size();

	//for (size_t i = 0; i < bwt.size(); ++i) {
	//	uint64_t sai = sa[i];
	//	if (i < 346652195 && i > 346652175) {
	//		std::cout << text.substr(sai - 1, 10) << std::endl;
	//	}
	//	if (sai == n - 1) {
	//		std::cout << "final char" << std::endl;
	//	}
	//	// cout << text.substr(sai, 10) << endl;

	//	// sai 应该在 0..n 范围内（n == text.size()）
	//	if (sai > n) {
	//		std::cerr << "❌ SA[" << i << "] = " << sai << " exceeds text length " << n << std::endl;
	//		passed = false;
	//		break;
	//	}
	//	// cout << text.substr(sai - 1, 10) << endl;
	//	// BWT[i] 应该是 text[(SA[i]-1+n)%n] 的字符
	//	char expected_char = text[(sai - 1 + n) % n];
	//	char expected_char2 = text[(sai + n) % n];

	//	if (bwt[i] != expected_char) {
	//		std::cerr << "❌ Mismatch at i=" << i
	//			<< ": BWT[i]=" << (int)(unsigned char)bwt[i]
	//			<< " vs expected=" << (int)(unsigned char)expected_char
	//				<< " (SA[i]=" << sai << ")\n";
	//			passed = false;
	//			break;
	//	}
	//}

	//if (passed) {
	//	std::cout << "✅ SA and BWT match perfectly!\n";
	//}
	//else {
	//	std::cerr << "❌ SA and BWT verification failed!\n";
	//}
	//
	
	return index_path;

}