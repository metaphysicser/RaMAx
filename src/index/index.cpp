
#include "index.h"


IndexManager::IndexManager(const FilePath work_dir, const uint_t thread_num) {
	this->work_dir = work_dir;
	this->index_dir = work_dir / INDEX_DIR;
	this->thread_num = thread_num;
}

FilePath IndexManager::buildIndex(const std::string prefix, FastaManager& fasta_manager, const IndexType index_type) {
	FilePath index_path = index_dir / prefix;

	if (!std::filesystem::exists(index_path)) {
		std::filesystem::create_directories(index_path);
	}

	// index_dir的路径加上prefix的前缀加上fasta_manager.fasta_path_的扩展名
	FilePath output_path = index_path / (prefix + std::filesystem::path(fasta_manager.fasta_path_).extension().string());

	FM_Index fm_index;
	spdlog::info("Indexing with prefix: {}, index path: {}", prefix, index_path.string());
	
	fm_index.buildIndexUsingBigBWT(fasta_manager.fasta_path_, output_path, thread_num);

	spdlog::info("Indexing finished, index path: {}", index_path.string());

	//FilePath parse_path = output_path;
	//parse_path += ".parse";

	//std::vector<uint32_t> parse = read_parse_file(parse_path);

	//const std::size_t subproblem_count(0);
	//const std::size_t max_context(0);

	//// 使用parse数据构建后缀数组
	//// CaPS_SA::Suffix_Array<uint32_t> sa(parse, parse.size(), 0, 0);

	//// 构建后缀数组
	//// sa.construct();

	//// const uint32_t* SA = sa.SA();
	//std::string text = "AAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGG$";

	//cout << text.size() << endl;
	///*CaPS_SA::Suffix_Array<uint32_t> suf_arr(text.c_str(), text.length(), subproblem_count, max_context);*/
	//CaPS_SA::Suffix_Array<uint32_t> suf_arr(text.c_str(), text.size(), subproblem_count, max_context);
	//suf_arr.construct();
	//const uint32_t* SA = suf_arr.SA();
	//const uint32_t* LCP = suf_arr.LCP();

	//// 打印sa
	//for (size_t i = 0; i < text.size(); ++i) {
	//	std::cout << SA[i] << " \n";
	//}
	//cout << "-------_______________-----------" << endl;
	//// 打印LCP
	//for (size_t i = 0; i < text.size(); ++i) {
	//	std::cout << LCP[i] << " \n";
	////}
	//// output_path = "/mnt/d/code/Another_project/Big-BWT/yeast.fasta";
	//FilePath bwt_path = output_path;
	//// FilePath bwt_path = "/mnt/d/code/Another_project/Big-BWT/yeast.fasta"
	//bwt_path += ".bwt";
	//std::ifstream bwt_file(bwt_path, std::ios::binary | std::ios::ate);
	///*	std::string bwt((std::istreambuf_iterator<char>(bwt_file)),
	//		std::istreambuf_iterator<char>())*/;

	//std::streamsize size = bwt_file.tellg();   // 获取文件大小
	//bwt_file.seekg(0, std::ios::beg);          // 回到文件开头

	//std::string bwt(size, '\0');        // 分配足够的空间
	//bwt_file.read(&bwt[0], size);           // 一次性读入
	//bwt.erase(0, 1);  // 从索引0开始，删除1个字符



	//FilePath sa_path = output_path;
	//sa_path += ".sa";
	//////		//// 遍历sa，打印出最大值
	//////uint64_t max_value = 0;
	//////for (size_t i = 0; i < sa.size(); ++i)
	//////{
	//////	if (sa[i] > max_value)
	//////	{
	//////		max_value = sa[i];
	//////	}
	//////}
	//std::vector<uint64_t> sa = read_sa(sa_path);

	////size_t eof_pos = bwt.find('\0');
	////bwt.erase(bwt.begin() + eof_pos);

	////// sa.insert(sa.begin() + eof_pos, sa.size());


	//std::string text = fasta_manager.concatRecords();

	//bool passed = true;
	//size_t n = text.size();

	//for (size_t i = 0; i < bwt.size(); ++i) {
	//	uint64_t sai = sa[i];
	//	if (i < 346652195 && i > 346652175) {
	//		std::cout << text.substr(sai-1, 10) << std::endl;
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

	

	
	return index_path;

}