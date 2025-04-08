
#include "index.h"

FM_Index::FM_Index() {

}

bool FM_Index::newScan(const FilePath& fasta_path, const FilePath& output_path, uint_t thread) {
	std::map<uint64_t, NewScan::word_stats> wordFreq;
	uint64_t totChar;

	FilePath occ_path = output_path;
	occ_path += ".";
	occ_path += EXTOCC;

	FilePath dict_path = output_path;
	dict_path += ".";
	dict_path += EXTDICT;

	FilePath parse_path = output_path;
	parse_path += ".";
	parse_path += EXTPARSE;

	if (std::filesystem::exists(occ_path) && std::filesystem::exists(dict_path) && std::filesystem::exists(parse_path)) {
		return true;
	}


	NewScan::Args args;
	args.inputFileName = fasta_path.string();
	args.outputFileName = output_path.string();
	args.w = WINDOW_SIZE;
	args.p = STOP_MODULUS;
	args.th = thread;
	args.is_fasta = true;
	args.SAinfo = true;
	args.verbose = 0;
	args.compress = false;


	totChar = NewScan::mt_process_file_fasta(args, wordFreq);
	uint64_t totDWord = wordFreq.size();

	std::vector<const std::string*> dictArray;
	dictArray.reserve(totDWord);

	uint64_t sumLen = 0;
	uint64_t totWord = 0;
	for (auto& x : wordFreq) {
		sumLen += x.second.str.size();
		totWord += x.second.occ;
		dictArray.push_back(&x.second.str);
	}
	sort(dictArray.begin(), dictArray.end(), NewScan::pstringCompare);

	NewScan::writeDictOcc(args, wordFreq, dictArray);
	dictArray.clear(); // reclaim memory

	NewScan::remapParse(args, wordFreq);

	return true;
}

BWTParse::sa_index_t* FM_Index::compute_SA(uint32_t* Text, long n, long k)
{
	BWTParse::sa_index_t* SA = (BWTParse::sa_index_t*)malloc(n * sizeof(*SA));
	if (SA == NULL) die("malloc failed  (SA)");
	// printf("Computing SA of size %ld over an alphabet of size %ld\n", n, k);
	int depth = sacak_int(Text, SA, n, k);
	if (depth < 0) {
		die("Error computing the SA");
	}

	return SA;
}

bool FM_Index::bwtParse(const FilePath& fasta_path, const FilePath& output_path, uint_t thread) {


	FilePath last_path = output_path;
	last_path += ".";
	last_path += EXTBWLST;

	FilePath sai_path = output_path;
	sai_path += ".";
	sai_path += "bwsai";

	FilePath ilist_path = output_path;
	ilist_path += ".";
	ilist_path += "ilist";

	if (std::filesystem::exists(last_path) || std::filesystem::exists(sai_path) || std::filesystem::exists(ilist_path)) {
		return true;
	}

	uint32_t* Text; // array of parsing symbols
	long n;         // length of Text[] (not including final 0 symbol)
	size_t s;
	BWTParse::Args args;
	args.SAinfo = true;
	args.th = thread;

	std::string filename_str = output_path.string();
	char* filename_char = new char[filename_str.length() + 1];  // Create a modifiable char array
	std::strcpy(filename_char, filename_str.c_str()); // Copy the string into the char array
	args.filename = filename_char;  // Assign to args.filename, which expects char*

	// Pass the address of n to read_parse
	Text = BWTParse::read_parse(filename_char, &n);

	long k = 0;
	for (long i = 0;i < n;i++) {
		if (Text[i] > k) k = Text[i];
	}

	BWTParse::sa_index_t* SA = compute_SA(Text, n + 1, k + 1);


	// load last file 
	uint8_t* last = load_last(&args, n);
	// load sa info file, if requested
	uint8_t* sa_info = load_sa_info(&args, n);
	FILE* lastout = open_aux_file(args.filename, EXTBWLST, "wb");
	FILE* sa_out = open_sa_out(&args);

	BWTParse::sa_index_t* BWTsa = SA;
	BWTsa[0] = Text[n - 1];
	if (fputc(last[n - 2], lastout) == EOF) die("bwlast output 1");
	if (args.SAinfo) get_and_write_myint(sa_info, n, n - 1, sa_out);

	for (long i = 1;i <= n;i++) {
		if (SA[i] == 0) {
			assert(i == 1);  // Text[0]=$abc... is the second lex word 
			BWTsa[i] = 0;   // eos in BWT, there is no phrase in D corresponding to this symbol so we write dummy values
			if (fputc(0, lastout) == EOF) die("bwlast output 2"); // dummy char 
			if (args.SAinfo) write_myint(0, sa_out); // dummy end of word position, it is never used an 0 does not appear elsewhere in sa_out
		}
		else {
			if (SA[i] == 1) {
				// BWT[i] = Text[0] = $abcd... = first word in the parsing where $ now plays the role of the EOS in the original text  
				if (fputc(last[n - 1], lastout) == EOF) die("bwlast output 3");
			}
			else { if (fputc(last[SA[i] - 2], lastout) == EOF) die("bwlast output 4"); }
			if (args.SAinfo) get_and_write_myint(sa_info, n, SA[i] - 1, sa_out); // ending position of BWT symbol in original text
			BWTsa[i] = Text[SA[i] - 1];
		}
	}
	if (fclose(lastout) != 0) die("bwlast close");

	free(last);
	if (args.SAinfo) {
		if (fclose(sa_out) != 0) die("sa_out close");
		free(sa_info);
	}



	// --- copy BWT to text array (symbol by symbol since sizeof could be different)
	uint32_t* BWT = Text;
	for (long i = 0;i <= n;i++)
		BWT[i] = BWTsa[i];

	// read # of occ of each char from file .occ
	uint32_t* occ = (uint32_t*)malloc((k + 1) * sizeof(*occ)); // extra space for the only occ of 0
	if (occ == NULL) die("malloc failed (OCC)");
	FILE* occin = open_aux_file(args.filename, "occ", "rb");
	s = fread(occ + 1, sizeof(*occ), k, occin);
	if (s != k) die("not enough occ data!");
	occ[0] = 1; // we know there is somewhere a 0 BWT entry 
	fclose(occin);
	// create F vector
	uint32_t* F = (uint32_t*)malloc((k + 1) * sizeof(*F));
	if (F == NULL) die("malloc failed (F)");
	// init F[] using occ[]
	F[0] = 0;
	for (int i = 1;i <= k;i++)
		F[i] = F[i - 1] + occ[i - 1];
	assert(F[k] + occ[k] == n + 1);
	// ----- compute inverse list overwriting SA
	uint32_t* IList = (uint32_t*)SA;
	for (long i = 0;i <= n;i++) {
		IList[F[BWT[i]]++] = i;
		occ[BWT[i]]--;
	}
	// ---check
	assert(IList[0] == 1); // EOF is in BWT[1] since P[0] = $xxx is the smallest word and appears once
	assert(BWT[IList[0]] == 0);
	for (long i = 0;i <= k;i++)
		assert(occ[i] == 0);
	// ---save Ilist   
	FILE* ilist = open_aux_file(args.filename, EXTILIST, "wb");
	s = fwrite(IList, sizeof(*IList), n + 1, ilist);
	if (s != n + 1) die("Ilist write");
	fclose(ilist);
	// dpfBWTeallocate
	free(F);
	free(occ);
	free(SA);
	free(Text);

	return true;
}

bool FM_Index::pfBWT(const FilePath& fasta_path, const FilePath& output_path, uint_t thread) {
	FilePath sa_path = output_path;
	sa_path += ".sa";


	FilePath bwt_path = output_path;
	bwt_path += ".bwt";

	if (std::filesystem::exists(sa_path) && std::filesystem::exists(bwt_path)) {
		return true;
	}

	PfBWT::Args args;

	std::string filename_str = output_path.string();
	char* filename_char = new char[filename_str.length() + 1];  // Create a modifiable char array
	std::strcpy(filename_char, filename_str.c_str()); // Copy the string into the char array
	args.filename = filename_char;  // Assign to args.filename, which expects char*

	//args.sampledSA |= START_RUN;
	//args.sampledSA |= END_RUN;
	args.SA = true;

	args.w = WINDOW_SIZE;
	args.th = thread;

	FILE* g = open_aux_file(args.filename, EXTDICT, "rb");
	fseek(g, 0, SEEK_END);

	long dsize = ftell(g);
	if (dsize < 0) die("ftell dictionary");
	if (dsize <= 1 + args.w) die("invalid dictionary file");

	uint8_t* d = new uint8_t[dsize];
	rewind(g);
	long e = fread(d, 1, dsize, g);
	if (e != dsize) die("fread");
	fclose(g);

	// read occ file
	g = open_aux_file(args.filename, EXTOCC, "rb");
	fseek(g, 0, SEEK_END);
	e = ftell(g);
	if (e < 0) die("ftell occ file");
	if (e % 4 != 0) die("invalid occ file");
	int dwords = e / 4;
	uint32_t* occ = new uint32_t[dwords + 1];  // dwords+1 since istart overwrites occ
	rewind(g);
	e = fread(occ, 4, dwords, g);
	if (e != dwords) die("fread 2");
	fclose(g);

	// read ilist file 
	g = open_aux_file(args.filename, EXTILIST, "rb");
	fseek(g, 0, SEEK_END);
	e = ftell(g);
	if (e < 0) die("ftell ilist file");
	if (e % 4 != 0) die("invalid ilist file");
	long psize = e / 4;
	if (psize > 0xFFFFFFFEL) die("More than 2^32 -2 words in the parsing");
	uint32_t* ilist = new uint32_t[psize];
	rewind(g);
	e = fread(ilist, 4, psize, g);
	if (e != psize) die("fread 3");
	fclose(g);
	assert(ilist[0] == 1); // EOF is in PBWT[1] 

	g = open_aux_file(args.filename, EXTBWLST, "rb");
	uint8_t* bwlast = new uint8_t[psize];
	e = fread(bwlast, 1, psize, g);
	if (e != psize) die("fread 4");
	fclose(g);

	// convert occ entries into starting positions inside ilist
 // ilist also contains the position of EOF but we don't care about it since it is not in dict 
	uint32_t last = 1; // starting position in ilist of the smallest dictionary word  
	for (long i = 0;i < dwords;i++) {
		uint32_t tmp = occ[i];
		occ[i] = last;
		last += tmp;
	}
	assert(last == psize);
	occ[dwords] = psize;
	// extra check: the smallest dictionary word is d0 =$.... that occurs once
	assert(occ[1] == occ[0] + 1);


	PfBWT::bwt(args, d, dsize, ilist, bwlast, psize, occ, dwords); // version not using threads

	
	delete[] bwlast;
	delete[] ilist;
	delete[] occ;
	delete[] d;

	return true;
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


bool FM_Index::buildIndexUsingBigBWT(const FilePath& fasta_path, const FilePath& output_path, uint_t thread) {
	// newscanNT.x file -w 10 -p 100 -t 32 -s
	newScan(fasta_path, output_path, thread);

	// parseBWT input -s -t 32 
	bwtParse(fasta_path, output_path, thread);

	// pfbwt -w 10 input -s -e 
	pfBWT(fasta_path, output_path, thread);

	spdlog::info("bigBWT finished.");

	
	//std::cout << "Max value in SA: " << max_value << std::endl;


	//std::unordered_set<char> dna{ 'A', 'C', 'G', 'T', 'N' };
	//std::unordered_set<char> seen;

	//for (size_t i = 0; i < bwt.size(); ++i)
	//{
	//	char c = bwt[i];
	//	if (dna.find(c) == dna.end() && seen.find(c) == seen.end())
	//	{
	//		seen.insert(c);
	//		std::cout << "Non-DNA character found: '"
	//			<< ((c >= 32 && c <= 126) ? c : '?')  // 打印可见字符或 '?'
	//			<< "' (ASCII: " << static_cast<int>(static_cast<unsigned char>(c))
	//			<< ") at position " << i << "\n";
	//	}
	//}

	//init_bigbwt(output_path.string(), thread);

	//spdlog::info("r-index finished.");
	//save_to_file(output_path.string());
	//spdlog::info("r-index save to " + output_path.string());
	return true;
}

std::vector<uint32_t> read_parse_file(const std::string& filename) {
	std::ifstream file(filename, std::ios::binary);
	if (!file) {
		throw std::runtime_error("Cannot open parse file: " + filename);
	}

	// 获取文件大小
	file.seekg(0, std::ios::end);
	std::streamsize size = file.tellg();
	file.seekg(0);

	// 校验大小是否为4字节整数的倍数
	if (size % sizeof(uint32_t) != 0) {
		throw std::runtime_error("Invalid parse file size: not multiple of 4");
	}

	size_t num_entries = size / sizeof(uint32_t);
	std::vector<uint32_t> parse(num_entries);

	// 读取内容
	file.read(reinterpret_cast<char*>(parse.data()), size);
	if (!file) {
		throw std::runtime_error("Error reading parse file: " + filename);
	}

	return parse;
}

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
	//}
	// output_path = "/mnt/d/code/Another_project/Big-BWT/yeast.fasta";
	FilePath bwt_path = output_path;
	// FilePath bwt_path = "/mnt/d/code/Another_project/Big-BWT/yeast.fasta"
	bwt_path += ".bwt";
	std::ifstream bwt_file(bwt_path, std::ios::binary | std::ios::ate);
	/*	std::string bwt((std::istreambuf_iterator<char>(bwt_file)),
			std::istreambuf_iterator<char>())*/;

	std::streamsize size = bwt_file.tellg();   // 获取文件大小
	bwt_file.seekg(0, std::ios::beg);          // 回到文件开头

	std::string bwt(size, '\0');        // 分配足够的空间
	bwt_file.read(&bwt[0], size);           // 一次性读入
	bwt.erase(0, 1);  // 从索引0开始，删除1个字符



	FilePath sa_path = output_path;
	sa_path += ".sa";
	////		//// 遍历sa，打印出最大值
	////uint64_t max_value = 0;
	////for (size_t i = 0; i < sa.size(); ++i)
	////{
	////	if (sa[i] > max_value)
	////	{
	////		max_value = sa[i];
	////	}
	////}
	std::vector<uint64_t> sa = read_sa(sa_path);

	//size_t eof_pos = bwt.find('\0');
	//bwt.erase(bwt.begin() + eof_pos);

	//// sa.insert(sa.begin() + eof_pos, sa.size());


	std::string text = fasta_manager.concatRecords();

	bool passed = true;
	size_t n = text.size();

	for (size_t i = 0; i < bwt.size(); ++i) {
		uint64_t sai = sa[i];
		if (i < 346652195 && i > 346652175) {
			std::cout << text.substr(sai-1, 10) << std::endl;
		}
		if (sai == n - 1) {
			std::cout << "final char" << std::endl;
		}
		// cout << text.substr(sai, 10) << endl;

		// sai 应该在 0..n 范围内（n == text.size()）
		if (sai > n) {
			std::cerr << "❌ SA[" << i << "] = " << sai << " exceeds text length " << n << std::endl;
			passed = false;
			break;
		}
		// cout << text.substr(sai - 1, 10) << endl;
		// BWT[i] 应该是 text[(SA[i]-1+n)%n] 的字符
		char expected_char = text[(sai - 1 + n) % n];
		char expected_char2 = text[(sai + n) % n];

		if (bwt[i] != expected_char) {
			std::cerr << "❌ Mismatch at i=" << i
				<< ": BWT[i]=" << (int)(unsigned char)bwt[i]
				<< " vs expected=" << (int)(unsigned char)expected_char
					<< " (SA[i]=" << sai << ")\n";
				passed = false;
				break;
		}
	}

	if (passed) {
		std::cout << "✅ SA and BWT match perfectly!\n";
	}
	else {
		std::cerr << "❌ SA and BWT verification failed!\n";
	}

	

	
	return index_path;

}