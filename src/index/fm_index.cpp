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

bool FM_Index::read_and_build_sampled_sa(const FilePath& sa_file_path)
{
	std::ifstream file(sa_file_path, std::ios::binary | std::ios::ate);
	if (!file) throw std::runtime_error("Cannot open SA file");

	std::streamsize size = file.tellg();
	file.seekg(0);

	if (size % 5 != 0) {
		throw std::runtime_error("SA file size not multiple of 5 bytes");
	}

	size_t count = size / 5;

	// Create a buffer to hold the raw SA data
	std::vector<char> buffer(size);
	file.read(buffer.data(), size);

	// Calculate the number of sampled suffixes based on the sampling rate of 32
	size_t sampled_count = (count + 31) / 32;

	// Resize sampled_sa to the correct number of elements based on sampled_count
	sampled_sa.resize(sampled_count);

	const uint8_t* ptr = reinterpret_cast<const uint8_t*>(buffer.data());

	size_t sampled_index = 0;

	for (size_t i = 0; i < count; i += 32) { // Step by 32 to sample every 32nd suffix
		ptr = reinterpret_cast<const uint8_t*>(buffer.data()) + (i * 5); // Move the pointer to the right position

		// Little-Endian interpretation:
		uint64_t val =
			(uint64_t(ptr[0])) |
			(uint64_t(ptr[1]) << 8) |
			(uint64_t(ptr[2]) << 16) |
			(uint64_t(ptr[3]) << 24) |
			(uint64_t(ptr[4]) << 32);

		// Ensure sampled_index does not exceed the size of sampled_sa
		
		sampled_sa[sampled_index] = val;
		sampled_index++;
		
	}

	return true;
}




bool FM_Index::buildIndexUsingBigBWT(const FilePath& fasta_path, const FilePath& output_path, uint_t thread) {
	// newscanNT.x file -w 10 -p 100 -t 32 -s
	newScan(fasta_path, output_path, thread);

	// parseBWT input -s -t 32 
	bwtParse(fasta_path, output_path, thread);

	// pfbwt -w 10 input -s -e 
	pfBWT(fasta_path, output_path, thread);

	spdlog::info("bigBWT finished.");

	FilePath sa_path = output_path;
	sa_path += ".sa";

	read_and_build_sampled_sa(sa_path);
	sdsl::util::bit_compress(sampled_sa);

	FilePath bwt_path = output_path;
	bwt_path += ".bwt";
	read_and_build_bwt(bwt_path);

	return true;
}

bool FM_Index::read_and_build_bwt(const FilePath& bwt_file_path) {
	// 打开文件，检查是否成功
	std::ifstream bwt_file(bwt_file_path, std::ios::binary | std::ios::ate);
	if (!bwt_file.is_open()) {
		throw std::runtime_error("Error: Cannot open file " + bwt_file_path);
	}

	// 获取文件大小，检查文件大小是否合理
	std::streamsize size = bwt_file.tellg();
	if (size <= 0) {
		throw std::runtime_error("Error: Invalid file size for " + bwt_file_path);
	}

	// 将文件指针重置到文件开头
	bwt_file.seekg(0, std::ios::beg);

	// 一次性将整个文件读入 std::string 中
	std::string bwt(size, '\0');
	if (!bwt_file.read(&bwt[0], size)) {
		throw std::runtime_error("Error: Failed to read file " + bwt_file_path);
	}

	// 关闭文件，以释放文件句柄
	bwt_file.close();

	// 根据需求删除索引 0 处的1个字符
	if (!bwt.empty()) {
		bwt.erase(0, 1);
	}
	else {
		throw std::runtime_error("Error: BWT string is empty after reading file " + bwt_file_path);
	}

	// 将 std::string 转换为 sdsl::int_vector<8>
	sdsl::int_vector<8> sdsl_bwt(bwt.size());
	for (size_t i = 0; i < bwt.size(); ++i) {
		sdsl_bwt[i] = static_cast<uint8_t>(bwt[i]);
	}

	// 及时释放原始 BWT 字符串占用的内存
	bwt.clear();
	bwt.shrink_to_fit();

	// 构造基于 wt_huff 的波列树（内部默认使用 sdsl::bit_vector）
	// this->wt_bwt 是 FM_Index 类中的成员变量，例如：
	// sdsl::wt_huff<sdsl::bit_vector> wt_bwt;
	sdsl::construct_im(this->wt_bwt, sdsl_bwt);
	//std::cout << "---- Access Operation Demo ----\n";
	//for (size_t i = 0; i < std::min<size_t>(10, sdsl_bwt.size()); ++i) {
	//	// 由于 wt 存储的是 uint8_t，这里转换为 char 进行输出
	//	std::cout << "wt.access(" << i << ") = " << static_cast<char>(wt_bwt[i]) << "\n";
	//}

	//// 演示 rank 操作
	//// 比如统计前 20 个字符中 'A' 出现的次数
	//uint8_t symbol = static_cast<uint8_t>('C');  // 可以换成其他字符如 'C', 'G', 'T'
	//size_t pos = 20;
	//size_t count = wt_bwt.rank(pos, symbol);
	//std::cout << "Rank of '" << static_cast<char>(symbol)
	//	<< "' in positions [0, " << pos << ") = " << count << "\n";

	//// 演示 select 操作
	//// 查找第 2 次出现 'A' 的位置，注意 select 索引从1开始计数
	//size_t kth = 2;
	//size_t pos_sel = wt_bwt.select(kth, symbol);
	//std::cout << "Select(" << kth << ", '" << static_cast<char>(symbol)
	//	<< "') returns position: " << pos_sel << "\n";

	//size_t raw_bwt_size = sdsl::size_in_bytes(sdsl_bwt);
	//std::cout << "Raw BWT size: " << raw_bwt_size << " bytes\n";

	//size_t wt_size = sdsl::size_in_bytes(wt_bwt);
	//std::cout << "Compressed WT (wt_huff) size: " << wt_size << " bytes\n";
	//
	//double compression_ratio = static_cast<double>(wt_size) / raw_bwt_size;
	//std::cout << "Compression Ratio: " << compression_ratio * 100.0 << " %\n";
	return true;
}

bool FM_Index::buildIndex(FastaManager& fasta_manager, FilePath output_path, bool fast_mode, uint_t thread) {

	if (fast_mode) {
		// TODO CaPS-SA
		return true;
	}
	else {
		
		buildIndexUsingBigBWT(fasta_manager.fasta_path_, output_path, thread);
	}
}