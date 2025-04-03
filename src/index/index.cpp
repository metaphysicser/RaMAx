
#include "index.h"

R_Index::R_Index() {

}

bool R_Index::newScan(const FilePath& fasta_path, const FilePath& output_path, uint_t thread) {
	map<uint64_t, NewScan::word_stats> wordFreq;
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

	if (std::filesystem::exists(occ_path) && std::filesystem::exists(dict_path)) {
		return true;
	}
	else {
		totChar = NewScan::mt_process_file_fasta(args, wordFreq);
		uint64_t totDWord = wordFreq.size();

		std::vector<const string*> dictArray;
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
	}

	
	if (std::filesystem::exists(parse_path)) {
		return true;
	}
	else {
		NewScan::remapParse(args, wordFreq);
	}
	

	return true;
}

BWTParse::sa_index_t* R_Index::compute_SA(uint32_t* Text, long n, long k)
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

bool R_Index::bwtParse(const FilePath& fasta_path, const FilePath& output_path, uint_t thread) {
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
	// deallocate
	free(F);
	free(occ);
	free(SA);
	free(Text);

	return true;
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

	R_Index r_index;
	// newscanNT.x file -w 10 -p 100 -t 32 -s
	r_index.newScan(fasta_manager.fasta_path_, output_path, thread_num);

	// parseBWT input -s -t 32 
	r_index.bwtParse(fasta_manager.fasta_path_, output_path, thread_num);

	// pfbwt -w 10 input -s -e 

	// R_Index r_index;

	
	return index_path;

}