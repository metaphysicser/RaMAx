#include "index.h"


FM_Index::FM_Index(SpeciesName species_name, FastaManager* fasta_manager, uint_t sample_rate):species_name(species_name), fasta_manager(fasta_manager) {
	this->sample_rate = sample_rate;
	this->total_size = fasta_manager->getConcatSeqLength();
	if (total_size < sample_rate) {
		this->sample_rate = 1;
	}
}

bool FM_Index::buildIndexUsingCaPS(uint_t thread_count)
{
	this->total_size += 1;
	std::string T = fasta_manager->concatRecords();
	size_t n = T.size();
	if (n == 0) return false;

	if (n <= std::numeric_limits<uint32_t>::max())
	{
		return buildIndexUsingCaPSImpl<uint32_t>(T, thread_count);
	}
	else
	{
		return buildIndexUsingCaPSImpl<uint64_t>(T, thread_count);
	}
}

template<typename index_t>
bool FM_Index::buildIndexUsingCaPSImpl(const std::string& T, uint_t thread_count)
{
	using SA_t = index_t;

	size_t n = T.size();

	std::string value = std::to_string(thread_count);
	setenv("PARLAY_NUM_THREADS", value.c_str(), 1);  // 1 means overwrite if exists

	// 构造 CaPS 后缀数组
	CaPS_SA::Suffix_Array<SA_t> suf_arr(T.c_str(), static_cast<SA_t>(n), 0, 0);
	suf_arr.construct();
	// 拷贝 SA_ 指针内容到 vector，方便后续操作（线程安全）
	std::vector<SA_t> SA(n);
	std::copy(suf_arr.SA(), suf_arr.SA() + n, SA.begin());

	// 构建 BWT 向量
	sdsl::int_vector<8> BWT(n);

	ThreadPool pool(thread_count);
	size_t chunk_size = (n + thread_count - 1) / thread_count;

	for (uint_t t = 0; t < thread_count; ++t) {
		size_t start = t * chunk_size;
		size_t end = std::min(start + chunk_size, n);

		pool.enqueue(
			[&T, &SA, &BWT, n, start, end]() {
				for (size_t i = start; i < end; ++i) {
					size_t si = static_cast<size_t>(SA[i]);
					BWT[i] = static_cast<uint8_t>(T[(si + n - 1) % n]);
				}
			}
		);
	}
	pool.waitAllTasksDone();

	// 构建采样 SA
	size_t sampled_count = (n + sample_rate - 1) / sample_rate;
	sampled_sa.resize(sampled_count);
	for (size_t i = 0, idx = 0; i < n; i += sample_rate, ++idx) {
		sampled_sa[idx] = static_cast<uint64_t>(SA[i]);
	}
	sdsl::util::bit_compress(sampled_sa);

	// 构建波形树
	sdsl::construct_im(this->wt_bwt, BWT);

	return true;
}



bool FM_Index::buildIndexUsingDivfsort(
	uint_t thread_count)
{
	// 1) Get the full concatenated sequence
	this->total_size += 1;
	std::string T = fasta_manager->concatRecords();
	size_t n = T.size();
	if (n == 0) return false;

	//// 2) Prepare the input buffer T[0..n-1] for divsufsort
	//std::vector<sauchar_t> T(n);
	//std::memcpy(T.data(), concat_seq.data(), n);
	const auto* Tptr = reinterpret_cast<const sauchar_t*>(T.data());

	// 3) Allocate space for SA[0..n-1]
	std::vector<saidx_t> SA(n);

	// 4) Build the suffix array
	if (divsufsort(Tptr, SA.data(), static_cast<saidx_t>(n)) < 0) {
		spdlog::error("divsufsort() failed");
		return false;
	}

	// 5) Allocate BWT buffer
	sdsl::int_vector<8> BWT(n);

	// 6) Parallel BWT construction:
	//    BWT[i] = T[( SA[i] + n - 1 ) % n]
	ThreadPool pool(thread_count);
	size_t chunk_size = (n + thread_count - 1) / thread_count;

	for (uint_t t = 0; t < thread_count; ++t) {
		size_t start = t * chunk_size;
		size_t end = std::min(start + chunk_size, n);

		pool.enqueue(
			[&T, &SA, &BWT, n, start, end]() {
				for (size_t i = start; i < end; ++i) {
					saidx_t si = SA[i];
					BWT[i] = static_cast<uint8_t>(
						T[(static_cast<size_t>(si) + n - 1) % n]
						);

				}
			}
		);
	}

	pool.waitAllTasksDone();

	size_t sampled_count = (n + sample_rate - 1) / sample_rate;

	sampled_sa.resize(sampled_count);
	for (size_t i = 0, idx = 0; i < n; i += sample_rate, ++idx) {
		sampled_sa[idx] = SA[i];
	}
	//SA.clear();             // 把 size() 设为 0
	//SA.shrink_to_fit();
	sdsl::util::bit_compress(sampled_sa);

	sdsl::construct_im(this->wt_bwt, BWT);

	return true;
}

bool FM_Index::buildIndex(FilePath output_path, bool fast_mode, uint_t thread) {

	if (fast_mode) {
		// TODO CaPS-SA
		buildIndexUsingCaPS(thread);
	}
	else {
		if (isFileSmallerThan(fasta_manager->fasta_path_, 1024)) {
			buildIndexUsingDivfsort(thread);
		}
		else {
			buildIndexUsingBigBWT(output_path, thread);
		}
		
	}
	
	size_t cumulative = 0;
	char2idx.fill(0xFF);
	count_array.fill(0);

	if (fasta_manager->has_n_in_fasta) {
		uint_t count = 0;
		// alpha_set 中包含 N，假设 alpha_set 已经按照字典序排序
		for (const auto& c : alpha_set) {
			// 计算字符 c 在完整 BWT 中的总出现次数
			size_t occ = wt_bwt.rank(wt_bwt.size(), c);
			// 对于当前字符 c，其C值即为前面所有字符出现次数的累积值
			count_array[static_cast<uint8_t>(c)] = cumulative;
			// 更新累积值
			cumulative += occ;
			char2idx[static_cast<uint8_t>(c)] = count;
			count++;
		}
	}
	else {
		// alpha_set_without_N 不包含 N，同样要求字母集合已排序
		uint_t count = 0;
		for (const auto& c : alpha_set_without_N) {
			size_t occ = wt_bwt.rank(wt_bwt.size(), c);
			count_array[static_cast<uint8_t>(c)] = cumulative;
			cumulative += occ;
			char2idx[static_cast<uint8_t>(c)] = count;
			count++;
		}
	}


	return true;
}

uint_t FM_Index::LF(uint_t pos) const {
	uint8_t c = wt_bwt[pos];
	return count_array[c] + wt_bwt.rank(pos, c);
}

uint_t FM_Index::getSA(size_t pos) const {
	uint_t steps = 0;
	size_t cur = pos;
	
		   /* 向上跳直到落在采样位置 (每 32 个一采样) */
	while ((cur) % sample_rate != 0) {
		cur = LF(cur);
		++steps;	
	}
	
		    /* 采样值 + 步数 (= text_length 时取模) */
	size_t sa_sample = sampled_sa[(cur) / sample_rate];
	size_t sa_val = sa_sample + steps;
	if (sa_val >= total_size) sa_val -= total_size;
	return sa_val;
	
}

SAInterval FM_Index::backwardExtend(const SAInterval& I, char c) {
	uint8_t ch = static_cast<uint8_t>(c);
	uint_t new_l = count_array[ch] + wt_bwt.rank(I.l, ch);
	uint_t new_r = count_array[ch] + wt_bwt.rank(I.r, ch);
	return { new_l, new_r };
}

AnchorPtrListVec FM_Index::findAnchors(ChrName query_chr, std::string query, SearchMode search_mode, Strand strand, uint_t query_offset, uint_t min_anchor_length, uint_t max_anchor_frequency)
{
	if (strand == FORWARD) {
		std::reverse(query.begin(), query.end());
	}
	else {
		for (char& ch : query) {
			ch = BASE_COMPLEMENT[static_cast<unsigned char>(ch)];
		}
	}
	AnchorPtrListVec anchor_ptr_list_vec;

	uint_t total_length = 0;
	uint_t query_length = query.length();
	uint_t last_pos = 0;
	while (total_length < query_length) {
		RegionVec region_vec;		
		uint_t match_length = findSubSeqAnchors(query.c_str() + total_length, query_length - total_length,
				region_vec, min_anchor_length, max_anchor_frequency);
		
		if (region_vec.size() > 0) {
			Region query_region(query_chr, total_length + query_offset, match_length);
			AnchorPtrList anchor_ptr_list;

			uint_t ref_end_pos = region_vec[0].start + match_length;
			if (search_mode == ACCURATE_SEARCH && ref_end_pos == last_pos) {
				total_length += match_length;
			}
			else {
				for (uint_t i = 0; i < region_vec.size(); i++) {
					Match match(region_vec[i], query_region, strand);
					Score_t score = caculateMatchScore(query.c_str() + total_length, match_length);
					Cigar_t cigar;
					cigar.push_back(cigarToInt('=', match_length));
					AnchorPtr p = std::make_shared<Anchor>(match, match_length, cigar, score);
					anchor_ptr_list.push_back(p);
				}
				anchor_ptr_list_vec.push_back(anchor_ptr_list);

				if (search_mode == FAST_SEARCH) {
					total_length += match_length;
				}
				else if (search_mode == ACCURATE_SEARCH) {
					total_length += min_anchor_length;
				}
			}
			last_pos = ref_end_pos;
		}
		else {
			total_length += 1;
		}
		
	}
	return anchor_ptr_list_vec;
}


uint_t FM_Index::findSubSeqAnchors(const char* query, uint_t query_length, RegionVec& region_vec, uint_t min_anchor_length, uint_t max_anchor_frequency)
{
	uint_t match_length = 0;
	SAInterval I = { 0, total_size-1 };
	SAInterval next_I = { 0, total_size-1 };

	while (match_length < query_length) {
		next_I = backwardExtend(I, query[match_length]);
		if (next_I.l == next_I.r) break;
		match_length++;
		I = next_I;
	}
	uint_t frequency = I.r - I.l;
	if (frequency > max_anchor_frequency || match_length < min_anchor_length) {
		return 1;
	}

	for (uint_t i = I.l; i < I.r; i++) {
		uint_t ref_pos = getSA(i);
		ChrName name = fasta_manager->getChrName(ref_pos, match_length);
		if (name.size() > 0) {
			region_vec.push_back(Region(name, ref_pos, match_length));
		}

	}
	return match_length;
}

bool FM_Index::saveToFile(const std::string& filename) const
{
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) return false;                       // 打开失败

    cereal::BinaryOutputArchive oar(ofs);         // 也可换 PortableBinary
    oar(*this);                                   // 调用上面 save()
    return static_cast<bool>(ofs);                // 检查流状态
}

bool FM_Index::loadFromFile(const std::string& filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) return false;

    cereal::BinaryInputArchive iar(ifs);
    iar(*this);                                   // 调用上面 load()
    return static_cast<bool>(ifs);
}