#include "index.h"
#include <sdsl/io.hpp>

// 构造函数：接收物种名、FastaManager 指针、采样率
FM_Index::FM_Index(SpeciesName species_name, SeqPro::ManagerVariant &fasta_manager,
                   uint_t sample_rate)
    : species_name(species_name), fasta_manager(fasta_manager) {
    this->sample_rate = sample_rate;
    
    // 若总序列长度小于采样率，则降采样率为 1，保证后续索引合法
    if (total_size < sample_rate) {
        this->sample_rate = 1;
    }
}

// 使用 CaPS 算法构建索引（适用于较大的序列）
bool FM_Index::buildIndexUsingCaPS(uint_t thread_count) {
    std::string T = std::visit([](auto&& manager_ptr) -> std::string {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if (!manager_ptr) {
            throw std::runtime_error("Manager pointer is null inside variant.");
        }
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
            return manager_ptr->concatAllSequences('\1');
        }
        else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
            return manager_ptr->concatAllSequencesSeparated('\1');
        }
        else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
        }, fasta_manager);
    size_t n = T.size();
    this->total_size = n;

    // 若总序列长度小于采样率，则降采样率为 1，保证后续索引合法
    if (total_size < sample_rate) {
        this->sample_rate = 1;
    }

    if (n == 0)
        return false;

    // 根据序列长度，选择使用 32 位或 64 位的索引类型
    if (n <= std::numeric_limits<uint32_t>::max()) {
        return buildIndexUsingCaPSImpl<uint32_t>(T, thread_count);
    } else {
        return buildIndexUsingCaPSImpl<uint64_t>(T, thread_count);
    }
}

// 模板函数：使用 CaPS 实现构建后缀数组 + BWT + 采样 SA + wavelet tree
template<typename index_t>
bool FM_Index::buildIndexUsingCaPSImpl(const std::string &T,
                                       uint_t thread_count) {
    using SA_t = index_t;
    size_t n = T.size() + 1;

    // 设置 parlay 线程数（用于 CaPS 并行）
    std::string value = std::to_string(thread_count);
    setenv("PARLAY_NUM_THREADS", value.c_str(), 1);

    // 构造后缀数组
    CaPS_SA::Suffix_Array<SA_t> suf_arr(T.c_str(), static_cast<SA_t>(n), 0, 0);
    suf_arr.construct();

    std::vector<SA_t> SA(n); // SA 向量
    std::copy(suf_arr.SA(), suf_arr.SA() + n, SA.begin());

    // 构建 BWT：BWT[i] = T[SA[i] - 1]（循环处理）
    sdsl::int_vector<8> BWT(n);
    ThreadPool pool(thread_count);
    size_t chunk_size = (n + thread_count - 1) / thread_count;

    for (uint_t t = 0; t < thread_count; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, n);
        pool.enqueue([&T, &SA, &BWT, n, start, end]() {
            for (size_t i = start; i < end; ++i) {
                size_t si = static_cast<size_t>(SA[i]);
                BWT[i] = static_cast<uint8_t>(T[(si + n - 1) % n]);
            }
        });
    }
    pool.waitAllTasksDone();

    // 构建采样 SA（每 sample_rate 个采一个）
    size_t sampled_count = (n + sample_rate - 1) / sample_rate;
    sampled_sa.resize(sampled_count);
    for (size_t i = 0, idx = 0; i < n; i += sample_rate, ++idx) {
        sampled_sa[idx] = static_cast<uint64_t>(SA[i]);
    }
    sdsl::util::bit_compress(sampled_sa); // 压缩空间

    // 构建 wavelet tree
    sdsl::construct_im(this->wt_bwt, BWT);

    return true;
}

// 使用 divsufsort 构建索引（适用于中小序列）
bool FM_Index::buildIndexUsingDivsufsort(uint_t thread_count) {
    std::string T = std::visit([](auto &&manager_ptr) -> std::string {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if (!manager_ptr) {
            throw std::runtime_error("Manager pointer is null inside variant.");
        }
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
            return manager_ptr->concatAllSequences('\1');
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
            return manager_ptr->concatAllSequencesSeparated('\1');
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, fasta_manager);
    size_t n = T.size();
    this->total_size = n;

    // 若总序列长度小于采样率，则降采样率为 1，保证后续索引合法
    if (total_size < sample_rate) {
        this->sample_rate = 1;
    }

    if (n == 0)
        return false;

    const auto *Tptr = reinterpret_cast<const sauchar_t *>(T.data());
    std::vector<saidx_t> SA(n); // 后缀数组

    if (divsufsort(Tptr, SA.data(), static_cast<saidx_t>(n)) < 0) {
        spdlog::error("divsufsort() failed");
        return false;
    }

    // 构建 BWT
    sdsl::int_vector<8> BWT(n);
    ThreadPool pool(thread_count);
    size_t chunk_size = (n + thread_count - 1) / thread_count;
    for (uint_t t = 0; t < thread_count; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, n);
        pool.enqueue([&T, &SA, &BWT, n, start, end]() {
            for (size_t i = start; i < end; ++i) {
                saidx_t si = SA[i];
                BWT[i] = static_cast<uint8_t>(T[(static_cast<size_t>(si) + n - 1) % n]);
            }
        });
    }
    pool.waitAllTasksDone();

    // 构建采样 SA
    size_t sampled_count = (n + sample_rate - 1) / sample_rate;
    sampled_sa.resize(sampled_count);
    for (size_t i = 0, idx = 0; i < n; i += sample_rate, ++idx) {
        sampled_sa[idx] = SA[i];
    }
    sdsl::util::bit_compress(sampled_sa);

    // 构建 wavelet tree
    sdsl::construct_im(this->wt_bwt, BWT);

    return true;
}

// 总调度函数：根据模式选择索引构建方式
bool FM_Index::buildIndex(FilePath output_path, bool fast_mode, uint_t thread) {
    std::string fasta_path_str;
    std::visit([&fasta_path_str](auto &&manager_ptr) {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
            fasta_path_str = manager_ptr->getFastaPath();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
            fasta_path_str = manager_ptr->getOriginalManager().getFastaPath();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, fasta_manager);
    if (isFileSmallerThan(fasta_path_str, 1024)) {
        buildIndexUsingDivsufsort(thread);
    } else {
        if (fast_mode) {
            // TODO 对齐divsufsort
            buildIndexUsingCaPS(thread);
        } else {
            buildIndexUsingBigBWT(output_path,
                                  thread); // 大文件处理方式（此处未实现）
        }
    }


    ////// 构建 C 表（前缀计数）和字母转下标映射表
    //size_t cumulative = 0;
    //char2idx.fill(0xFF);
    //count_array.fill(0);

    //bool has_ambiguous_bases = std::visit([](auto &&manager_ptr) -> bool {
    //    using PtrType = std::decay_t<decltype(manager_ptr)>;
    //    if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
    //        return manager_ptr->hasAmbiguousBasesAll();
    //    } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager> >) {
    //        return manager_ptr->getOriginalManager().hasAmbiguousBasesAll();
    //    } else {
    //        throw std::runtime_error("Unhandled manager type in variant.");
    //    }
    //}, fasta_manager);
    //if (has_ambiguous_bases) {
    //    uint_t count = 0;
    //    for (const auto &c: alpha_set) {
    //        size_t occ = wt_bwt.rank(wt_bwt.size(), c);
    //        count_array[static_cast<uint8_t>(c)] = cumulative;
    //        cumulative += occ;
    //        char2idx[static_cast<uint8_t>(c)] = count;
    //        count++;
    //    }
    //} else {
    //    uint_t count = 0;
    //    for (const auto &c: alpha_set_without_N) {
    //        size_t occ = wt_bwt.rank(wt_bwt.size(), c);
    //        count_array[static_cast<uint8_t>(c)] = cumulative;
    //        cumulative += occ;
    //        char2idx[static_cast<uint8_t>(c)] = count;
    //        count++;
    //    }
    //}
    ///---------- 1. 判断是否需要把 'N' 纳入字母表 ----------
    bool has_ambiguous_bases = std::visit([](auto&& manager_ptr) -> bool {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>)
            return manager_ptr->hasAmbiguousBasesAll();          // 原始 FASTA
        else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>)
            return manager_ptr->getOriginalManager().hasAmbiguousBasesAll(); // Masked 版本
        else
            throw std::runtime_error("Unhandled manager type in variant.");
        }, fasta_manager);

    // ---------- 2. 选定实际使用的字母集合 ----------
    const std::span<const char> alph =
        has_ambiguous_bases
        ? std::span<const char>(alpha_set.data(), alpha_set.size())
        : std::span<const char>(alpha_set_without_N.data(), alpha_set_without_N.size());
    // alpha_set           = {'\0','\1','A','C','G','T','N'};
    // alpha_set_without_N = {'\0','\1','A','C','G','T'};

    // ---------- 3. 初始化数据结构 ----------
    size_t cumulative = 0;          // C[c] 的前缀计数游标
    char2idx.fill(0xFF);            // 0xFF 表示“未映射”
    count_array.fill(0);            // C 表全部先清零

    // ---------- 4. 顺序遍历 alph，填充 C 表和 char→rank 映射 ----------
    uint_t rank = 0;                // 连续的字符 rank：0,1,2,...
    for (char c : alph) {
        // BWT 中字符 c 出现的总次数
        size_t occ = wt_bwt.rank(wt_bwt.size(), c);

        // C[c]：小于 c 的所有字符出现总数
        count_array[static_cast<uint8_t>(c)] = cumulative;
        cumulative += occ;          // 更新前缀和

        // char2idx：字符 → rank，下游 wavelet‑tree 查询要用
        char2idx[static_cast<uint8_t>(c)] = rank++;
    }

    return true;
}

// FM 索引的 LF 映射操作：获取前一位置
uint_t FM_Index::LF(uint_t pos) const {
    uint8_t c = wt_bwt[pos];
    return count_array[c] + wt_bwt.rank(pos, c);
}

// 根据位置 pos 向前跳转，恢复 SA 值
uint_t FM_Index::getSA(size_t pos) const {
    uint_t steps = 0;
    size_t cur = pos;
    while ((cur) % sample_rate != 0) {
        cur = LF(cur);
        ++steps;
    }
    size_t sa_sample = sampled_sa[cur / sample_rate];
    size_t sa_val = sa_sample + steps;
    if (sa_val >= total_size)
        sa_val -= total_size;
    return sa_val;
}

// 反向扩展（Backward Search 的核心）：区间 I 按字符 c 向前扩展
SAInterval FM_Index::backwardExtend(const SAInterval &I, char c) {
    uint8_t ch = static_cast<uint8_t>(c);
    uint_t new_l = count_array[ch] + wt_bwt.rank(I.l, ch);
    uint_t new_r = count_array[ch] + wt_bwt.rank(I.r, ch);
    return {new_l, new_r};
}

// 根据查询模式调度查找函数
MatchVec2DPtr FM_Index::findAnchors(ChrName query_chr, std::string query,
                                    SearchMode search_mode, Strand strand,
                                    bool allow_MEM, uint_t query_offset,
                                    uint_t min_anchor_length,
                                    bool allow_short_mum,
                                    uint_t max_anchor_frequency,
                                    sdsl::int_vector<0> &ref_global_cache,
                                    SeqPro::Length sampling_interval) {
    if (search_mode == FAST_SEARCH) {
        return findAnchorsFast(query_chr, query, strand, allow_MEM, query_offset,
                               min_anchor_length, allow_short_mum, max_anchor_frequency, ref_global_cache, sampling_interval);
    } else if (search_mode == ACCURATE_SEARCH) {
        return findAnchorsAccurate(query_chr, query, strand, allow_MEM,
                                   query_offset, min_anchor_length, allow_short_mum,
                                   max_anchor_frequency, ref_global_cache, sampling_interval);
    } else if (search_mode == MIDDLE_SEARCH) {
        return findAnchorsMiddle(query_chr, query, strand, allow_MEM, query_offset,
                                 min_anchor_length, allow_short_mum, max_anchor_frequency, ref_global_cache, sampling_interval);
    } else {
        throw std::invalid_argument("Invalid search mode");
    }
}

// 快速查找模式：逐段匹配 query 子串
MatchVec2DPtr FM_Index::findAnchorsMiddle(ChrName query_chr, std::string query, Strand strand, bool allow_MEM,
                                        uint_t query_offset, uint_t min_anchor_length, bool allow_short_mum, uint_t max_anchor_frequency,
                                        sdsl::int_vector<0> &ref_global_cache,
                                        SeqPro::Length sampling_interval) {
    // query 字符串反向或互补（FM 索引默认支持反向搜索）
    if (strand == FORWARD) {
        std::reverse(query.begin(), query.end());
    } else {
        //std::reverse(query.begin(), query.end());
        for (char &ch: query) {
            ch = BASE_COMPLEMENT[static_cast<unsigned char>(ch)];
        }
    }

    MatchVec2DPtr anchor_ptr_list_vec = std::make_shared<MatchVec2D>();

    uint_t total_length = 0;
    uint_t query_length = query.length();

    // 主循环：从前往后逐段查找 anchor
    while (total_length < query_length) {
        RegionVec region_vec;
        uint_t match_length = findSubSeqAnchors(
            query.c_str() + total_length, query_length - total_length, allow_MEM,
            region_vec, min_anchor_length, allow_short_mum, max_anchor_frequency, ref_global_cache, sampling_interval);


        if (!region_vec.empty()) {
            Region query_region;
            if (strand == FORWARD) {
                query_region = Region(query_chr, query_length - match_length - total_length + query_offset,
                                      match_length);
            } else {
                query_region = Region(query_chr, total_length + query_offset, match_length);
            }

            MatchVec anchor_ptr_list;
            anchor_ptr_list.reserve(region_vec.size());
            for (uint_t i = 0; i < region_vec.size(); i++) {
                // Match match(region_vec[i], query_region, strand);
                Match match(region_vec[i].chr_name, region_vec[i].start, query_chr, query_region.start, match_length, strand);
                /*Score_t score = caculateMatchScore(query.c_str() + total_length,
                match_length); Cigar_t cigar; cigar.push_back(cigarToInt('=',
                match_length));*/
                // Anchor p = Anchor(match, match_length, cigar, score);
                anchor_ptr_list.push_back(match);
            }

            anchor_ptr_list_vec->push_back(anchor_ptr_list);
        }

        total_length += match_length;
    }
    anchor_ptr_list_vec->shrink_to_fit();
    return anchor_ptr_list_vec;
}

// 中速模式：跳跃推进，减少重复匹配
MatchVec2DPtr FM_Index::findAnchorsFast(ChrName query_chr, std::string query,
    Strand strand, bool allow_MEM,
    uint_t query_offset,
    uint_t min_anchor_length,
    bool allow_short_mum,
    uint_t max_anchor_frequency,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval) {
    if (strand == FORWARD) {
        std::reverse(query.begin(), query.end());
    }
    else {
        for (char& ch : query) {
            ch = BASE_COMPLEMENT[static_cast<unsigned char>(ch)];
        }
    }

    MatchVec2DPtr anchor_ptr_list_vec = std::make_shared<MatchVec2D>();

    uint_t total_length = 0;
    uint_t query_length = query.length();
    uint_t last_pos = 0;

    while (total_length < query_length) {
        RegionVec region_vec;
        uint_t match_length = findSubSeqAnchors(
            query.c_str() + total_length, query_length - total_length, allow_MEM,
            region_vec, min_anchor_length, allow_short_mum, max_anchor_frequency, ref_global_cache, sampling_interval);

        if (!region_vec.empty()) {
            Region query_region;
            if (strand == FORWARD) {
                query_region = Region(query_chr, query_length - match_length - total_length + query_offset,
                    match_length);
            }
            else {
                query_region = Region(query_chr, total_length + query_offset, match_length);
            }
            MatchVec anchor_ptr_list;
			anchor_ptr_list.reserve(region_vec.size());
            uint_t ref_end_pos = region_vec[0].start + match_length;

            // 仅当与上一段区域不重复时才添加
            if (ref_end_pos != last_pos) {
                for (uint_t i = 0; i < region_vec.size(); i++) {
                    Match match(region_vec[i].chr_name, region_vec[i].start, query_chr, query_region.start, match_length, strand);
                    /* Score_t score = caculateMatchScore(query.c_str() + total_length,
                     match_length); Cigar_t cigar; cigar.push_back(cigarToInt('=',
                     match_length)); Anchor p = Anchor(match, match_length, cigar,
                     score);*/
                    anchor_ptr_list.push_back(match);
                }
                anchor_ptr_list_vec->push_back(anchor_ptr_list);
            }
            last_pos = ref_end_pos;
        }
        total_length += std::min(min_anchor_length, match_length == 1 ? min_anchor_length : match_length);
    }
    anchor_ptr_list_vec->shrink_to_fit();

    return anchor_ptr_list_vec;
}

// -------------------------------------------------------------
// Middle-Search：先 fast，后区间二分补搜
// -------------------------------------------------------------
MatchVec2DPtr FM_Index::findAnchorsAccurate(ChrName query_chr,
                                            std::string query, Strand strand,
                                            bool allow_MEM, uint_t query_offset,
                                            uint_t min_anchor_length,
                                            bool allow_short_mum,
                                            uint_t max_anchor_frequency,
                                            sdsl::int_vector<0> &ref_global_cache,
                                            SeqPro::Length sampling_interval) {
    MatchVec2DPtr out = std::make_shared<MatchVec2D>();

    uint_t query_length = query.length();

    /* ---------- STEP-0：先跑 Fast 拿到顶层 MUM ---------- */
    MatchVec2DPtr fast_mums =
            findAnchorsFast(query_chr, query, strand, allow_MEM, query_offset,
                            min_anchor_length, max_anchor_frequency, allow_short_mum, ref_global_cache, sampling_interval);

    if (fast_mums->empty())
        return out;

    /* ---------- 方向处理，与 fast 保持一致 ---------- */
    if (strand == FORWARD)
        std::reverse(query.begin(), query.end());
    else
        for (char &ch: query)
            ch = BASE_COMPLEMENT[static_cast<unsigned char>(ch)];

    out->reserve(fast_mums->size());

    uint_t count;
    for (const auto &lst: *fast_mums) {
        bool left_is_mum = (lst.size() == 1);

        const Match &a = lst.front();
        // start = query_length - match_length - total_length + query_offset
        // total_length = query_length - match_length - start + query_offset
        uint_t L = query_length - a.match_len() - a.qry_start +
                   query_offset;
        uint_t n = a.match_len();
        uint_t R = L + n;

        out->push_back(lst); // fast 结果入库

        // 右端位置（R-1）再做一次搜索，作为 bisect 的 right 端点
        uint_t right_pos = R - 1;
        RegionVec regs;
        uint_t right_len = findSubSeqAnchors(
            query.c_str() + right_pos, query_length - right_pos, allow_MEM, regs,
            min_anchor_length, allow_short_mum, max_anchor_frequency, ref_global_cache, sampling_interval);

        if (!regs.empty()) {
            Region query_region;
            if (strand == FORWARD) {
                query_region = Region(query_chr, query_length - right_len - right_pos + query_offset, right_len);
            } else {
                query_region = Region(query_chr, query_offset + right_pos, right_len);
            }

            MatchVec anchor_ptr_list;
            anchor_ptr_list.reserve(regs.size());
            for (uint_t i = 0; i < regs.size(); i++) {
                Match match(regs[i].chr_name, regs[i].start, query_chr, query_region.start, right_len, strand);
                /*Score_t score = caculateMatchScore(query.c_str() + right_pos,
                right_len); Cigar_t cigar; cigar.push_back(cigarToInt('=', right_len));
                Anchor p = Anchor(match, right_len, cigar, score);*/
                anchor_ptr_list.push_back(match);
            }

            out->push_back(anchor_ptr_list);
        }

        bool right_is_mum = (regs.size() == 1);

        // 调用成员函数递归补搜
        bisectAnchors(query, query_chr, strand, allow_MEM, query_offset,
                      query_length, min_anchor_length, allow_short_mum, max_anchor_frequency,
                      {L, n, left_is_mum}, {right_pos, right_len, right_is_mum},
                      *out, ref_global_cache, sampling_interval);
        count++;
    }
    out->shrink_to_fit();
    return out;
}

// uint_t FM_Index::findSubSeqAnchorsFast(const char* query, uint_t
// query_length, RegionVec& region_vec, uint_t min_anchor_length, uint_t
// max_anchor_frequency)
//{
//	uint_t left = 0;
//	uint_t right = 0;
//	uint64_t kmer_index = 0;
//
//	if (kmer_size <= query_length) {
//		bool encode_success = encode_kmer(query, kmer_size, kmer_index);
//		if (encode_success) {
//			left = kmer_table_left[kmer_index];
//			right = kmer_table_right[kmer_index];
//		}
//	}
//
//	uint_t match_length = 0;
//	SAInterval I = { 0, total_size - 1 };
//	SAInterval next_I = { 0, total_size - 1 };
//
//	if (left != 0 || right != 0) {
//		match_length = kmer_size;
//		I = { left, right };
//		next_I = { left, right };
//	}
//
//	while (match_length < query_length) {
//		next_I = backwardExtend(I, query[match_length]);
//		if (next_I.l == next_I.r) break;
//		match_length++;
//		I = next_I;
//	}
//	uint_t frequency = I.r - I.l;
//	if (frequency > max_anchor_frequency || match_length <
// min_anchor_length) { 		return 1;
//	}
//
//	for (uint_t i = I.l; i < I.r; i++) {
//		uint_t ref_pos = getSA(i);
//		ChrName name = fasta_manager->getChrName(ref_pos, match_length);
//		if (name.size() > 0) {
//			region_vec.push_back(Region(name, ref_pos,
// match_length));
//		}
//
//	}
//	return match_length;
// }

// -------------------------------------------------------------
// 递归二分：把 (left.pos, right.pos) 区间彻底搜索干净
// -------------------------------------------------------------
void FM_Index::bisectAnchors(const std::string &query, ChrName query_chr,
                             Strand strand, bool allow_MEM, uint_t query_offset,
                             uint_t query_length, uint_t min_len, bool allow_short_mum,
                             uint_t max_freq, const MUMInfo &left,
                             const MUMInfo &right, MatchVec2D &out,
                             sdsl::int_vector<0> &ref_global_cache,
                             uint_t sampling_interval) {
    if (right.pos <= left.pos + 1)
        return; // 区间不足 1bp

    uint_t mid = left.pos + (right.pos - left.pos) / 2;

    RegionVec regs;
    uint_t mid_len = findSubSeqAnchors(query.c_str() + mid, query_length - mid,
                                       allow_MEM, regs, min_len, allow_short_mum, max_freq, ref_global_cache, sampling_interval);

    bool mid_is_mum = (regs.size() == 1);
    bool same_as_left = (left.pos + left.len == mid + mid_len);
    bool same_as_right = (right.pos + right.len == mid + mid_len);
    // 如果mid_is_mum是false，写入，如果是true，判断same_as_left和same_as_right都为false则写入，否则不写入
    /* ---- 若找到匹配就写结果 ---- */
    if ((regs.size() == 1 || allow_MEM) && !mid_is_mum || (!same_as_left && !same_as_right)) {
        Region qreg;
        if (strand == FORWARD) {
            qreg = Region(query_chr, query_length - mid_len - mid + query_offset, mid_len);
        } else {
            qreg = Region(query_chr, query_offset + mid, mid_len);
        }

        MatchVec lst;
        lst.reserve(regs.size());
        for (auto const &rg: regs) {
            Match match(rg.chr_name, rg.start, query_chr, qreg.start, mid_len, strand);
            /*Score_t sc = caculateMatchScore(query.c_str() + mid, mid_len);
            Cigar_t cg = { cigarToInt('=', mid_len) };*/
            lst.emplace_back(match);
        }
        out.emplace_back(std::move(lst));
    }

    /* ---- 递归左侧 ---- */
    if (!(left.is_mum && same_as_left)) {
        bisectAnchors(query, query_chr, strand, allow_MEM, query_offset,
                      query_length, min_len, allow_short_mum, max_freq, left,
                      {mid, mid_len, mid_is_mum}, out, ref_global_cache, sampling_interval);
    }

    /* ---- 递归右侧 ---- */
    if (!(right.is_mum && same_as_right)) {
        bisectAnchors(query, query_chr, strand, allow_MEM, query_offset,
                      query_length, min_len, allow_short_mum, max_freq, {mid, mid_len, mid_is_mum},
                      right, out, ref_global_cache, sampling_interval);
    }
}

uint_t FM_Index::findSubSeqAnchors(const char *query, uint_t query_length,
                                   bool allow_MEM, RegionVec &region_vec,
                                   uint_t min_anchor_length,
                                   bool allow_short_mum,
                                   uint_t max_anchor_frequency,
                                   sdsl::int_vector<0> &ref_global_cache,
                                   SeqPro::Length sampling_interval) {
    uint_t match_length = 0;
    SAInterval I = {0, total_size - 1};
    SAInterval next_I = {0, total_size - 1};

    while (match_length < query_length) {
        next_I = backwardExtend(I, query[match_length]);
        if (next_I.l == next_I.r) break;
        match_length++;
        I = next_I;
    }
    uint_t frequency = I.r - I.l;
    if (frequency > max_anchor_frequency) {
        return 1;
    }
    if (allow_MEM == false && frequency > 1) {
        return 1;
    }

	if (match_length < min_anchor_length && allow_short_mum == false) {
		return 1;
	}
    //if (frequency > max_anchor_frequency) {
    //	return 1;
    //}
    if (true) {
        region_vec.reserve(region_vec.size() + frequency);
        // 使用传入的采样间隔参数
        auto cache_size = ref_global_cache.size();

        for (uint_t i = I.l; i < I.r; i++) {
            uint_t ref_global_pos = getSA(i);

            std::visit([&](auto &&manager_ptr) {
                SeqPro::SequenceId seq_id = SeqPro::SequenceIndex::INVALID_ID;
                SeqPro::Position local_pos = 0;
                // TODO 用的时候得加个验证看你要转的global坐标和这条染色体的起始全局坐标的差是否小于这条染色体的长度，如果大于的话得切到下一条染色体
                // TODO 尝试使用缓存快速获取序列ID，避免二分搜索
                // if (cache_size > 0) {
                //     auto cache_index = ref_global_pos / sampling_interval;

                //     if (cache_index < cache_size) {
                //         auto candidate_seq_id = ref_global_cache[cache_index];

                //         if (candidate_seq_id != SeqPro::SequenceIndex::INVALID_ID) {
                //             // 获取候选序列的信息进行快速验证
                //             using PtrType = std::decay_t<decltype(manager_ptr)>;
                //             const SeqPro::SequenceInfo *candidate_info = nullptr;

                //             if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                //                 candidate_info = manager_ptr->getIndex().getSequenceInfo(candidate_seq_id);
                //             } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                //                 SeqPro::MaskedSequenceManager> >) {
                //                 candidate_info = manager_ptr->getOriginalManager().getIndex().getSequenceInfo(
                //                     candidate_seq_id);
                //             }

                //             // 快速验证：检查全局坐标是否在该序列范围内
                //             if (candidate_info &&
                //                 ref_global_pos >= candidate_info->global_start_pos &&
                //                 ref_global_pos < candidate_info->global_start_pos + candidate_info->length) {
                //                 // 缓存命中！直接计算局部坐标
                //                 seq_id = candidate_seq_id;
                //                 //local_pos = ref_global_pos - candidate_info->global_start_pos;
                                
                //                 if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                //                     auto [id, pos] = manager_ptr->globalToLocal(ref_global_pos);
                //                     local_pos = pos;
                //                     region_vec.emplace_back(candidate_info->name, local_pos, match_length);
                //                 }
                //                 else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                //                     SeqPro::MaskedSequenceManager> >) {
                //                     auto [fallback_seq_name, fallback_local_pos] = manager_ptr->globalToLocalSeparated(ref_global_pos);
                //                     region_vec.emplace_back(fallback_seq_name, fallback_local_pos, match_length);
                //                 }
                                
                //             }
                //         }
                //     }
                // }

                // 如果缓存未命中，回退到原始的二分搜索
                if (seq_id == SeqPro::SequenceIndex::INVALID_ID) {
                    using PtrType = std::decay_t<decltype(manager_ptr)>;
                    if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                        const SeqPro::SequenceInfo* candidate_info = nullptr;
                        auto [fallback_seq_id, fallback_local_pos] = manager_ptr->globalToLocal(ref_global_pos);
                        candidate_info = manager_ptr->getIndex().getSequenceInfo(fallback_seq_id);
                        region_vec.emplace_back(candidate_info->name, fallback_local_pos, match_length);
                        //region_vec.emplace_back(candidate_info->name, ref_global_pos, match_length);

                    }
                    else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                        SeqPro::MaskedSequenceManager> >) {
                        auto [fallback_seq_name, fallback_local_pos] = manager_ptr->globalToLocalSeparated(ref_global_pos);
                        region_vec.emplace_back(fallback_seq_name, fallback_local_pos, match_length);

                        //region_vec.emplace_back(fallback_seq_name, ref_global_pos, match_length);
                    }
                }

                if (seq_id == SeqPro::SequenceIndex::INVALID_ID) {
                    return;
                }
            }, fasta_manager);
        }
    }

    return match_length;
}

bool FM_Index::saveToFile(const std::string &filename) const {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) return false;

    // 保存基本信息
    ofs.write(reinterpret_cast<const char *>(&sample_rate), sizeof(sample_rate));
    ofs.write(reinterpret_cast<const char *>(&total_size), sizeof(total_size));
    ofs.write(reinterpret_cast<const char *>(&alpha_set), sizeof(alpha_set));
    ofs.write(reinterpret_cast<const char *>(&alpha_set_without_N), sizeof(alpha_set_without_N));
    ofs.write(reinterpret_cast<const char *>(&count_array), sizeof(count_array));
    ofs.write(reinterpret_cast<const char *>(&char2idx), sizeof(char2idx));

    // 使用SDSL的store_to_file函数
    std::string sa_file = filename + ".sa";
    std::string wt_file = filename + ".wt";

    if (!sdsl::store_to_file(sampled_sa, sa_file)) {
        return false;
    }
    if (!sdsl::store_to_file(wt_bwt, wt_file)) {
        return false;
    }

    return static_cast<bool>(ofs);
}

bool FM_Index::loadFromFile(const std::string &filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) return false;

    // 加载基本信息
    ifs.read(reinterpret_cast<char *>(&sample_rate), sizeof(sample_rate));
    ifs.read(reinterpret_cast<char *>(&total_size), sizeof(total_size));
    ifs.read(reinterpret_cast<char *>(&alpha_set), sizeof(alpha_set));
    ifs.read(reinterpret_cast<char *>(&alpha_set_without_N), sizeof(alpha_set_without_N));
    ifs.read(reinterpret_cast<char *>(&count_array), sizeof(count_array));
    ifs.read(reinterpret_cast<char *>(&char2idx), sizeof(char2idx));

    // 使用SDSL的load_from_file函数
    std::string sa_file = filename + ".sa";
    std::string wt_file = filename + ".wt";

    if (!sdsl::load_from_file(sampled_sa, sa_file)) {
        return false;
    }
    if (!sdsl::load_from_file(wt_bwt, wt_file)) {
        return false;
    }

    return static_cast<bool>(ifs);
}

// bool FM_Index::encode_kmer(const char* kmer, uint_t length, uint64_t code) {
//	code = 0;
//	for (uint_t i = 0; i < length; ++i) {
//		code <<= 2;
//		switch (kmer[i]) {
//		case 'A': code |= 0; break;
//		case 'C': code |= 1; break;
//		case 'G': code |= 2; break;
//		case 'T': code |= 3; break;
//		default:
//			return false;
//		}
//	}
//	return true;
// }
//
//
// void FM_Index::build_kmer_table(uint_t k, uint_t thread_count) {
//	const size_t table_size = 1ULL << (2 * k);  // 4^k
//	// std::vector<SAInterval> kmer_table;
//	kmer_table_left.resize(table_size);
//	kmer_table_right.resize(table_size);
//	const char base_table[4] = { 'A', 'C', 'G', 'T' };
//
//	ThreadPool pool(thread_count);
//
//	// 分段处理，每个任务处理一段 k-mer code 区间
//	const size_t block_size = (table_size + thread_count - 1) /
// thread_count;  // 可调参数 	for (size_t block_start = 0; block_start <
// table_size; block_start += block_size) { 		size_t block_end =
// std::min(block_start + block_size, table_size);
//
//		pool.enqueue([=, &base_table]() {
//			std::string kmer(k, 'A');  // 每线程本地缓冲，避免共享
//
//			for (uint64_t code = block_start; code < block_end;
//++code) { 				uint64_t temp = code;
//
//				// 解码整数 → k-mer 字符串
//				for (int i = k - 1; i >= 0; --i) {
//					kmer[i] = base_table[temp & 0b11];
//					temp >>= 2;
//				}
//
//				SAInterval interval = { 0, total_size - 1 };
//				SAInterval notfind_interval = { 0, 0 };
//				bool find = true;
//				for (int i = k - 1; i >= 0 && !interval.empty();
//--i) { 					interval =
// backwardExtend(interval, kmer[i]); 					if
//(interval.r == interval.l) { 						find =
// false;
//						// kmer_table[code] =
// notfind_interval;
// kmer_table_left[code] = 0;
// kmer_table_right[code] = 0; 						break;
//					}
//				}
//
//				// 空区间表示未命中
//				if (find) {
//					kmer_table_left[code] = interval.l;
//					kmer_table_right[code] = interval.r;
//				}
//
//			}
//			});
//	}
//
//	pool.waitAllTasksDone();
//
//	sdsl::util::bit_compress(kmer_table_left);
//	sdsl::util::bit_compress(kmer_table_right);
//
//	return;
// }
