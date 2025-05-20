#ifndef INDEX_H
#define INDEX_H

#include "data_process.h"
#include "anchor.h"
#include "newscan.hpp"
#include "bwtparse.hpp"
#include "pfbwt.hpp"
extern "C" {
#include "gsa/gsacak_wrapper.h"
}
#include "CaPS-SA/Suffix_Array.hpp"
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_huff.hpp>
#include <sdsl/util.hpp>
#include <sdsl/suffix_arrays.hpp>
#include "divsufsort64.h"
#include "divsufsort.h"
#include <cstdlib>

// ----------------------------------------------------------------------
//  Constants
// ----------------------------------------------------------------------
#define WINDOW_SIZE   10u
#define STOP_MODULUS 100u

#define FMINDEX_EXTESION "fmidx"

// ----------------------------------------------------------------------
//  Type aliases
// ----------------------------------------------------------------------
using SampledSA_t = sdsl::int_vector<0>;
// using KmerTable_t = sdsl::int_vector<0>;
using WtHuff_t = sdsl::wt_huff<sdsl::bit_vector>;

// ----------------------------------------------------------------------
//  Suffix‑array 区间
// ----------------------------------------------------------------------
struct SAInterval {
    uint_t l{ 0 };
    uint_t r{ 0 };

    bool empty()    const noexcept { return l == r; }
    bool isUnique() const noexcept { return (r - l) == 1; }
};

enum SearchMode {
    FAST_SEARCH,
    MIDDLE_SEARCH,
    ACCURATE_SEARCH
};

inline const char* SearchModeToString(SearchMode mode) {
    switch (mode) {
    case FAST_SEARCH:     return "fast";
    case MIDDLE_SEARCH:   return "middle";
    case ACCURATE_SEARCH: return "accurate";
    default:              return "unknown";
    }
}

// ----------------------------------------------------------------------
//  FM-Index 类
// ----------------------------------------------------------------------
class FM_Index {
public:
    // -------- 构造 & 析构 --------
    FM_Index() = default;
    FM_Index(SpeciesName species_name, FastaManager* fasta_manager, uint_t sample_rate = 32);

    // -------- 索引构建 --------
    bool buildIndex(FilePath        output_path,
        bool            fast_mode,
        uint_t          thread);

    bool buildIndexUsingBigBWT(
        const FilePath& output_path,
        uint_t          thread);

    bool buildIndexUsingDivsufsort(
        uint_t          thread);

    bool buildIndexUsingCaPS(
        uint_t          thread);

    template<typename index_t>
    bool buildIndexUsingCaPSImpl(const std::string& T, uint_t thread_count);


    bool read_and_build_sampled_sa(const FilePath& sa_file_path);
    bool read_and_build_bwt(const FilePath& bwt_file_path);

    bool newScan(const FilePath& fasta_path,
        const FilePath& output_path,
        uint_t          thread);

    bool bwtParse(const FilePath& fasta_path,
        const FilePath& output_path,
        uint_t          thread);

    bool pfBWT(const FilePath& fasta_path,
        const FilePath& output_path,
        uint_t          thread);

    // -------- SA 操作 --------
    BWTParse::sa_index_t* compute_SA(uint32_t* Text, long n, long k);
    uint_t getSA(uint_t pos)               const;
    uint_t LF(uint_t pos)                  const;
    SAInterval backwardExtend(const SAInterval& I, char c);
    AnchorPtrListVec findAnchors(ChrName query_chr, std::string query, SearchMode search_mode, Strand strand, uint_t query_offset, uint_t min_anchor_length, uint_t max_anchor_frequency);
    AnchorPtrListVec findAnchorsFast(ChrName query_chr, std::string query, Strand strand, uint_t query_offset, uint_t min_anchor_length, uint_t max_anchor_frequency);
    AnchorPtrListVec findAnchorsAccurate(ChrName query_chr, std::string query, Strand strand, uint_t query_offset, uint_t min_anchor_length, uint_t max_anchor_frequency);
    AnchorPtrListVec findAnchorsMiddle(ChrName query_chr, std::string query, Strand strand, uint_t query_offset, uint_t min_anchor_length, uint_t max_anchor_frequency);
    struct MUMInfo {
        uint_t pos;     // query 上的起点
        uint_t len;     // 匹配长度
        bool   is_mum;  // true=MUM, false=MEM
    };

    void bisectAnchors(const std::string& query,
        ChrName   query_chr,
        Strand    strand,
        uint_t    query_offset,
        uint_t    query_length,
        uint_t    min_anchor_length,
        uint_t    max_anchor_frequency,
        const MUMInfo& left,
        const MUMInfo& right,
        AnchorPtrListVec& out);
    uint_t findSubSeqAnchorsFast(const char* query, uint_t query_length, RegionVec& region_vec, uint_t min_anchor_length, uint_t max_anchor_frequency);
    uint_t findSubSeqAnchors(const char* query, uint_t query_length, RegionVec& region_vec, uint_t min_anchor_length, uint_t max_anchor_frequency);
    
    bool saveToFile(const std::string& filename) const;
    bool loadFromFile(const std::string& filename);

   /* bool encode_kmer(const char* kmer, uint_t length, uint64_t code);

    void build_kmer_table(uint_t k, uint_t thread_count);*/




    // -------- 工具函数 --------
    template <size_t N>
    static std::array<char, N> repositionNullAfter(const std::array<char, N>& arr, char c) {
        std::vector<char> temp;
        temp.reserve(N);
        for (char ch : arr) {
            if (ch != '\0') {
                temp.push_back(ch);
            }
        }
        auto it = std::find(temp.begin(), temp.end(), c);
        if (it == temp.end()) {
            throw std::runtime_error("Character not found in array");
        }
        size_t pos = std::distance(temp.begin(), it);
        temp.insert(temp.begin() + pos + 1, '\0');

        if (temp.size() != N) {
            throw std::runtime_error("Unexpected size mismatch after reordering");
        }
        std::array<char, N> result;
        std::copy(temp.begin(), temp.end(), result.begin());
        return result;
    }

    template <class Archive>
    void save(Archive& ar) const
    {
        ar(sample_rate,
            total_size,
            alpha_set,
            alpha_set_without_N,
            count_array,
            char2idx,
            sampled_sa,     // SDSL 的 int_vector 序列化
            wt_bwt         // SDSL 的 wavelet tree 序列化
        );
    }


    template <class Archive>
    void load(Archive& ar)
    {
        ar(sample_rate,
            total_size,
            alpha_set,
            alpha_set_without_N,
            count_array,
            char2idx,
            sampled_sa,
            wt_bwt
        );
    }


    SpeciesName species_name;
    FastaManager* fasta_manager{ nullptr };
    std::array<char, 6> alpha_set{ '\0','A','C','G','N','T' };
    std::array<char, 5> alpha_set_without_N{ '\0','A','C','G','T' };
    SampledSA_t       sampled_sa;
    WtHuff_t          wt_bwt;
    std::array<uint_t, 256>   count_array{};
    uint_t              sample_rate{ 32 };
    uint_t              total_size{ 0 };
    std::array<uint8_t, 256>  char2idx{};
    //uint_t kmer_size;
    //KmerTable_t kmer_table_left;
    //KmerTable_t kmer_table_right;
};

#endif