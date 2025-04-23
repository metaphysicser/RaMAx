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


// ----------------------------------------------------------------------
//  Constants
// ----------------------------------------------------------------------
#define WINDOW_SIZE   10u
#define STOP_MODULUS 100u

// ----------------------------------------------------------------------
//  Type aliases
// ----------------------------------------------------------------------
using IndexPathMap = std::unordered_map<SpeciesName, FilePath>;
using SampledSA_t = sdsl::int_vector<0>;
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

// ----------------------------------------------------------------------
//  FM-Index 类
// ----------------------------------------------------------------------
class FM_Index {
public:
    // -------- 构造 & 析构 --------
    FM_Index() = default;
    FM_Index(SpeciesName species_name, FastaManager* fasta_manager, uint_t sample_rate = 32);

    // -------- 索引构建 --------
    bool buildIndex(FastaManager& fasta_manager,
        FilePath        output_path,
        bool            fast_mode,
        uint_t          thread);

    bool buildIndexUsingBigBWT(const FilePath& fasta_path,
        const FilePath& output_path,
        uint_t          thread);

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
	AnchorPtrVec findAnchors(ChrName query_chr, std::string query, Strand strand, uint_t query_offset, uint_t min_anchor_length, uint_t max_anchor_frequency);
    uint_t findSubSeqAnchors(const char* query, uint_t query_length, RegionVec& region_vec, uint_t min_anchor_length, uint_t max_anchor_frequency);




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
};

// ----------------------------------------------------------------------
//  索引管理类
// ----------------------------------------------------------------------
class IndexManager {
public:
    FilePath work_dir;
    FilePath index_dir;
    uint_t thread_num;
    IndexManager(const FilePath work_dir, const uint_t thread_num);
    FilePath buildIndex(const std::string prefix, FastaManager& ref_fasta_manager);
};
#endif