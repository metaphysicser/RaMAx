#ifndef INDEX_H
#define INDEX_H

#include "data_process.h"
#include "newscan.hpp"
#include "bwtparse.hpp"
#include "pfbwt.hpp"
extern "C" {
#include "gsa/gsacak_wrapper.h"
}
#include "CaPS-SA/Suffix_Array.hpp"
// #include <sdsl/csa_sampling_strategy.hpp>  // sa_order_sa_sampling
//#include <sdsl/bit_vectors.hpp>          // bit_vector
//#include <sdsl/wt_int.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_huff.hpp>
#include <sdsl/util.hpp>
#include <sdsl/suffix_arrays.hpp>
//#include <sdsl/bit_vectors.hpp>


#define WINDOW_SIZE 10
#define STOP_MODULUS 100

using IndexPathMap = std::unordered_map<Species, FilePath>;

using SampledSAType = sdsl::int_vector<0>;
using WtHuffType = sdsl::wt_huff<sdsl::bit_vector>;

struct SAInterval { uint_t l; uint_t r; };
inline bool empty(const SAInterval& I) { return I.l == I.r; }
inline bool isUnique(const SAInterval& I) { return I.r - I.l == 1; }


class FM_Index {
public:
	FastaManager* fasta_manager;
	std::array<char, 6> alpha_set = { '\0', 'A', 'C', 'G', 'N', 'T'};
	std::array<char, 5> alpha_set_without_N = { '\0', 'A', 'C', 'G', 'T' };
	SampledSAType sampled_sa;
	WtHuffType wt_bwt;
	std::array<uint_t, 256> count_array;
	uint_t sample_rate;
	uint_t total_size;
	std::array<uint8_t, 256> char2idx;

	FM_Index();
	FM_Index(FastaManager* fasta_manager, uint_t sample_rate = 32);
	bool buildIndex(FastaManager& fasta_manager, FilePath output_path, bool fast_mode,  uint_t thread);
	bool buildIndexUsingBigBWT(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
	bool read_and_build_sampled_sa(const FilePath& sa_file_path);
	bool read_and_build_bwt(const FilePath& bwt_file_path);
	// bool build_sampled_sa(const FullSAType& full_sa, const std::string& output_path, uint64_t sample_rate);
	bool newScan(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
	bool bwtParse(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
	bool pfBWT(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
	BWTParse::sa_index_t* compute_SA(uint32_t* Text, long n, long k);
	uint_t getSA(uint_t pos) const;
	uint_t LF(uint_t pos) const;

	SAInterval backwardExtend(const SAInterval& I, char c);

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
};


class IndexManager {
public:
	FilePath work_dir;
	FilePath index_dir;
	uint_t thread_num;
	IndexManager(const FilePath work_dir, const uint_t thread_num);
	FilePath buildIndex(const std::string prefix, FastaManager& ref_fasta_manager);
};


#endif