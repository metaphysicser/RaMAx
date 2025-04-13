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

enum IndexType {
	RIndexType,
	FMIndexType,
};

class FM_Index {
public:
	SampledSAType sampled_sa;
	sdsl::wt_huff<sdsl::bit_vector> wt_bwt;

	FM_Index();
	bool buildIndex(FastaManager& fasta_manager, FilePath output_path, bool fast_mode,  uint_t thread);
	bool buildIndexUsingBigBWT(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
	bool read_and_build_sampled_sa(const FilePath& sa_file_path);
	bool read_and_build_bwt(const FilePath& bwt_file_path);
	// bool build_sampled_sa(const FullSAType& full_sa, const std::string& output_path, uint64_t sample_rate);
	bool newScan(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
	bool bwtParse(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
	bool pfBWT(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
	BWTParse::sa_index_t* compute_SA(uint32_t* Text, long n, long k);
};

class IndexManager {
public:
	FilePath work_dir;
	FilePath index_dir;
	uint_t thread_num;
	IndexManager(const FilePath work_dir, const uint_t thread_num);
	FilePath buildIndex(const std::string prefix, FastaManager& ref_fasta_manager, const IndexType index_type);
};


#endif