#ifndef INDEX_H
#define INDEX_H

#include "data_process.h"
#include "newscan.hpp"
#include "bwtparse.hpp"
#include "pfbwt.hpp"
extern "C" {
#include "gsa/gsacak.h"
}
#include "CaPS-SA/Suffix_Array.hpp"
#define WINDOW_SIZE 10
#define STOP_MODULUS 100
using IndexPathMap = std::unordered_map<Species, FilePath>;

enum IndexType {
	RIndexType,
	FMIndexType,
};

class FM_Index {
public:
	FM_Index();
	bool buildIndexUsingBigBWT(const FilePath& fasta_path, const FilePath& output_path, uint_t thread);
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