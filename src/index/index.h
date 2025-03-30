#ifndef INDEX_H
#define INDEX_H

#include "config.hpp"

using IndexPathMap = std::unordered_map<Species, FilePath>;

enum IndexType {
	R_Index,
	FM_Index
};

class IndexManager {
public:
	FilePath work_dir;
	FilePath index_dir;
	uint_t thread_num;
	IndexManager(const FilePath work_dir, const uint_t thread_num);
	FilePath buildIndex(const std::string prefix, const FilePath reference_path, const IndexType index_type);
};


#endif