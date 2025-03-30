
#include "index.h"

IndexManager::IndexManager(const FilePath work_dir, const uint_t thread_num) {
	this->work_dir = work_dir;
	this->index_dir = work_dir / INDEX_DIR;
	this->thread_num = thread_num;
}

FilePath IndexManager::buildIndex(const std::string prefix, const FilePath reference_path, const IndexType index_type) {
	FilePath index_path = index_dir / prefix;
	return index_path;

}