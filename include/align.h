#ifndef ALIGN_H
#define ALIGN_H
#include "config.hpp"
#include "index.h"

class PairRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    FastaManager ref_fasta_manager;
    uint_t thread_num;
    PairRareAligner(const FilePath work_dir, const uint_t thread_num);
    FilePath buildIndex(const std::string prefix, const FilePath fasta_path);
    void alignQueryFile(const std::string prefix, const FilePath fasta_path);
};
#endif


