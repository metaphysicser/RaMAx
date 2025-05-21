#ifndef RARE_ALIGNER_H
#define RARE_ALIGNER_H
#include "config.hpp"
#include "index.h"
#include "threadpool.h"

using SpeciesPathMap = std::unordered_map<SpeciesName, FilePath>;

// 此类暂停开发
class PairRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    FastaManager ref_fasta_manager;
    FM_Index ref_index;

    uint_t chunk_size;
    uint_t overlap_size;
    uint_t min_anchor_length;
    uint_t max_anchor_frequency;

	AnchorPtrListVec unique_anchors;
    AnchorPtrListVec repeat_anchors;

    uint_t thread_num;
    PairRareAligner(const FilePath work_dir, const uint_t thread_num, uint_t chunk_size, uint_t overlap_size, uint_t min_anchor_length, uint_t max_anchor_frequency);
    FilePath buildIndex(const std::string prefix, const FilePath fasta_path, bool fast_build = false);

    AnchorPtrListVec findQueryFileAnchor(const std::string prefix, const FilePath fasta_path, SearchMode search_mode);
	void filterAnchors(AnchorPtrListVec& anchors);

};

// 这个类用来完成多基因组比对，tqz负责完成整体框架，zpl负责starAlignGroup
class MultipeRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    FastaManager ref_fasta_manager;
    FM_Index ref_index;

    uint_t chunk_size;
    uint_t overlap_size;
    uint_t min_anchor_length;
    uint_t max_anchor_frequency;

    AnchorPtrListVec unique_anchors;
    AnchorPtrListVec repeat_anchors;

    uint_t thread_num;
    MultipeRareAligner(const FilePath work_dir, const uint_t thread_num, uint_t chunk_size, uint_t overlap_size, uint_t min_anchor_length, uint_t max_anchor_frequency);
    FilePath buildIndex(SpeciesPathMap species_path_map, bool fast_build = false);

    // 这个函数就是星比对的入口
    void starAlignGroup(SpeciesPathMap species_path_map, FilePath hal_path);

};
#endif


