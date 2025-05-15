#ifndef RARE_ALIGNER_H
#define RARE_ALIGNER_H
#include "config.hpp"
#include "index.h"
#include "threadpool.h"

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
#endif


