#ifndef RARE_ALIGNER_H
#define RARE_ALIGNER_H
#include "config.hpp"
#include "index.h"
#include "threadpool.h"

// 暂时废弃
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

    uint_t thread_num;
    PairRareAligner(const FilePath work_dir, const uint_t thread_num, uint_t chunk_size, uint_t overlap_size, uint_t min_anchor_length, uint_t max_anchor_frequency);
    FilePath buildIndex(const std::string prefix, const FilePath fasta_path, bool fast_build = false);

    MatchVec3DPtr findQueryFileAnchor(const std::string prefix, FastaManager& query_fasta_manager, SearchMode search_mode, bool allow_MEM);

	void clusterPairSpeciesAnchors(MatchVec3DPtr& anchors, FastaManager& query_fasta_manager);

};

// 多基因组比对核心调度类
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

    void starAlignGroup(SpeciesPathMap species_path_map, FilePath hal_path);

};
#endif


