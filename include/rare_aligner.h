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

    uint_t thread_num;
    PairRareAligner(const FilePath work_dir, const uint_t thread_num, uint_t chunk_size, uint_t overlap_size, uint_t min_anchor_length, uint_t max_anchor_frequency);
    FilePath buildIndex(const std::string prefix, const FilePath fasta_path, bool fast_build = false);

    MatchVec3DPtr findQueryFileAnchor(const std::string prefix, FastaManager& query_fasta_manager, SearchMode search_mode, bool allow_MEM);

	void filterPairSpeciesAnchors(MatchVec3DPtr& anchors, FastaManager& query_fasta_manager);

};

// 多基因组比对核心调度类
class MultipleRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    SpeciesPathMap species_path_map;
    NewickParser newick_tree;

    uint_t chunk_size;
    uint_t overlap_size;
    uint_t min_anchor_length;
    uint_t max_anchor_frequency;

    uint_t thread_num;
    uint_t group_id;
    uint_t round_id;

    // 构造函数声明：注意名称必须与类名完全一致
    MultipleRareAligner(
        const FilePath& work_dir,
        SpeciesPathMap& species_path_map,
        NewickParser& newick_tree,
        uint_t thread_num,
        uint_t chunk_size,
        uint_t overlap_size,
        uint_t min_anchor_length,
        uint_t max_anchor_frequency
    );

    void starAlignment(uint_t tree_root);

    SpeciesMatchVec3DPtrMapPtr alignMultipleQuerys(SpeciesName ref_name, SpeciesFastaManagerMap& species_fasta_manager_map, SearchMode search_mode, bool fast_build, bool allow_MEM);

    MatchVec3DPtr MultipleRareAligner::alignSingleQuerys(
        const std::string& prefix,
        FM_Index& ref_index,
        FastaManager& query_fm,
        SearchMode         search_mode,
        bool               allow_MEM,
        ThreadPool& pool);
    
   
   
    
};

#endif


