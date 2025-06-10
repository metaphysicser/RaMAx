#ifndef RARE_ALIGNER_H
#define RARE_ALIGNER_H
#include "config.hpp"
#include "index.h"
#include "threadpool.h"

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

    void starAlignment(
        uint_t tree_root,
        SearchMode                 search_mode,
        bool                       fast_build,
        bool                       allow_MEM);

    SpeciesMatchVec3DPtrMapPtr alignMultipleGenome(SpeciesName ref_name, SpeciesFastaManagerMap& species_fasta_manager_map, SearchMode search_mode, bool fast_build, bool allow_MEM);

    void filterMultipeSpeciesAnchors(SpeciesName ref_name, SpeciesFastaManagerMap& species_fasta_manager_map, SpeciesMatchVec3DPtrMapPtr species_match_map);
    
   
   
    
};


class PairRareAligner {
public:
    FilePath work_dir;
    FilePath index_dir;

    FastaManager* ref_fasta_manager_ptr;
    FM_Index ref_index;

    uint_t chunk_size;
    uint_t overlap_size;
    uint_t min_anchor_length;
    uint_t max_anchor_frequency;

    uint_t group_id;
    uint_t round_id;

    uint_t thread_num;
    PairRareAligner(const FilePath work_dir, const uint_t thread_num, uint_t chunk_size, uint_t overlap_size, uint_t min_anchor_length, uint_t max_anchor_frequency);

    // 新增构造函数：从 MultipleRareAligner 初始化
    PairRareAligner(const MultipleRareAligner& mra)
        : work_dir(mra.work_dir),
        index_dir(mra.index_dir),
        chunk_size(mra.chunk_size),
        overlap_size(mra.overlap_size),
        min_anchor_length(mra.min_anchor_length),
        max_anchor_frequency(mra.max_anchor_frequency),
        thread_num(mra.thread_num)
    {
        this->group_id = mra.group_id;
        this->round_id = mra.round_id;
    }

    MatchVec3DPtr alignPairGenome(SpeciesName query_name, FastaManager& query_fasta_manager, SearchMode search_mode, bool allow_MEM);
    FilePath buildIndex(const std::string prefix, FastaManager& ref_fasta_manager_, bool fast_build);


    MatchVec3DPtr findQueryFileAnchor(const std::string prefix, FastaManager& query_fasta_manager, SearchMode search_mode, bool allow_MEM, ThreadPool& pool);

    void filterPairSpeciesAnchors(MatchVec3DPtr& anchors, FastaManager& query_fasta_manager);

};


#endif


