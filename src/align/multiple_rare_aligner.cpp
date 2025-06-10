#include "rare_aligner.h"

MultipleRareAligner::MultipleRareAligner(
    const FilePath& work_dir_,       // 与声明中的类型、顺序一致
    SpeciesPathMap& species_path_map_,
    NewickParser& newick_tree_,
    uint_t thread_num_,              // 同理
    uint_t chunk_size_,
    uint_t overlap_size_,
    uint_t min_anchor_length_,
    uint_t max_anchor_frequency_
)
    : work_dir(work_dir_),                                  // 初始化成员
    index_dir(work_dir_ / INDEX_DIR),
    species_path_map(species_path_map_),
    newick_tree(newick_tree_),
    chunk_size(chunk_size_),
    overlap_size(overlap_size_),
    min_anchor_length(min_anchor_length_),
    max_anchor_frequency(max_anchor_frequency_),
    thread_num(thread_num_)
{
    // 确保工作目录存在
    if (!std::filesystem::exists(work_dir)) {
        std::filesystem::create_directories(work_dir);
        spdlog::info("Created work directory: {}", work_dir.string());
    }

    // 确保索引目录存在（默认放在 work_dir/index）
    if (!std::filesystem::exists(index_dir)) {
        std::filesystem::create_directories(index_dir);
        spdlog::info("Created index directory: {}", index_dir.string());
    }

    this->group_id = 0;
    this->round_id = 0;

}

void MultipleRareAligner::starAlignment(
    uint_t tree_root,
    SearchMode                 search_mode,
    bool                       fast_build,
    bool                       allow_MEM)
{
    std::vector leaf_vec = newick_tree.orderLeavesGreedyMinSum(tree_root);
	uint_t leaf_num = leaf_vec.size();

    // TODO 完成迭代的星比对
    for (uint_t i = 0; i < 1; i++) {
    /*for (uint_t i = 0; i < leaf_num; i++) {*/
		uint_t ref_id = leaf_vec[i];
		SpeciesName ref_name = newick_tree.getNodes()[ref_id].name;
				
        SpeciesFastaManagerMap species_fasta_manager_map;
        for (uint_t j = i; j < leaf_num; j++) {
            uint_t query_id = leaf_vec[j];
			SpeciesName query_name = newick_tree.getNodes()[query_id].name;
			FilePath query_fasta_path = species_path_map[query_name];
			FastaManager query_fasta_manager(query_fasta_path, getFaiIndexPath(query_fasta_path));
            species_fasta_manager_map.emplace(
                query_name,                      // key
                std::move(query_fasta_manager)   // value：右值，触发 move 构造
            );
            
        }

        SpeciesMatchVec3DPtrMapPtr mathc_ptr = alignMultipleGenome(
            ref_name, species_fasta_manager_map,
            ACCURATE_SEARCH, fast_build, allow_MEM
        );

    }
    return;
}

SpeciesMatchVec3DPtrMapPtr MultipleRareAligner::alignMultipleGenome(
    SpeciesName                ref_name,
    SpeciesFastaManagerMap& species_fasta_manager_map,
    SearchMode                 search_mode,
    bool                       fast_build,
    bool                       allow_MEM)
{
    /* ---------- 0. 合法性检查 ---------- */
    if (!species_fasta_manager_map.count(ref_name))
        throw std::runtime_error("[alignMultipleQuerys] reference species not found: " + ref_name);

    if (species_fasta_manager_map.size() <= 1) {
        spdlog::warn("[alignMultipleQuerys] only reference genome present, nothing to align.");
        return std::make_shared<SpeciesMatchVec3DPtrMap>();
    }

    /* ---------- 1. 结果目录与缓存文件 ---------- */
    FilePath result_dir = work_dir / RESULT_DIR
        / ("group_" + std::to_string(group_id))
        / ("round_" + std::to_string(round_id));
    std::filesystem::create_directories(result_dir);

    FilePath anchor_file = result_dir / (ref_name + "_"
        + SearchModeToString(search_mode) + "." + ANCHOR_EXTENSION);

    /* ---------- 2. 如果已存在结果文件直接读取 ---------- */
    if (std::filesystem::exists(anchor_file)) {
        spdlog::info("[alignMultipleQuerys] Load from {}", anchor_file.string());
        auto mp = std::make_shared<SpeciesMatchVec3DPtrMap>();
        if (loadSpeciesMatchMap(anchor_file, mp))
            return mp;
        // 如果读取失败则继续重新计算
    }

    /* ---------- 3. 准备参考基因组索引 ---------- */
    FilePath ref_index_path = index_dir / ref_name;
    std::filesystem::create_directories(ref_index_path);

    PairRareAligner pra(*this);
	pra.buildIndex(ref_name, species_fasta_manager_map[ref_name], fast_build);
	spdlog::info("[alignMultipleQuerys] reference index built for {}.", ref_name);

    /* ---------- 4. 创建共享线程池 ---------- */
    ThreadPool shared_pool(thread_num);

    /* ---------- 5. 为每个 query 物种启动异步任务 ---------- */
    std::unordered_map<SpeciesName, std::future<MatchVec3DPtr>> fut_map;

    for (auto& kv : species_fasta_manager_map) {
        SpeciesName  sp = kv.first;
        if (sp == ref_name) continue;           // 跳过参考

        std::string   prefix = ref_name + "_vs_" + sp;
        FastaManager& fm = kv.second;

        fut_map.emplace(
            sp,
            std::async(std::launch::async,
                [&pra, prefix, &fm, search_mode, allow_MEM, &shared_pool]() -> MatchVec3DPtr {
                    // 统一走单物种比对逻辑，公用 shared_pool
                    return pra.findQueryFileAnchor(prefix, fm, search_mode, allow_MEM, shared_pool);
                })
        );
    }

    shared_pool.waitAllTasksDone();

    /* ---------- 6. 收集所有结果 ---------- */
    auto result_map = std::make_shared<SpeciesMatchVec3DPtrMap>();

    for (auto& kv : fut_map) {
        const SpeciesName& sp = kv.first;
        try {
            MatchVec3DPtr mv3 = kv.second.get();       // 等待并获取
            (*result_map)[sp] = std::move(mv3);
            spdlog::info("[alignMultipleQuerys] {} aligned.", sp);
        }
        catch (const std::exception& e) {
            spdlog::error("[alignMultipleQuerys] {} failed: {}", sp, e.what());
        }
    }

    /* ---------- 7. 保存到文件 ---------- */
    // saveSpeciesMatchMap(anchor_file, result_map);
    spdlog::info("[alignMultipleQuerys] all species done. Saved to {}", anchor_file.string());

    return result_map;
}

void MultipleRareAligner::filterMultipeSpeciesAnchors(
    SpeciesName                ref_name,
    SpeciesFastaManagerMap& species_fasta_manager_map, 
    SpeciesMatchVec3DPtrMapPtr species_match_map
    ) {
    ThreadPool shared_pool(thread_num);
}

