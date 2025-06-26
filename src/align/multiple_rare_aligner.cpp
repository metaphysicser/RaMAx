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
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    uint_t tree_root,
    SearchMode                 search_mode,
    bool                       fast_build,
    bool                       allow_MEM,
    bool                       mask_mode,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval)
{
    std::vector leaf_vec = newick_tree.orderLeavesGreedyMinSum(tree_root);
	uint_t leaf_num = leaf_vec.size();

    // TODO 完成迭代的星比对
    for (uint_t i = 0; i < 1; i++) {
    /*for (uint_t i = 0; i < leaf_num; i++) {*/
		uint_t ref_id = leaf_vec[i];
		SpeciesName ref_name = newick_tree.getNodes()[ref_id].name;
        std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant> species_fasta_manager_map;
        for (uint_t j = i; j < leaf_num; j++) {
            uint_t query_id = leaf_vec[j];
			SpeciesName query_name = newick_tree.getNodes()[query_id].name;
            auto query_fasta_manager = seqpro_managers.at(query_name);
            species_fasta_manager_map.emplace(query_name, query_fasta_manager);
        }

        SpeciesMatchVec3DPtrMapPtr match_ptr = alignMultipleGenome(
            ref_name, species_fasta_manager_map,
            ACCURATE_SEARCH, fast_build, allow_MEM, ref_global_cache, sampling_interval
        );

		filterMultipeSpeciesAnchors(
			ref_name, species_fasta_manager_map, match_ptr
		);

    }
    return;
}

SpeciesMatchVec3DPtrMapPtr MultipleRareAligner::alignMultipleGenome(
    SpeciesName                ref_name,
    std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fasta_manager_map,
    SearchMode                 search_mode,
    bool                       fast_build,
    bool                       allow_MEM,
    sdsl::int_vector<0>& ref_global_cache,
    SeqPro::Length sampling_interval)
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
	pra.buildIndex(ref_name, *species_fasta_manager_map[ref_name], fast_build);
	spdlog::info("[alignMultipleQuerys] reference index built for {}.", ref_name);

    /* ---------- 4. 创建共享线程池 ---------- */
    ThreadPool shared_pool(thread_num);

    /* ---------- 5. 为每个 query 物种启动异步任务 ---------- */
    std::unordered_map<SpeciesName, std::future<MatchVec3DPtr>> fut_map;

    for (auto& kv : species_fasta_manager_map) {
        SpeciesName  sp = kv.first;
        if (sp == ref_name) continue;           // 跳过参考

        std::string   prefix = ref_name + "_vs_" + sp;
        auto& fm = kv.second;

        fut_map.emplace(
            sp,
            std::async(std::launch::async,
                [&pra, prefix, &fm, search_mode, allow_MEM, &shared_pool, &ref_global_cache, sampling_interval]() -> MatchVec3DPtr {
                    // 在多基因组对齐中，每个查询物种使用独立的空缓存（暂时禁用缓存功能）
                    return pra.findQueryFileAnchor(prefix, *fm, search_mode, allow_MEM, shared_pool, ref_global_cache, sampling_interval);
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
    SpeciesName                       ref_name,
    std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fm_map,
    SpeciesMatchVec3DPtrMapPtr        species_match_map)
{
    if (!species_match_map || species_match_map->empty()) return;


    ThreadPool pool(thread_num);                      // 全局线程池

    /*-------------- 预分配表 ----------------------------------------*/
    std::unordered_map<SpeciesName, MatchByStrandByQueryRefPtr> unique_map;
    std::unordered_map<SpeciesName, MatchByStrandByQueryRefPtr> repeat_map;
    SpeciesClusterMap cluster_map;                    // 最终结果

    /*========================= Phase-1  : group =====================*/
    for (auto it = species_match_map->begin();
        it != species_match_map->end(); ++it)
    {
        const SpeciesName  species = it->first;
        if (species == ref_name) continue;

        MatchVec3DPtr mv3_ptr = it->second;
        auto& qfm = species_fm_map.at(species);
        auto& rfm = species_fm_map.at(ref_name);

        auto u_ptr = std::make_shared<MatchByStrandByQueryRef>();
        auto r_ptr = std::make_shared<MatchByStrandByQueryRef>();
        unique_map[species] = u_ptr;
        repeat_map[species] = r_ptr;

        // 把所有需要的变量显式捕获
        pool.enqueue(
            [&mv3_ptr,
            &u_ptr,
            &r_ptr,
            &qfm,
            &rfm,
            &pool]()
            {
                groupMatchByQueryRef(mv3_ptr,
                    u_ptr,
                    r_ptr,
                    *rfm,
                    *qfm,
                    pool);          // 仍用同一池
            });
    }
    pool.waitAllTasksDone();                          // —— Phase-1 完

    // 一步到位，最快释放
    for (auto it = species_match_map->begin();
        it != species_match_map->end(); ++it)
    {
        it->second.reset();   // 释放 MatchVec3D 对象占用的全部内存
    }
    species_match_map->clear();         // 清掉 map 自身节点


    /*========================= Phase-2  : sort ======================*/
    for (auto it = unique_map.begin(); it != unique_map.end(); ++it) {
        MatchByStrandByQueryRefPtr u_ptr = it->second;
        pool.enqueue([&u_ptr, &pool]() {
            sortMatchByQueryStart(u_ptr, pool);
            });
    }
    for (auto it = repeat_map.begin(); it != repeat_map.end(); ++it) {
        MatchByStrandByQueryRefPtr r_ptr = it->second;
        pool.enqueue([&r_ptr, &pool]() {
            sortMatchByQueryStart(r_ptr, pool);
            });
    }
    pool.waitAllTasksDone();                          // —— Phase-2 完

    /*========================= Phase-3  : cluster ===================*/
    using Fut = std::future<ClusterVecPtrByStrandByQueryRefPtr>;
    std::unordered_map<SpeciesName, Fut> fut_map;

    for (auto it = unique_map.begin(); it != unique_map.end(); ++it) {
        const SpeciesName  species = it->first;
        auto               u_ptr = it->second;
        auto               r_ptr = repeat_map.at(species);

        fut_map.emplace(
            species,
            pool.enqueue(
                [&u_ptr, &r_ptr, &pool]() -> ClusterVecPtrByStrandByQueryRefPtr {
                    return clusterAllChrMatch(u_ptr, r_ptr, pool);
                })
        );
    }

    // pool.waitAllTasksDone();
    // 收集 future
    for (auto it = fut_map.begin(); it != fut_map.end(); ++it) {
        const SpeciesName& species = it->first;
        ClusterVecPtrByStrandByQueryRefPtr clus_ptr = it->second.get();
        cluster_map.emplace(species, std::move(clus_ptr));
    }
    pool.waitAllTasksDone();                          // —— Phase-3 完

    return;

}
