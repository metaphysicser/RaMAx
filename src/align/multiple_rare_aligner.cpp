#include <sequence_utils.h>

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
    SeqPro::Length sampling_interval)
{
    std::vector leaf_vec = newick_tree.orderLeavesGreedyMinSum(tree_root);
	uint_t leaf_num = leaf_vec.size();
    // 初始化Ref缓存
    sdsl::int_vector<0> ref_global_cache;
    // TODO 完成迭代的星比对
    for (uint_t i = 0; i < 1; i++) {
    /*for (uint_t i = 0; i < leaf_num; i++) {*/
        // 使用工具函数构建缓存
        SequenceUtils::buildRefGlobalCache(seqpro_managers[newick_tree.getNodes()[leaf_vec[i]].name], sampling_interval, ref_global_cache);

		uint_t ref_id = leaf_vec[i];
		SpeciesName ref_name = newick_tree.getNodes()[ref_id].name;
        std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant> species_fasta_manager_map;
        for (uint_t j = i; j < leaf_num; j++) {
            uint_t query_id = leaf_vec[j];
			SpeciesName query_name = newick_tree.getNodes()[query_id].name;
            auto query_fasta_manager = seqpro_managers.at(query_name);
            species_fasta_manager_map.emplace(query_name, query_fasta_manager);
        }

        // 创建共享线程池，供比对和过滤过程共同使用
        ThreadPool shared_pool(thread_num);
        
        SpeciesMatchVec3DPtrMapPtr match_ptr = alignMultipleGenome(
            ref_name, species_fasta_manager_map,
            ACCURATE_SEARCH, fast_build, allow_MEM, ref_global_cache, sampling_interval
        );
        
        // 使用同一个线程池进行过滤比对结果，获取cluster数据
		SpeciesClusterMapPtr cluster_map = filterMultipeSpeciesAnchors(
			ref_name, species_fasta_manager_map, match_ptr, shared_pool
		);
        
        // 创建多基因组图
        RaMesh::RaMeshMultiGenomeGraph multi_graph;
        
        // 并行构建多个比对结果图，共用线程池
        constructMultipleGraphsByGreedy(
            ref_name, *cluster_map, multi_graph, shared_pool
        );
        
        spdlog::info("Multiple genome graphs constructed for reference: {}", ref_name);

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


SpeciesClusterMapPtr MultipleRareAligner::filterMultipeSpeciesAnchors(
    SpeciesName                       ref_name,
    std::unordered_map<SpeciesName, SeqPro::SharedManagerVariant>& species_fm_map,
    SpeciesMatchVec3DPtrMapPtr        species_match_map,
    ThreadPool&                       shared_pool)
{
    if (!species_match_map || species_match_map->empty()) {
        return std::make_shared<SpeciesClusterMap>();
    }

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
        shared_pool.enqueue(
            [mv3_ptr,
            u_ptr,
            r_ptr,
            qfm,
            rfm,
            &shared_pool]() mutable
            {
                groupMatchByQueryRef(mv3_ptr,
                    u_ptr,
                    r_ptr,
                    *rfm,
                    *qfm,
                    shared_pool);          // 仍用同一池
            });
    }
    shared_pool.waitAllTasksDone();                          // —— Phase-1 完

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
        shared_pool.enqueue([u_ptr, &shared_pool]() mutable {
            sortMatchByQueryStart(u_ptr, shared_pool);
            });
    }
    for (auto it = repeat_map.begin(); it != repeat_map.end(); ++it) {
        MatchByStrandByQueryRefPtr r_ptr = it->second;
        shared_pool.enqueue([r_ptr, &shared_pool]() mutable {
            sortMatchByQueryStart(r_ptr, shared_pool);
            });
    }
    shared_pool.waitAllTasksDone();                          // —— Phase-2 完

    /*========================= Phase-3  : cluster ===================*/
    using Fut = std::future<ClusterVecPtrByStrandByQueryRefPtr>;
    std::unordered_map<SpeciesName, Fut> fut_map;

    for (auto it = unique_map.begin(); it != unique_map.end(); ++it) {
        const SpeciesName  species = it->first;
        auto               u_ptr = it->second;
        auto               r_ptr = repeat_map.at(species);

        fut_map.emplace(
            species,
            shared_pool.enqueue(
                [u_ptr, r_ptr, &shared_pool]() mutable -> ClusterVecPtrByStrandByQueryRefPtr {
                    return clusterAllChrMatch(u_ptr, r_ptr, shared_pool);
                })
        );
    }

    // shared_pool.waitAllTasksDone();
    // 收集 future
    for (auto it = fut_map.begin(); it != fut_map.end(); ++it) {
        const SpeciesName& species = it->first;
        ClusterVecPtrByStrandByQueryRefPtr clus_ptr = it->second.get();
        cluster_map.emplace(species, std::move(clus_ptr));
    }
    shared_pool.waitAllTasksDone();                          // —— Phase-3 完

    // 返回cluster数据用于后续处理
    auto cluster_map_ptr = std::make_shared<SpeciesClusterMap>(std::move(cluster_map));
    return cluster_map_ptr;

}

/* ============================================================= *
 *  多基因组比对类的成员函数：并行构建多个比对结果图，共用一个线程池
 * ============================================================= */
/**
 * @brief  并行构建多个比对结果图，基于贪婪算法处理多个物种的cluster数据
 * 
 * @param ref_name             参考物种名称
 * @param species_cluster_map  所有物种的cluster数据映射
 * @param graph                多基因组图对象
 * @param shared_pool          共享线程池
 * @param min_span             最小跨度阈值
 */
void MultipleRareAligner::constructMultipleGraphsByGreedy(
    SpeciesName ref_name,
    const SpeciesClusterMap& species_cluster_map,
    RaMesh::RaMeshMultiGenomeGraph& graph,
    ThreadPool& shared_pool,
    uint_t min_span)
{
    if (species_cluster_map.empty()) {
        spdlog::warn("[constructMultipleGraphsByGreedy] Empty cluster map, nothing to process.");
        return;
    }

    spdlog::info("[constructMultipleGraphsByGreedy] Processing {} species clusters", 
                species_cluster_map.size());

    /* ---------- 1. 为每个物种并行处理cluster数据 ---------- */
    std::vector<std::future<void>> species_futures;
    species_futures.reserve(species_cluster_map.size());

    for (const auto& [species_name, cluster_ptr] : species_cluster_map) {
        if (species_name == ref_name) continue;  // 跳过参考物种

        // 为每个物种启动异步任务
        auto species_future = std::async(std::launch::async, 
            [this, ref_name, species_name, cluster_ptr, &graph, &shared_pool, min_span]() {
                try {
                    if (!cluster_ptr || cluster_ptr->empty()) {
                        spdlog::warn("[constructMultipleGraphsByGreedy] Empty cluster data for species: {}", 
                                   species_name);
                        return;
                    }

                    // 使用PairRareAligner的贪婪算法构建图
                    PairRareAligner pra(*this);
                    pra.ref_name = ref_name;

                    // 为每个chromosome的cluster数据并行处理
                    std::vector<std::future<void>> chr_futures;
                    for (const auto& strand_data : *cluster_ptr) {
                        for (const auto& query_ref_data : strand_data) {
                            for (const auto& cluster_vec : query_ref_data) {
                                if (cluster_vec && !cluster_vec->empty()) {
                                    // 将集合转换为向量以便处理
                                    auto cluster_vec_ptr = std::make_shared<MatchClusterVec>(
                                        cluster_vec->begin(), cluster_vec->end());
                                    
                                    chr_futures.emplace_back(
                                        std::async(std::launch::async, 
                                            [&pra, species_name, cluster_vec_ptr, &graph, min_span]() {
                                                pra.constructGraphByGreedy(species_name, cluster_vec_ptr, 
                                                                         graph, min_span);
                                            })
                                    );
                                }
                            }
                        }
                    }

                    // 等待所有chromosome处理完成
                    for (auto& chr_future : chr_futures) {
                        chr_future.wait();
                    }

                    spdlog::info("[constructMultipleGraphsByGreedy] Species {} processed successfully", 
                               species_name);
                }
                catch (const std::exception& e) {
                    spdlog::error("[constructMultipleGraphsByGreedy] Error processing species {}: {}", 
                                species_name, e.what());
                }
            });

        species_futures.emplace_back(std::move(species_future));
    }

    /* ---------- 2. 等待所有物种处理完成 ---------- */
    for (auto& future : species_futures) {
        try {
            future.wait();
        }
        catch (const std::exception& e) {
            spdlog::error("[constructMultipleGraphsByGreedy] Error waiting for species processing: {}", 
                        e.what());
        }
    }

    // 确保共享线程池中的所有任务都完成
    shared_pool.waitAllTasksDone();

    spdlog::info("[constructMultipleGraphsByGreedy] All species graphs constructed successfully");
}
