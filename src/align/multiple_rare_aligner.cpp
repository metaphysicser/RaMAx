#include <sequence_utils.h>
#include <algorithm>

#include "rare_aligner.h"
#include "anchor.h"  // 包含 UnionFind 定义
#include "SeqPro.h"  // 包含 SeqPro 相关定义
#include "ramesh.h"  // 包含 RaMesh 图结构定义

// 辅助函数：根据CIGAR字符串计算query区间对应关系
namespace {

    /**
     * @brief 获取query segment在ref上的映射区间
     * @param query_segment query segment对象
     * @param ref_block_ptr segment所属的block，用于获取ref anchor
     * @param ref_name 参考物种名称
     * @return pair<ref_start, ref_end> 映射到ref上的区间，如果无法映射返回{0,0}
     */
    std::pair<uint_t, uint_t> getRefMappedInterval(
        const RaMesh::SegPtr& query_segment,
        const RaMesh::BlockPtr& ref_block_ptr,
        const SpeciesName& ref_name) {

        if (!query_segment || !ref_block_ptr) {
            return {0, 0};
        }

        // 查找对应的ref segment
        std::shared_lock block_lock(ref_block_ptr->rw);

        // 在block的anchors中查找ref segment
        for (const auto& [species_chr_pair, segment] : ref_block_ptr->anchors) {
            const auto& [species_name, chr_name] = species_chr_pair;
            if (species_name == ref_name && segment && segment->isSegment()) {
                // 找到了ref segment，返回其区间
                return {segment->start, segment->start + segment->length};
            }
        }

        spdlog::warn("[getRefMappedInterval] Cannot find corresponding ref segment");
        return {0, 0};
    }

    /**
     * @brief 根据CIGAR字符串计算ref区间对应的query区间
     * @param cigar 原始segment的CIGAR字符串
     * @param target_ref_start 目标ref区间的起始位置
     * @param target_ref_end 目标ref区间的结束位置
     * @param original_ref_start 原始segment在ref上的起始位置
     * @param original_query_start 原始segment在query上的起始位置
     * @return pair<query_start, query_length> 对应的query区间
     */
    std::pair<uint_t, uint_t> calculateQueryInterval(
        const Cigar_t& cigar,
        uint_t target_ref_start, uint_t target_ref_end,
        uint_t original_ref_start, uint_t original_query_start) {

        // 如果CIGAR为空，使用简单的线性映射
        if (cigar.empty()) {
            uint_t ref_offset = target_ref_start - original_ref_start;
            uint_t length = target_ref_end - target_ref_start;
            return {original_query_start + ref_offset, length};
        }

        uint_t current_ref_pos = original_ref_start;
        uint_t current_query_pos = original_query_start;
        uint_t target_query_start = 0;
        uint_t target_query_end = 0;
        bool found_start = false;
        bool found_end = false;

        // 遍历CIGAR操作
        for (const auto& cigar_unit : cigar) {
            uint_t length = cigar_unit >> 4;  // 高28位是长度
            uint_t op_code = cigar_unit & 0xF; // 低4位是操作码

            switch (op_code) {
                case 0x0: // M - match/mismatch
                case 0x7: // = - exact match
                case 0x8: // X - mismatch
                {
                    // 这些操作同时消耗ref和query位置
                    uint_t ref_end_this_op = current_ref_pos + length;

                    // 检查target_ref_start是否在这个操作范围内
                    if (!found_start && current_ref_pos <= target_ref_start && target_ref_start < ref_end_this_op) {
                        uint_t offset = target_ref_start - current_ref_pos;
                        target_query_start = current_query_pos + offset;
                        found_start = true;
                    }

                    // 检查target_ref_end是否在这个操作范围内
                    if (!found_end && current_ref_pos <= target_ref_end && target_ref_end <= ref_end_this_op) {
                        uint_t offset = target_ref_end - current_ref_pos;
                        target_query_end = current_query_pos + offset;
                        found_end = true;
                        break;
                    }

                    current_ref_pos += length;
                    current_query_pos += length;
                    break;
                }
                case 0x1: // I - insertion (query only)
                {
                    // 插入操作只消耗query位置，不影响ref坐标映射
                    current_query_pos += length;
                    break;
                }
                case 0x2: // D - deletion (ref only)
                {
                    // 删除操作只消耗ref位置
                    uint_t ref_end_this_op = current_ref_pos + length;

                    // 检查删除区间是否包含我们的目标区间
                    if (!found_start && current_ref_pos <= target_ref_start && target_ref_start < ref_end_this_op) {
                        target_query_start = current_query_pos;
                        found_start = true;
                    }

                    if (!found_end && current_ref_pos <= target_ref_end && target_ref_end <= ref_end_this_op) {
                        target_query_end = current_query_pos;
                        found_end = true;
                        break;
                    }

                    current_ref_pos += length;
                    // query位置不变
                    break;
                }
                default:
                    // 其他CIGAR操作暂时按match处理
                    current_ref_pos += length;
                    current_query_pos += length;
                    break;
            }

            if (found_start && found_end) break;
        }

        // 如果没有找到精确的映射，使用线性近似
        if (!found_start || !found_end) {
            spdlog::warn("[calculateQueryInterval] Cannot precisely map CIGAR interval, using linear approximation");
            uint_t ref_offset = target_ref_start - original_ref_start;
            uint_t length = target_ref_end - target_ref_start;
            return {original_query_start + ref_offset, length};
        }

        uint_t query_length = target_query_end - target_query_start;
        return {target_query_start, query_length};
    }

    /**
     * @brief 创建分割后的segment CIGAR字符串
     * @param original_cigar 原始CIGAR
     * @param ref_start_offset ref区间在原始segment中的偏移
     * @param ref_length ref区间长度
     * @return 分割后的CIGAR字符串
     */
    Cigar_t extractCigarSubsequence(
        const Cigar_t& original_cigar,
        uint_t ref_start_offset, uint_t ref_length) {

        if (original_cigar.empty()) {
            // 如果原始CIGAR为空，创建一个简单的match操作
            return Cigar_t{cigarToInt('M', ref_length)};
        }

        Cigar_t result_cigar;
        uint_t current_ref_pos = 0;
        uint_t target_ref_start = ref_start_offset;
        uint_t target_ref_end = ref_start_offset + ref_length;

        for (const auto& cigar_unit : original_cigar) {
            uint_t length = cigar_unit >> 4;
            uint_t op_code = cigar_unit & 0xF;

            switch (op_code) {
                case 0x0: case 0x7: case 0x8: case 0x2: // M, =, X, D - 消耗ref位置
                {
                    uint_t ref_end_this_op = current_ref_pos + length;

                    // 检查这个操作是否与目标区间重叠
                    if (current_ref_pos < target_ref_end && ref_end_this_op > target_ref_start) {
                        uint_t overlap_start = std::max(current_ref_pos, target_ref_start);
                        uint_t overlap_end = std::min(ref_end_this_op, target_ref_end);
                        uint_t overlap_length = overlap_end - overlap_start;

                        if (overlap_length > 0) {
                            result_cigar.push_back(cigarToInt(
                                (op_code == 0x0) ? 'M' :
                                (op_code == 0x7) ? '=' :
                                (op_code == 0x8) ? 'X' : 'D',
                                overlap_length));
                        }
                    }

                    current_ref_pos += length;
                    break;
                }
                case 0x1: // I - 插入操作，不消耗ref位置
                {
                    // 插入操作在ref坐标范围内时保留
                    if (current_ref_pos >= target_ref_start && current_ref_pos < target_ref_end) {
                        result_cigar.push_back(cigarToInt('I', length));
                    }
                    break;
                }
                default:
                    // 其他操作暂时跳过
                    break;
            }
        }

        // 如果结果为空，创建一个默认的match操作
        if (result_cigar.empty()) {
            result_cigar.push_back(cigarToInt('M', ref_length));
        }

        return result_cigar;
    }

} // anonymous namespace

/**
 * @brief 确保 SeqPro manager 是 MaskedSequenceManager 类型
 * @param manager_variant 当前的 manager variant
 * @return 指向 MaskedSequenceManager 的指针
 */
SeqPro::MaskedSequenceManager* ensureMaskedManager(SeqPro::SharedManagerVariant& manager_variant) {
    auto& variant = *manager_variant;
    
    // 检查当前类型
    if (std::holds_alternative<std::unique_ptr<SeqPro::MaskedSequenceManager>>(variant)) {
        // 已经是 MaskedSequenceManager
        return std::get<std::unique_ptr<SeqPro::MaskedSequenceManager>>(variant).get();
    } 
    else if (std::holds_alternative<std::unique_ptr<SeqPro::SequenceManager>>(variant)) {
        // 需要转换为 MaskedSequenceManager
        auto seq_manager = std::move(std::get<std::unique_ptr<SeqPro::SequenceManager>>(variant));
        auto masked_manager = std::make_unique<SeqPro::MaskedSequenceManager>(std::move(seq_manager));
        auto* result_ptr = masked_manager.get();
        
        // 替换 variant 中的内容
        variant = std::move(masked_manager);
        
        return result_ptr;
    }
    
    throw std::runtime_error("Invalid SeqPro manager variant type");
}

/**
 * @brief 从比对结果图中提取已比对的区间，并作为遮蔽区间添加到对应的 SeqPro manager 中
 * @param graph 比对结果图
 * @param seqpro_managers SeqPro manager 映射
 * @param ref_name 参考物种名称
 */
void addAlignedRegionsAsMask(
    const RaMesh::RaMeshMultiGenomeGraph& graph,
    std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    const SpeciesName& ref_name) {
    
    if (graph.blocks.empty()) {
        spdlog::info("[addAlignedRegionsAsMask] No blocks to process for masking");
        return;
    }
    
    spdlog::info("[addAlignedRegionsAsMask] Extracting aligned regions as mask intervals from {} blocks", 
                 graph.blocks.size());
    
    // 按物种和染色体分组收集区间
    std::unordered_map<SpeciesName, std::unordered_map<ChrName, std::vector<SeqPro::MaskInterval>>> 
        species_chr_intervals;
    
    size_t total_intervals = 0;
    size_t valid_blocks = 0;
    
    // 遍历所有 blocks，提取 segment 区间
    for (const auto& weak_block : graph.blocks) {
        auto block_ptr = weak_block.lock();
        if (!block_ptr) continue;
        
        valid_blocks++;
        std::shared_lock block_lock(block_ptr->rw);
        
        // 处理该 block 中的所有 anchors
        for (const auto& [species_chr_pair, segment] : block_ptr->anchors) {
            const auto& [species_name, chr_name] = species_chr_pair;
            
            // 只处理有效的 segment
            if (!segment || !segment->isSegment() || segment->length == 0) {
                continue;
            }
            
            // 检查该物种是否在 seqpro_managers 中
            if (seqpro_managers.find(species_name) == seqpro_managers.end()) {
                continue;
            }
            
            // 创建遮蔽区间（使用原始坐标）
            SeqPro::MaskInterval interval(segment->start, segment->start + segment->length);
            species_chr_intervals[species_name][chr_name].push_back(interval);
            total_intervals++;
        }
    }
    
    spdlog::info("[addAlignedRegionsAsMask] Collected {} intervals from {} valid blocks across {} species", 
                 total_intervals, valid_blocks, species_chr_intervals.size());
    
    // 为每个物种批量添加遮蔽区间
    for (auto& [species_name, chr_intervals] : species_chr_intervals) {
        try {
            // 确保该物种的 manager 是 MaskedSequenceManager
            auto* masked_manager = ensureMaskedManager(seqpro_managers[species_name]);
            
            size_t species_total_intervals = 0;
            
            // 按染色体处理区间
            for (auto& [chr_name, intervals] : chr_intervals) {
                if (intervals.empty()) continue;
                
                // 构造序列名（假设格式为染色体名）
                std::string seq_name = chr_name;
                
                                // 检查序列是否存在
                if (masked_manager->getSequenceId(seq_name) == SeqPro::SequenceIndex::INVALID_ID) {
                    spdlog::warn("[addAlignedRegionsAsMask] Sequence not found: {}:{}, skipping", 
                                species_name, seq_name);
                    continue;
                }
                
                // 批量添加区间（segment中的坐标是遮蔽后的坐标，需要转换为原始坐标）
                masked_manager->addMaskIntervals(seq_name, intervals);
                species_total_intervals += intervals.size();
                
                spdlog::debug("[addAlignedRegionsAsMask] Added {} intervals for {}:{}", 
                             intervals.size(), species_name, seq_name);
            }
            
            // 定案该物种的所有遮蔽区间
            masked_manager->finalizeMaskIntervals();
            
            spdlog::info("[addAlignedRegionsAsMask] Successfully added {} mask intervals for species {}", 
                        species_total_intervals, species_name);
        }
        catch (const std::exception& e) {
            spdlog::error("[addAlignedRegionsAsMask] Error processing species {}: {}", 
                         species_name, e.what());
        }
    }
    
    spdlog::info("[addAlignedRegionsAsMask] Mask interval addition completed for all species");
}

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

std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> MultipleRareAligner::starAlignment(
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
    uint_t tree_root,
    SearchMode                 search_mode,
    bool                       fast_build,
    bool                       allow_MEM,
    bool                       mask_mode,
    SeqPro::Length sampling_interval,
    uint_t min_span)
{
    std::vector leaf_vec = newick_tree.orderLeavesGreedyMinSum(tree_root);
	uint_t leaf_num = leaf_vec.size();
    // 初始化Ref缓存
    sdsl::int_vector<0> ref_global_cache;
    // 创建共享线程池，供比对和过滤过程共同使用
    ThreadPool shared_pool(thread_num);

    // 创建当前迭代的多基因组图
    auto multi_graph = std::make_unique<RaMesh::RaMeshMultiGenomeGraph>();
    for (uint_t i = 0; i < 1; i++) {
    // for (uint_t i = 0; i < leaf_num; i++) {
        // 使用工具函数构建缓存
        spdlog::info("build ref global cache for {}", newick_tree.getNodes()[leaf_vec[i]].name);
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

        spdlog::info("align multiple genome for {}", ref_name);
        SpeciesMatchVec3DPtrMapPtr match_ptr = alignMultipleGenome(
            ref_name, species_fasta_manager_map,
            ACCURATE_SEARCH, fast_build, allow_MEM, ref_global_cache, sampling_interval
        );

        // 使用同一个线程池进行过滤比对结果，获取cluster数据
        spdlog::info("filter multiple species anchors for {}", ref_name);
		SpeciesClusterMapPtr cluster_map = filterMultipeSpeciesAnchors(
			ref_name, species_fasta_manager_map, match_ptr, shared_pool
		);

        // 并行构建多个比对结果图，共用线程池
        spdlog::info("construct multiple genome graphs for {}", ref_name);
        constructMultipleGraphsByGreedy(
           seqpro_managers, ref_name, *cluster_map, *multi_graph, shared_pool, min_span
        );
#ifdef _DEBUG_
        multi_graph->verifyGraphCorrectness(true);
#endif // _DEBUG_

        spdlog::info("merge multiple genome graphs for {}", ref_name);


        mergeMultipleGraphs(ref_name, *multi_graph, shared_pool);

#ifdef _DEBUG_
        multi_graph->verifyGraphCorrectness(true);
#endif // _DEBUG_

        // 将当前轮次的比对结果作为遮蔽区间添加到 SeqPro managers 中
        // 这样后续轮次就不会重复比对已经成功比对的区间
        spdlog::info("Adding aligned regions as mask intervals for {}", ref_name);
        try {
            addAlignedRegionsAsMask(*multi_graph, seqpro_managers, ref_name);
            spdlog::info("Successfully added mask intervals for round with reference {}", ref_name);
        }
        catch (const std::exception& e) {
            spdlog::error("Failed to add mask intervals for {}: {}", ref_name, e.what());
        }
    }

    // 在所有迭代完成后，进行最终的图正确性验证
    spdlog::default_logger()->flush();


    return std::move(multi_graph);

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
    round_id++;
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

    /*========================= Phase-3  :  ===================*/
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
std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers,
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

    // 【修复】：添加互斥锁保护graph的并发访问
    std::mutex graph_mutex;

    /* ---------- 1. 为每个物种并行处理cluster数据 ---------- */
    std::vector<std::future<void>> species_futures;
    species_futures.reserve(species_cluster_map.size());

    for (const auto& [species_name, cluster_ptr] : species_cluster_map) {
        if (species_name == ref_name) continue;  // 跳过参考物种

        // 为每个物种启动异步任务
        auto species_future = std::async(std::launch::async,
            [this, ref_name, species_name, cluster_ptr, &graph, &shared_pool, &graph_mutex, min_span,&seqpro_managers]() {
                try {
                    if (!cluster_ptr || cluster_ptr->empty()) {
                        spdlog::warn("[constructMultipleGraphsByGreedy] Empty cluster data for species: {}",
                                   species_name);
                        return;
                    }

                    // 使用PairRareAligner的贪婪算法构建图
                    PairRareAligner pra(*this);
                    pra.ref_name = ref_name;
                    // 【修复】：设置ref_seqpro_manager，避免空指针
                    pra.ref_seqpro_manager = &(*seqpro_managers.at(ref_name));

                    //// 【修复】：避免过度并行化，改为串行处理chromosome数据
                    //// 每个物种内部串行处理，避免graph的深度并发访问
                    //for (const auto& strand_data : *cluster_ptr) {
                    //    for (const auto& query_ref_data : strand_data) {
                    //        for (const auto& cluster_vec : query_ref_data) {
                    //            if (cluster_vec && !cluster_vec->empty()) {
                    //                // 将集合转换为向量以便处理
                    //                auto cluster_vec_ptr = std::make_shared<MatchClusterVec>(
                    //                    cluster_vec->begin(), cluster_vec->end());

                    //                // 【修复】：使用互斥锁保护graph访问
                    //                {
                    //                    std::lock_guard<std::mutex> lock(graph_mutex);
                    //                    pra.constructGraphByGreedy(species_name, *seqpro_managers[species_name],cluster_vec_ptr,
                    //                                             graph, min_span);
                    //                }
                    //            }
                    //        }
                    //    }
                    //}
                    {
                        std::lock_guard<std::mutex> lock(graph_mutex);
                        pra.constructGraphByGreedy(species_name, *seqpro_managers[species_name], cluster_ptr,
                                                                      graph, min_span);
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


void MultipleRareAligner::mergeMultipleGraphs(
    const SpeciesName& ref_name,
    RaMesh::RaMeshMultiGenomeGraph& graph,
    ThreadPool& shared_pool)
{
#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Starting multi-genome graph merging, reference species: {}, total {} blocks",
                ref_name, graph.blocks.size());
#endif

    if (graph.blocks.empty()) {
#ifdef _DEBUG_
        spdlog::warn("[mergeMultipleGraphs] No blocks to merge");
#endif
        return;
    }

    // 性能计时器
    auto start_time = std::chrono::high_resolution_clock::now();

    // ==================== 阶段一：识别重叠的Block关联集 ====================
#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Phase 1: Identifying overlapping block association sets...");
#endif
    auto phase1_start = std::chrono::high_resolution_clock::now();

    // 轻量级区间信息结构
    struct SegmentInterval {
        RaMesh::SegPtr segment;
        uint_t start;
        uint_t end;
        size_t block_index;
        RaMesh::WeakBlock weak_block;  // 保存block的弱引用

        SegmentInterval(RaMesh::SegPtr seg, size_t idx, RaMesh::WeakBlock wb)
            : segment(seg), start(seg->start), end(seg->start + seg->length),
              block_index(idx), weak_block(wb) {}
    };

    // 按染色体分组存储ref-segments
    using SpeciesChrKey = std::pair<SpeciesName, ChrName>;
    using ChrSegmentMap = std::unordered_map<SpeciesChrKey, std::vector<SegmentInterval>,
                                           RaMesh::SpeciesChrPairHash>;

    // 第1步：收集所有ref-segments，并预留容量
    ChrSegmentMap chr_segments;
    chr_segments.reserve(graph.blocks.size() / 10); // 预估染色体数量

    for (size_t i = 0; i < graph.blocks.size(); ++i) {
        auto block_ptr = graph.blocks[i].lock();
        if (!block_ptr) continue;

        std::shared_lock block_lock(block_ptr->rw);
        for (const auto& [species_chr_pair, segment] : block_ptr->anchors) {
            const auto& [species_name, chr_name] = species_chr_pair;
            if (segment && segment->isSegment() && species_name == ref_name) {
                SpeciesChrKey key = species_chr_pair;
                chr_segments[key].emplace_back(segment, i, graph.blocks[i]);
            }
        }
    }

    // 统计总segment数量
    size_t total_segments = 0;
    for (const auto& [key, segments] : chr_segments) {
        total_segments += segments.size();
    }

#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Step 1 completed: Extracted {} ref-segments from {} blocks, distributed across {} chromosomes",
                total_segments, graph.blocks.size(), chr_segments.size());
#endif

    // 第二步：并行排序每个染色体的segments
    if (!chr_segments.empty()) {
        std::vector<std::future<void>> sort_futures;
        sort_futures.reserve(chr_segments.size());

        for (auto& [key, segments] : chr_segments) {
            if (segments.size() <= 1) continue;

            // 只对较大的染色体进行并行排序
            if (segments.size() > 1000) {
                auto sort_future = shared_pool.enqueue([&segments]() {
                    std::sort(segments.begin(), segments.end(),
                             [](const SegmentInterval& a, const SegmentInterval& b) {
                                 if (a.start != b.start) return a.start < b.start;
                                 return a.end < b.end;
                             });
                });
                sort_futures.emplace_back(std::move(sort_future));
            } else {
                // 小染色体直接排序
                std::sort(segments.begin(), segments.end(),
                         [](const SegmentInterval& a, const SegmentInterval& b) {
                             if (a.start != b.start) return a.start < b.start;
                             return a.end < b.end;
                         });
            }
        }

        for (auto& future : sort_futures) {
            future.get();
        }
    }

    // 第三步：构建重叠关系图
#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Step 3: Building overlap relationship graph...");
#endif

    // 使用并查集来跟踪Block之间的连通关系
    UnionFind block_union_find(graph.blocks.size());

    // 统计重叠检测的计数
    std::atomic<size_t> overlap_count(0);
    std::atomic<size_t> total_pairs_checked(0);

    // 【修复】：收集所有重叠关系，避免并发unite操作
    std::mutex overlaps_mutex;
    std::vector<std::pair<size_t, size_t>> all_overlaps;

    // 并行检测每个染色体内的重叠关系
    std::vector<std::future<void>> overlap_futures;
    overlap_futures.reserve(chr_segments.size());

    for (const auto& [chr_key, segments] : chr_segments) {
        if (segments.size() <= 1) continue;

        // 只对有足够segments的染色体进行并行处理
        if (segments.size() > 100) {
            auto overlap_future = shared_pool.enqueue([&segments,
                                                      &overlap_count, &total_pairs_checked,
                                                      &all_overlaps, &overlaps_mutex]() {
                size_t local_overlap_count = 0;
                size_t local_pairs_checked = 0;
                std::vector<std::pair<size_t, size_t>> local_overlaps;
                local_overlaps.reserve(segments.size() * 2); // 预分配

                // 扫描线算法优化重叠检测
                for (size_t i = 0; i < segments.size(); ++i) {
                    const auto& seg_i = segments[i];

                    for (size_t j = i + 1; j < segments.size(); ++j) {
                        const auto& seg_j = segments[j];

                        if (seg_j.start >= seg_i.end) {
                            break; // 后续不可能重叠
                        }

                        local_pairs_checked++;

                        if (seg_i.start < seg_j.end && seg_j.start < seg_i.end) {
                            local_overlaps.emplace_back(seg_i.block_index, seg_j.block_index);
                            local_overlap_count++;
                        }
                    }
                }

                // 【修复】：批量添加到全局重叠列表而不是直接unite
                if (!local_overlaps.empty()) {
                    std::lock_guard<std::mutex> lock(overlaps_mutex);
                    all_overlaps.insert(all_overlaps.end(),
                                      local_overlaps.begin(), local_overlaps.end());
                }

                overlap_count += local_overlap_count;
                total_pairs_checked += local_pairs_checked;
            });

            overlap_futures.emplace_back(std::move(overlap_future));
        } else {
            // 小染色体直接处理
            for (size_t i = 0; i < segments.size(); ++i) {
                const auto& seg_i = segments[i];

                for (size_t j = i + 1; j < segments.size(); ++j) {
                    const auto& seg_j = segments[j];

                    if (seg_j.start >= seg_i.end) break;

                    total_pairs_checked++;

                    if (seg_i.start < seg_j.end && seg_j.start < seg_i.end) {
                        // 【修复】：直接添加到重叠列表
                        std::lock_guard<std::mutex> lock(overlaps_mutex);
                        all_overlaps.emplace_back(seg_i.block_index, seg_j.block_index);
                        overlap_count++;
                    }
                }
            }
        }
    }

    for (auto& future : overlap_futures) {
        future.get();
    }

    // 串行执行所有unite操作，确保线程安全
#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Starting serial execution of {} unite operations...", all_overlaps.size());
#endif
    for (const auto& [block_i, block_j] : all_overlaps) {
        block_union_find.unite(block_i, block_j);
    }

    auto phase1_end = std::chrono::high_resolution_clock::now();
    auto phase1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(phase1_end - phase1_start);

#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Phase 1 completed ({} ms): Checked {} segment pairs, found {} overlap relationships",
                phase1_duration.count(), total_pairs_checked.load(), overlap_count.load());
#endif

    // 第四步：划分连通分量
    std::unordered_map<int_t, std::vector<size_t>> components;
    components.reserve(graph.blocks.size() / 2); // 预估连通分量数量

    for (size_t i = 0; i < graph.blocks.size(); ++i) {
        int_t root = block_union_find.find(i);
        components[root].push_back(i);
    }

    // 统计连通分量信息
    size_t num_components = components.size();
    size_t largest_component_size = 0;
    size_t singleton_components = 0;

    for (const auto& [root, block_indices] : components) {
        if (block_indices.size() == 1) {
            singleton_components++;
        }
        largest_component_size = std::max(largest_component_size, block_indices.size());
    }

#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Connected component statistics:");
    spdlog::info("[mergeMultipleGraphs]   - Total blocks: {}", graph.blocks.size());
    spdlog::info("[mergeMultipleGraphs]   - Number of components: {}", num_components);
    spdlog::info("[mergeMultipleGraphs]   - Singleton blocks: {}", singleton_components);
    spdlog::info("[mergeMultipleGraphs]   - Largest component size: {}", largest_component_size);
#endif

    // ==================== 阶段二：并行处理每个Block关联集 ====================
#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Phase 2: Starting parallel processing of {} block association sets...", num_components);
#endif
    auto phase2_start = std::chrono::high_resolution_clock::now();

    // 存储重建后的新blocks
    std::mutex new_blocks_mutex;
    std::vector<RaMesh::BlockPtr> new_blocks;
    new_blocks.reserve(total_segments); // 预分配，避免重分配

    // 统计处理进度
    std::atomic<size_t> processed_components(0);
    std::atomic<size_t> created_blocks(0);

    // 为每个连通分量启动并行任务
    std::vector<std::future<void>> component_futures;

    // 根据连通分量大小决定处理策略
    std::vector<std::pair<size_t, std::vector<size_t>>> sorted_components;
    for (const auto& [root, indices] : components) {
        sorted_components.emplace_back(indices.size(), indices); // 不使用std::move，避免破坏原始数据
    }

    // 按大小排序，大的先处理
    std::sort(sorted_components.rbegin(), sorted_components.rend());

    size_t component_id = 0;
    for (const auto& [comp_size, block_indices] : sorted_components) {
        component_id++;

        // 单独的Block直接处理
        if (comp_size == 1) {
            auto block_ptr = graph.blocks[block_indices[0]].lock();
            if (block_ptr) {
                std::lock_guard<std::mutex> lock(new_blocks_mutex);
                new_blocks.push_back(block_ptr);
                processed_components++;
            }
            continue;
        }

        // 大连通分量并行处理
        auto component_future = shared_pool.enqueue(
            [component_id, &block_indices, &graph, &chr_segments, ref_name,
             &new_blocks, &new_blocks_mutex, &processed_components, &created_blocks, &sorted_components]() {
            try {
                // spdlog::debug("[mergeMultipleGraphs] 处理连通分量 {}: {} 个blocks",
                //              component_id, block_indices.size());

                // 第4步：计算基本区间
                std::set<uint_t> coordinates;
                std::vector<RaMesh::BlockPtr> valid_blocks;
                valid_blocks.reserve(block_indices.size());

                // 收集该连通分量中所有有效blocks的ref-segment坐标
                for (size_t block_idx : block_indices) {
                    auto block_ptr = graph.blocks[block_idx].lock();
                    if (!block_ptr) continue;

                    valid_blocks.push_back(block_ptr);

                    std::shared_lock block_lock(block_ptr->rw);
                    for (const auto& [species_chr_pair, segment] : block_ptr->anchors) {
                        const auto& [species_name, chr_name] = species_chr_pair;
                        if (segment && segment->isSegment() && species_name == ref_name) {
                            coordinates.insert(segment->start);
                            coordinates.insert(segment->start + segment->length);
                        }
                    }
                }

                if (coordinates.size() < 2) {
#ifdef _DEBUG_
                    spdlog::debug("[mergeMultipleGraphs] Component {} has insufficient coordinate points, skipping", component_id);
#endif
                    processed_components++;
                    return;
                }

                // 将坐标转换为有序向量
                std::vector<uint_t> sorted_coords(coordinates.begin(), coordinates.end());

                // 第5步：批量处理基本区间
                std::vector<RaMesh::BlockPtr> component_new_blocks;
                component_new_blocks.reserve(sorted_coords.size() - 1);

                // 避免过度并行化，大连通分量也使用串行处理以防止线程争用
                // 所有连通分量都使用串行处理基本区间，避免在并行任务内部再创建并行任务
                for (size_t i = 0; i < sorted_coords.size() - 1; ++i) {
                    uint_t interval_start = sorted_coords[i];
                    uint_t interval_end = sorted_coords[i + 1];

                    auto new_block = MultipleRareAligner::processElementaryInterval(interval_start, interval_end,
                                                             valid_blocks, ref_name);
                    if (new_block && !new_block->anchors.empty()) {
                        component_new_blocks.push_back(new_block);
                    }
                }

                // 批量添加到全局结果
                if (!component_new_blocks.empty()) {
                    std::lock_guard<std::mutex> lock(new_blocks_mutex);
                    new_blocks.insert(new_blocks.end(),
                                    component_new_blocks.begin(),
                                    component_new_blocks.end());
                    created_blocks += component_new_blocks.size();
                }

                processed_components++;

                // 定期报告进度
#ifdef _DEBUG_
                if (processed_components % 1000 == 0) {
                    spdlog::debug("[mergeMultipleGraphs] Progress: {} / {} components processed",
                               processed_components.load(), sorted_components.size());
                }
#endif
            }
            catch (const std::exception& e) {
                spdlog::error("[mergeMultipleGraphs] Error processing component {}: {}",
                            component_id, e.what());
                processed_components++;
            }
        });

        component_futures.emplace_back(std::move(component_future));
    }

    // 等待所有连通分量处理完成
    for (auto& future : component_futures) {
        future.get();
    }

    auto phase2_end = std::chrono::high_resolution_clock::now();
    auto phase2_duration = std::chrono::duration_cast<std::chrono::milliseconds>(phase2_end - phase2_start);

#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Phase 2 completed ({} ms): Rebuilt {} new blocks",
                phase2_duration.count(), new_blocks.size());
#endif

    // ==================== 阶段三：原子化图结构更新 ====================
#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Phase 3: Starting atomic graph structure update...");
#endif
    auto phase3_start = std::chrono::high_resolution_clock::now();

    // 收集要被替换的原始block的索引
    std::set<size_t> blocks_to_remove;
    for (const auto& [comp_size, block_indices] : sorted_components) {
        if (comp_size > 1) {
            blocks_to_remove.insert(block_indices.begin(), block_indices.end());
        }
    }

    // 保存原始blocks数量，用于最终统计
    size_t original_blocks_count = graph.blocks.size();

    // 在全局写锁保护下执行原子化更新
    size_t processed_segments = 0;
    {
        std::unique_lock<std::shared_mutex> global_lock(graph.rw);

        // 更新Block列表
        auto& global_blocks = graph.blocks;
        size_t removed_blocks = 0;

        // 创建新的blocks列表
        std::vector<RaMesh::WeakBlock> new_global_blocks;
        new_global_blocks.reserve(global_blocks.size() - blocks_to_remove.size() + new_blocks.size());

        // 保留不在移除列表中的block
        for (size_t i = 0; i < global_blocks.size(); ++i) {
            if (blocks_to_remove.count(i) == 0) {
                if (!global_blocks[i].expired()) {
                    new_global_blocks.push_back(global_blocks[i]);
                }
            } else {
                removed_blocks++;
            }
        }

        // 添加所有新生成的Block
        for (auto& new_block : new_blocks) {
            new_global_blocks.emplace_back(RaMesh::WeakBlock(new_block));
        }

        // 替换全局blocks列表
        global_blocks = std::move(new_global_blocks);

#ifdef _DEBUG_
        spdlog::debug("[mergeMultipleGraphs] Removed {} original blocks, final block count: {}",
                    removed_blocks, global_blocks.size());

        // 重建GenomeEnd链表
        spdlog::debug("[mergeMultipleGraphs] Starting chromosome segment chain reconstruction...");
#endif

        // 收集所有有效blocks中的segments
        std::unordered_map<std::pair<SpeciesName, ChrName>, std::vector<RaMesh::SegPtr>,
                          RaMesh::SpeciesChrPairHash> all_segments_by_chr;

        for (const auto& weak_block : global_blocks) {
            auto block = weak_block.lock();
            if (!block) continue;

            std::shared_lock block_lock(block->rw);
            for (const auto& [species_chr_pair, segment] : block->anchors) {
                if (segment && segment->isSegment()) {
                    all_segments_by_chr[species_chr_pair].push_back(segment);
                }
            }
        }

        // 为每个染色体重建完整的segment链表
        for (auto& [species_chr_key, segments] : all_segments_by_chr) {
            const auto& [species_name, chr_name] = species_chr_key;

            auto species_graph_it = graph.species_graphs.find(species_name);
            if (species_graph_it == graph.species_graphs.end()) continue;

            auto& species_graph = species_graph_it->second;
            auto chr_end_it = species_graph.chr2end.find(chr_name);
            if (chr_end_it == species_graph.chr2end.end()) continue;

            auto& genome_end = chr_end_it->second;

            // 排序segments
            std::sort(segments.begin(), segments.end(),
                [](RaMesh::SegPtr a, RaMesh::SegPtr b) {
                    return a->start < b->start;
                });

            // 去重：移除完全相同的segments（相同的起始位置和长度）
            // 先移除length==0的无效segment
            segments.erase(std::remove_if(segments.begin(), segments.end(),
                [](RaMesh::SegPtr s){ return s->length == 0; }), segments.end());

            if (segments.empty()) {
                continue; // 没有有效segment可处理
            }

            auto last = std::unique(segments.begin(), segments.end(),
                [](RaMesh::SegPtr a, RaMesh::SegPtr b) {
                    return a->start == b->start && a->length == b->length;
                });
            segments.erase(last, segments.end());

            // 验证没有重叠
            bool has_overlap = false;
            for (size_t i = 1; i < segments.size(); ++i) {
                if (segments[i]->start < segments[i-1]->start + segments[i-1]->length) {
                    spdlog::error("[mergeMultipleGraphs] Detected overlapping segments: [{}, {}) and [{}, {})",
                                segments[i-1]->start, segments[i-1]->start + segments[i-1]->length,
                                segments[i]->start, segments[i]->start + segments[i]->length);
                    has_overlap = true;
                }
            }

            if (has_overlap) {
                spdlog::error("[mergeMultipleGraphs] Chromosome {}:{} has overlapping segments, skipping reconstruction",
                            species_name, chr_name);
                continue;
            }

            // 重建链表
            if (!segments.empty()) {
                try {
                    // 重要：先清理segments的prev/next指针，避免残留的旧指针导致链表不一致
                    for (auto& seg : segments) {
                        seg->primary_path.prev.store(nullptr, std::memory_order_relaxed);
                        seg->primary_path.next.store(nullptr, std::memory_order_relaxed);
                    }

                    // 清空现有链表及采样表
                    genome_end.clearAllSegments();

                    // 链接所有segments
                    RaMesh::Segment::linkChain(segments);

                    // 插入链表
                    genome_end.spliceSegmentChain(segments,
                                                segments.front()->start,
                                                segments.back()->start + segments.back()->length);

                    processed_segments += segments.size();
                }
                catch (const std::exception& e) {
                    spdlog::error("[mergeMultipleGraphs] Failed to rebuild linked list: {}", e.what());
                }
            }
        }
    } // 释放全局写锁

    auto phase3_end = std::chrono::high_resolution_clock::now();
    auto phase3_duration = std::chrono::duration_cast<std::chrono::milliseconds>(phase3_end - phase3_start);
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(phase3_end - start_time);

    // 最终统计
    size_t final_block_count = 0;
    {
        std::shared_lock<std::shared_mutex> verify_lock(graph.rw);
        for (const auto& weak_block : graph.blocks) {
            if (!weak_block.expired()) {
                final_block_count++;
            }
        }
    }

#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Phase 3 completed ({} ms)", phase3_duration.count());
    spdlog::info("[mergeMultipleGraphs] Merge completion statistics:");
    spdlog::info("[mergeMultipleGraphs]   - Original blocks: {}", original_blocks_count);
    spdlog::info("[mergeMultipleGraphs]   - Final blocks: {}", final_block_count);
    spdlog::info("[mergeMultipleGraphs]   - Processed segments: {}", processed_segments);
    spdlog::info("[mergeMultipleGraphs]   - Total time: {} ms", total_duration.count());
#endif

    spdlog::info("Multi-genome graph merging completed successfully.");
    
#ifdef _DEBUG_
    spdlog::info("[mergeMultipleGraphs] Multi-genome graph merging completed!");
#endif
}

// 辅助函数：处理基本区间
RaMesh::BlockPtr MultipleRareAligner::processElementaryInterval(
    uint_t interval_start, uint_t interval_end,
    const std::vector<RaMesh::BlockPtr>& valid_blocks,
    const SpeciesName& ref_name)
{
    try {
        // 确定覆盖当前基本区间的源Block集合
        std::vector<RaMesh::BlockPtr> covering_blocks;
        covering_blocks.reserve(valid_blocks.size());

        for (auto& block_ptr : valid_blocks) {
            std::shared_lock block_lock(block_ptr->rw);

            // 检查该block的ref-segment是否覆盖当前基本区间
            for (const auto& [species_chr_pair, segment] : block_ptr->anchors) {
                const auto& [species_name, chr_name] = species_chr_pair;
                if (segment && segment->isSegment() && species_name == ref_name) {
                    uint_t seg_start = segment->start;
                    uint_t seg_end = segment->start + segment->length;

                    // 检查segment是否覆盖基本区间
                    if (seg_start <= interval_start && seg_end >= interval_end) {
                        covering_blocks.push_back(block_ptr);
                        break;
                    }
                }
            }
        }

        if (covering_blocks.empty()) {
            return nullptr;
        }

        // 创建新的Block
        RaMesh::BlockPtr new_block = RaMesh::Block::create(covering_blocks.size() * 2);

        // 处理每个覆盖的源Block
        for (auto& source_block : covering_blocks) {
            std::shared_lock source_lock(source_block->rw);

            // 处理该block中的所有segments
            for (const auto& [species_chr_pair, source_segment] : source_block->anchors) {
                if (!source_segment || !source_segment->isSegment()) continue;

                const auto& [species_name, chr_name] = species_chr_pair;

                if (species_name == ref_name) {
                    // ref-segment：直接映射基本区间
                    uint_t seg_start = source_segment->start;
                    uint_t seg_end = source_segment->start + source_segment->length;

                    if (seg_start <= interval_start && seg_end >= interval_end) {
                        uint_t new_length = interval_end - interval_start;
                        if (new_length == 0) {
                            continue; // 跳过零长度片段，避免后续重叠
                        }
                        // 创建新的ref-segment片段
                        auto new_segment = RaMesh::Segment::create(
                            interval_start,
                            new_length,
                            source_segment->strand,
                            Cigar_t{},  // ref segment通常没有CIGAR
                            source_segment->align_role,
                            source_segment->seg_role,
                            new_block
                        );

                        // 注册到新block
                        {
                            std::unique_lock new_lock(new_block->rw);
                            // 检查是否已存在相同的anchor，避免重复
                            if (new_block->anchors.find(species_chr_pair) == new_block->anchors.end()) {
                                new_block->anchors[species_chr_pair] = new_segment;
                            }
                        }
                    }
                } else {
                    // query-segment：需要映射到ref坐标系统
                    auto [ref_start, ref_end] = getRefMappedInterval(source_segment, source_block, ref_name);

                    if (ref_start == 0 && ref_end == 0) {
                        continue; // 无法获取映射
                    }

                    // 检查query segment的ref映射是否覆盖基本区间
                    if (ref_start <= interval_start && ref_end >= interval_end) {
                        // 计算query区间
                        auto [query_start, query_length] = calculateQueryInterval(
                            source_segment->cigar,
                            interval_start, interval_end,
                            ref_start,
                            source_segment->start
                        );

                        // 跳过零长度query segment
                        if (query_length == 0) {
                            continue;
                        }

                        // 提取对应的CIGAR子序列
                        uint_t ref_offset = interval_start - ref_start;
                        Cigar_t new_cigar = extractCigarSubsequence(
                            source_segment->cigar,
                            ref_offset,
                            interval_end - interval_start
                        );

                        auto new_segment = RaMesh::Segment::create(
                            query_start,
                            query_length,
                            source_segment->strand,
                            new_cigar,
                            source_segment->align_role,
                            source_segment->seg_role,
                            new_block
                        );

                        {
                            std::unique_lock new_lock(new_block->rw);
                            // 检查是否已存在相同的anchor，避免重复
                            if (new_block->anchors.find(species_chr_pair) == new_block->anchors.end()) {
                                new_block->anchors[species_chr_pair] = new_segment;
                            }
                        }
                    }
                }
            }
        }

        return new_block;
    }
    catch (const std::exception& e) {
        spdlog::error("[processElementaryInterval] Error processing interval [{}, {}): {}",
                    interval_start, interval_end, e.what());
        return nullptr;
    }
}
