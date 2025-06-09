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

void MultipleRareAligner::starAlignment(uint_t tree_root)
{
    std::vector leaf_vec = newick_tree.orderLeavesGreedyMinSum(tree_root);
	uint_t leaf_num = leaf_vec.size();

    for (uint_t i = 0; i < leaf_num; i++) {
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

    }
    return;
}

SpeciesMatchVec3DPtrMapPtr MultipleRareAligner::alignMultipleQuerys(
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
        auto mp = std::make_shared<SpeciesMatchVec3DPtrMap>();
        if (loadSpeciesMatchMap(anchor_file, mp))
            return mp;
        // 如果读取失败则继续重新计算
    }

    /* ---------- 3. 准备参考基因组索引 ---------- */
    FilePath ref_index_path = index_dir / ref_name;
    std::filesystem::create_directories(ref_index_path);

    FilePath ref_fasta = species_fasta_manager_map[ref_name].fasta_path_;
    FilePath out_fa = ref_index_path / (ref_name + ref_fasta.extension().string());
    FilePath idx_file = ref_index_path / (ref_name + "." + FMINDEX_EXTESION);

    spdlog::info("Indexing {} ...", ref_name);
    FM_Index ref_index(ref_name, &species_fasta_manager_map[ref_name]);

    if (!std::filesystem::exists(idx_file)) {
        ref_index.buildIndex(out_fa, fast_build, thread_num);
        ref_index.saveToFile(idx_file.string());
    }
    else {
        ref_index.loadFromFile(idx_file.string());
    }
    spdlog::info("Index ready: {}", idx_file.string());

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
                [this, &ref_index, prefix, &fm, search_mode, allow_MEM, &shared_pool]() -> MatchVec3DPtr {
                    // 统一走单物种比对逻辑，公用 shared_pool
                    return this->alignSingleQuerys(prefix, ref_index, fm, search_mode, allow_MEM, shared_pool);
                })
        );
    }

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
    saveSpeciesMatchMap(anchor_file, result_map);
    spdlog::info("[alignMultipleQuerys] all species done. Saved to {}", anchor_file.string());

    return result_map;
}



MatchVec3DPtr MultipleRareAligner::alignSingleQuerys(
    const std::string& prefix,
    FM_Index& ref_index,
    FastaManager& query_fm,
    SearchMode         search_mode,
    bool               allow_MEM,
    ThreadPool& pool)        // ← 外部共享线程池
{
    /* ---------- 1. 结果文件路径，与多基因组保持同一目录 ---------- */
    FilePath result_dir = work_dir / RESULT_DIR
        / ("group_" + std::to_string(group_id))
        / ("round_" + std::to_string(round_id));
    std::filesystem::create_directories(result_dir);

    FilePath anchor_file = result_dir /
        (prefix + "_" + SearchModeToString(search_mode) + "." + ANCHOR_EXTENSION);

    /* ---------- 2. 若文件已存在，直接加载 ---------- */
    if (std::filesystem::exists(anchor_file)) {
        auto mv3 = std::make_shared<MatchVec3D>();
        if (loadMatchVec3D(anchor_file, mv3))
            return mv3;
        // 读取失败则重新计算
    }

    /* ---------- 3. FASTA 分块 ---------- */
    RegionVec chunks = query_fm.preAllocateChunks(chunk_size, overlap_size);
    if (chunks.empty()) {
        spdlog::warn("[alignSingleQuerys] no chunks, prefix {}", prefix);
        return std::make_shared<MatchVec3D>();
    }

    /* ---------- 4. 把块任务丢进共享线程池 ---------- */
    size_t num_chunks_per_task =
        std::max<size_t>(1, (chunks.size() + thread_num - 1) / thread_num);

    std::vector<std::future<MatchVec2DPtr>> futures;
    futures.reserve((chunks.size() + num_chunks_per_task - 1) / num_chunks_per_task * 2);

    for (size_t i = 0; i < chunks.size(); i += num_chunks_per_task) {
        RegionVec grp(chunks.begin() + i,
            chunks.begin() + std::min(i + num_chunks_per_task, chunks.size()));

        // forward
        futures.emplace_back(pool.enqueue(
            [this, grp, &ref_index, &query_fm, search_mode, allow_MEM]() -> MatchVec2DPtr {
                auto g = std::make_shared<MatchVec2D>();
                for (const auto& ck : grp) {
                    std::string seq = query_fm.getSubSequence(ck.chr_name, ck.start, ck.length);
                    auto mv = ref_index.findAnchors(
                        ck.chr_name, seq,
                        search_mode,
                        Strand::FORWARD,
                        allow_MEM,
                        ck.start,
                        0,
                        max_anchor_frequency);
                    g->insert(g->end(), mv->begin(), mv->end());
                }
                return g;
            }));

        // reverse
        futures.emplace_back(pool.enqueue(
            [this, grp, &ref_index, &query_fm, allow_MEM]() -> MatchVec2DPtr {
                auto g = std::make_shared<MatchVec2D>();
                for (const auto& ck : grp) {
                    std::string seq = query_fm.getSubSequence(ck.chr_name, ck.start, ck.length);
                    auto mv = ref_index.findAnchors(
                        ck.chr_name, seq,
                        MIDDLE_SEARCH,
                        Strand::REVERSE,
                        allow_MEM,
                        ck.start,
                        min_anchor_length,
                        max_anchor_frequency);
                    g->insert(g->end(), mv->begin(), mv->end());
                }
                return g;
            }));
    }

    /* ---------- 5. 收集结果 ---------- */
    auto mv3 = std::make_shared<MatchVec3D>();
    mv3->reserve(futures.size());

    for (auto& fut : futures) {
        MatchVec2DPtr part = fut.get();          // 等待任务完成
        mv3->emplace_back(std::move(*part));     // Move 到结果
    }

    /* ---------- 6. 保存到文件 ---------- */
    saveMatchVec3D(anchor_file, mv3);
    return mv3;
}

