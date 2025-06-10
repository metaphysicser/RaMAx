#include "anchor.h"
#include "data_process.h" 

// 构造函数（注释掉了）：初始化参考物种、查询物种及其锚点集合，并重建 R 树。
//PairGenomeAnchor::PairGenomeAnchor(SpeciesName ref,
//    SpeciesName query,
//    AnchorPtrListVec init_anchors)
//    : ref_species(ref)
//    , query_species(query)
//    , anchors(std::move(init_anchors))
//{
//    rebuild(); // 构建 R 树索引
//}

// 重建 R 树索引：将 anchors 中所有 Anchor 加入 R 树中。
//
//void PairGenomeAnchor::rebuild() {
//    // 1) 清空已有的 R 树索引
//    anchor_rtree.clear();
//
//    // 2) 遍历 anchors 的每个 list，插入 AnchorPtr 到 R-tree 中
//    for (auto const& list_group : anchors) {
//        for (auto const& ap : list_group) {
//            anchor_rtree.insert({ make_box(*ap), ap });
//        }
//    }
//}

// 将一个 Anchor 插入 R 树中
//void PairGenomeAnchor::insertAnchorRtree(AnchorPtr ap) {
//    anchor_rtree.insert({ make_box(*ap), ap });
//}

// 从 R 树中移除一个 Anchor
//void PairGenomeAnchor::removeAnchorRtree(AnchorPtr ap) {
//    // 从 R 树中删除 Anchor
//    anchor_rtree.remove({ make_box(*ap), ap });
//}

// 查询与给定区域相交的所有锚点
//std::vector<AnchorPtr> PairGenomeAnchor::query(const AnchorBox& region) const {
//    // 在 R 树中查找与 region 相交的 AnchorEntry
//    std::vector<AnchorEntry> hits;
//    anchor_rtree.query(bgi::intersects(region),
//        std::back_inserter(hits));
//
//    // 提取命中的 AnchorPtr 到输出向量中
//    std::vector<AnchorPtr> out;
//    out.reserve(hits.size());
//    for (auto const& e : hits) {
//        out.push_back(e.second);
//    }
//    return out;
//}




//void groupMatchByQueryRef(MatchVec3DPtr& anchors,
//    MatchByQueryRef& unique_anchors,
//    MatchByQueryRef& repeat_anchors,
//    FastaManager& ref_fasta_manager,
//    FastaManager& query_fasta_manager,
//    uint_t thread_num) {
//    //----------------------------------------------------------------
//    // 0) 初始化二维目标矩阵
//    //----------------------------------------------------------------
//    const uint_t ref_chr_cnt = ref_fasta_manager.idx_map.size();
//    const uint_t query_chr_cnt = query_fasta_manager.idx_map.size();
//
//    // resize 外层：query 维度
//    unique_anchors.resize(query_chr_cnt);
//    repeat_anchors.resize(query_chr_cnt);
//
//    // resize 内层：ref 维度
//    for (uint_t q = 0; q < query_chr_cnt; ++q) {
//        unique_anchors[q].resize(ref_chr_cnt);   // AnchorsByRef
//        repeat_anchors[q].resize(ref_chr_cnt);   // AnchorsByRef
//    }
//
//    //----------------------------------------------------------------
//    // 1) 行级互斥锁：每个 query-chr 一把
//    //----------------------------------------------------------------
//    std::vector<std::mutex> rowLocks(query_chr_cnt);
//
//    //----------------------------------------------------------------
//    // 2) 并行遍历 3D 数据，每个 slice 启动一个任务
//    //----------------------------------------------------------------
//    ThreadPool pool(thread_num);
//
//    for (auto& slice : *anchors)              // slice: std::vector<AnchorVec>
//    {
//        pool.enqueue([&, &slice]() {
//
//            if (slice.empty()) return;
//
//            // 2a) 获取 slice 的 query-chr 索引（所有 vec 同一 chr）
//            const Match& first = slice.front().front();
//            auto itQ = query_fasta_manager.idx_map.find(
//                first.query_region.chr_name);
//            if (itQ == query_fasta_manager.idx_map.end()) return;
//            const uint_t qIdx = itQ->second;
//
//            // 2b) 遍历 slice 中每个 AnchorVec
//            for (auto& vec : slice) {
//                if (vec.empty()) continue;
//
//                const ChrName& rName = vec.front().ref_region.chr_name;
//                auto itR = ref_fasta_manager.idx_map.find(rName);
//                if (itR == ref_fasta_manager.idx_map.end()) continue;
//                const uint_t rIdx = itR->second;
//
//                /* ---- 临界区：只锁当前 query 行 ---- */
//                {
//                    std::lock_guard<std::mutex> lk(rowLocks[qIdx]);
//                    MatchVec& target =
//                        (vec.size() == 1)
//                        ? unique_anchors[qIdx][rIdx]
//                        : repeat_anchors[qIdx][rIdx];
//
//                        target.insert(target.end(),
//                            std::make_move_iterator(vec.begin()),
//                            std::make_move_iterator(vec.end()));
//                }
//            }
//            });
//    }
//
//    pool.waitAllTasksDone();
//
//    //----------------------------------------------------------------
//    // 3) 释放原始 3D 数据以节省内存
//    //----------------------------------------------------------------
//    anchors->clear();
//    anchors->shrink_to_fit();
//}

//--------------------------------------------------------------------
// 适配 4‑D 结构的聚类分发
//--------------------------------------------------------------------
void groupMatchByQueryRef(MatchVec3DPtr& anchors,
    MatchByStrandByQueryRefPtr unique_anchors,
    MatchByStrandByQueryRefPtr repeat_anchors,
    FastaManager& ref_fasta_manager,
    FastaManager& query_fasta_manager,
	ThreadPool& pool)
{
    //------------------------------------------------------------
    // 0) 初始化 4‑D 输出矩阵  [strand][query][ref]
    //------------------------------------------------------------
    constexpr uint_t STRAND_CNT = 2;            // FORWARD=0, REVERSE=1

    const uint_t ref_chr_cnt = static_cast<uint_t>(ref_fasta_manager.idx_map.size());
    const uint_t query_chr_cnt = static_cast<uint_t>(query_fasta_manager.idx_map.size());

    auto initTarget = [&](MatchByStrandByQueryRefPtr& tgt) {
        (*tgt).resize(STRAND_CNT);
        for (uint_t s = 0; s < STRAND_CNT; ++s) {
            (*tgt)[s].resize(query_chr_cnt);              // [query]
            for (uint_t q = 0; q < query_chr_cnt; ++q)
                (*tgt)[s][q].resize(ref_chr_cnt);         // [ref]
        }
        };

    initTarget(unique_anchors);
    initTarget(repeat_anchors);

    //------------------------------------------------------------
    // 1) 行级互斥锁: 针对 (strand, query) 组合各一把锁
    //------------------------------------------------------------
    std::vector<std::mutex> rowLocks(STRAND_CNT * query_chr_cnt);
    auto rowLockIdx = [&](uint_t strand, uint_t qIdx) {
        return strand * query_chr_cnt + qIdx;
        };

    //------------------------------------------------------------
    // 2) 并行遍历 3‑D 数据，slice 级别开任务
    //------------------------------------------------------------

    for (auto& slice : *anchors)
    {
        pool.enqueue([&, slice]() {              // ← 值捕获
            if (slice.empty()) return;

            const Match& first = slice.front().front();

            // 1. strand 索引归一化
            uint_t sIdx = (first.strand == REVERSE ? 1u : 0u);

            auto itQ = query_fasta_manager.idx_map.find(first.query_region.chr_name);
            if (itQ == query_fasta_manager.idx_map.end()) return;
            uint_t qIdx = itQ->second;

            for (auto& vec : slice)
            {
                if (vec.empty()) continue;
                auto itR = ref_fasta_manager.idx_map.find(vec.front().ref_region.chr_name);
                if (itR == ref_fasta_manager.idx_map.end()) continue;
                uint_t rIdx = itR->second;

                std::lock_guard<std::mutex> lk(rowLocks[rowLockIdx(sIdx, qIdx)]);
                MatchVec& tgt = (vec.size() == 1)
                    ? (*unique_anchors)[sIdx][qIdx][rIdx]
                    : (*repeat_anchors)[sIdx][qIdx][rIdx];
                    tgt.insert(tgt.end(),
                        std::make_move_iterator(vec.begin()),
                        std::make_move_iterator(vec.end()));
            }
            });
    }


    pool.waitAllTasksDone();

    //------------------------------------------------------------
    // 3) 释放原始 3‑D 数据以节省内存
    //------------------------------------------------------------
    anchors->clear();
    anchors->shrink_to_fit();
}

//void sortMatchByRefStart(MatchByQueryRefPtr& anchors, ThreadPool& pool) {
//
//    for (auto& row : *anchors) {
//        for (auto& vec : row) {
//            pool.enqueue([&vec]() {
//                if (vec.size() == 0) return;
//                std::sort(vec.begin(), vec.end(), [](const Match& a, const Match& b) {
//                    return (a.ref_region.start < b.ref_region.start) || (a.ref_region.start == b.ref_region.start && a.query_region.start < b.query_region.start);
//                    });
//                });
//        }
//    }
//
//    pool.waitAllTasksDone();
//}

// anchors[strand][query][ref] -- MatchByStrandByQueryRef
void sortMatchByQueryStart(MatchByStrandByQueryRefPtr& anchors, ThreadPool& pool)
{

    for (auto& strand_layer : *anchors)          // 第一维：Strand
        for (auto& query_row : strand_layer)    // 第二维：query-chr
            for (auto& vec : query_row)         // 第三维：ref-chr
            {
                // 用指针值捕获，避免 vec 变量在下一轮循环被复用后悬空
                pool.enqueue([v = &vec]() {
                    if (v->empty()) return;
                    std::sort(v->begin(), v->end(),
                        [](const Match& a, const Match& b)
                        { return a.query_region.start < b.query_region.start; });
                    });
            }

    pool.waitAllTasksDone();
}


AnchorPtrVec findNonOverlapAnchors(const AnchorVec& anchors)
{
    AnchorPtrVec uniques;
    const size_t n = anchors.size();
    if (n == 0) return uniques;

    Coord_t max_end_so_far = 0;  // 记录到当前位置前的最大 end

    for (size_t i = 0; i < n; ++i) {
        const Anchor& curr = anchors[i];
        Coord_t curr_start = curr.match.ref_region.start;
        Coord_t curr_end = curr_start + curr.match.ref_region.length;

        bool overlap_left = (i > 0 && max_end_so_far > curr_start);
        bool overlap_right = (i + 1 < n &&
            anchors[i + 1].match.ref_region.start < curr_end);

        if (!overlap_left && !overlap_right)
            uniques.emplace_back(std::make_shared<Anchor>(curr));

        // 更新左侧最大 end
        if (curr_end > max_end_so_far) max_end_so_far = curr_end;
    }
    return std::move(uniques);
}






