#include "anchor.h"
#include "data_process.h" 

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
    // 0) 初始化输出矩阵 [strand][query][ref]
    //------------------------------------------------------------
    constexpr uint_t STRAND_CNT = 2;                         // 0=FWD,1=REV
    const uint_t ref_chr_cnt = static_cast<uint_t>(ref_fasta_manager.idx_map.size());
    const uint_t query_chr_cnt = static_cast<uint_t>(query_fasta_manager.idx_map.size());

    auto initTarget = [&](MatchByStrandByQueryRefPtr& tgt) {
        tgt->resize(STRAND_CNT);
        for (uint_t s = 0; s < STRAND_CNT; ++s) {
            (*tgt)[s].resize(query_chr_cnt);
            for (uint_t q = 0; q < query_chr_cnt; ++q)
                (*tgt)[s][q].resize(ref_chr_cnt);
        }
        };
    initTarget(unique_anchors);
    initTarget(repeat_anchors);

    //------------------------------------------------------------
    // 1) 行级互斥锁：放到堆上，用 shared_ptr 延长生命周期
    //------------------------------------------------------------
    auto rowLocks = std::make_shared<std::vector<std::mutex>>(STRAND_CNT * query_chr_cnt);

    //------------------------------------------------------------
    // 2) 并行遍历 3-D 数据，slice 级别开任务
    //------------------------------------------------------------
    for (auto& slice : *anchors) {
        // 按值捕获 slice、rowLocks、unique_anchors 等
        pool.enqueue(
            [slice,                                               // 值
            rowLocks,                                            // 值 (shared_ptr)
            unique_anchors, repeat_anchors,                      // 值 (shared_ptr)
            &ref_fasta_manager, &query_fasta_manager,            // 引用
            query_chr_cnt] () mutable
            {
                if (slice.empty()) return;

                const Match& first = slice.front().front();
                uint_t sIdx = (first.strand == REVERSE ? 1u : 0u);

                auto itQ = query_fasta_manager.idx_map.find(first.query_region.chr_name);
                if (itQ == query_fasta_manager.idx_map.end()) return;
                uint_t qIdx = itQ->second;

                for (auto& vec : slice) {
                    if (vec.empty()) continue;

                    auto itR = ref_fasta_manager.idx_map.find(vec.front().ref_region.chr_name);
                    if (itR == ref_fasta_manager.idx_map.end()) continue;
                    uint_t rIdx = itR->second;

                    // 计算锁下标并加锁
                    uint_t lockIdx = sIdx * query_chr_cnt + qIdx;
                    std::lock_guard<std::mutex> lk((*rowLocks)[lockIdx]);

                    MatchVec& tgt = (vec.size() == 1)
                        ? (*unique_anchors)[sIdx][qIdx][rIdx]
                        : (*repeat_anchors)[sIdx][qIdx][rIdx];

                        tgt.insert(tgt.end(),
                            std::make_move_iterator(vec.begin()),
                            std::make_move_iterator(vec.end()));
                }
            });
    }
    // **不再在这里 wait；由调用者在外部 pool.waitAllTasksDone() 同步**
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







