#include "anchor.h"

#include <SeqPro.h>

#include "data_process.h"

////--------------------------------------------------------------------
//// 适配 4‑D 结构的聚类分发
////--------------------------------------------------------------------
//void groupMatchByQueryRef(MatchVec3DPtr& anchors,
//    MatchByStrandByQueryRefPtr unique_anchors,
//    MatchByStrandByQueryRefPtr repeat_anchors,
//    SeqPro::ManagerVariant& ref_fasta_manager,
//    SeqPro::ManagerVariant& query_fasta_manager,
//	ThreadPool& pool)
//{
//    //------------------------------------------------------------
//    // 0) 初始化输出矩阵 [strand][query][ref]
//    //------------------------------------------------------------
//    constexpr uint_t STRAND_CNT = 2;                         // 0=FWD,1=REV
//    const uint_t ref_chr_cnt = std::visit([](auto&& manager) {
//        using PtrType = std::decay_t<decltype(manager)>;
//        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
//            return manager->getSequenceCount();
//        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//            return manager->getSequenceCount();
//        } else {
//            throw std::runtime_error("Unhandled manager type in variant.");
//        }
//    }, ref_fasta_manager);
//
//    const uint_t query_chr_cnt = std::visit([](auto&& manager) {
//        using PtrType = std::decay_t<decltype(manager)>;
//        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
//            return manager->getSequenceCount();
//        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//            return manager->getSequenceCount();
//        } else {
//            throw std::runtime_error("Unhandled manager type in variant.");
//        }
//    }, query_fasta_manager);
//    auto initTarget = [&](MatchByStrandByQueryRefPtr& tgt) {
//        tgt->resize(STRAND_CNT);
//        for (uint_t s = 0; s < STRAND_CNT; ++s) {
//            (*tgt)[s].resize(query_chr_cnt);
//            for (uint_t q = 0; q < query_chr_cnt; ++q)
//                (*tgt)[s][q].resize(ref_chr_cnt);
//        }
//        };
//    initTarget(unique_anchors);
//    initTarget(repeat_anchors);
//
//    //------------------------------------------------------------
//    // 1) 行级互斥锁：放到堆上，用 shared_ptr 延长生命周期
//    //------------------------------------------------------------
//    auto rowLocks = std::make_shared<std::vector<std::mutex>>(STRAND_CNT * query_chr_cnt);
//
//    //------------------------------------------------------------
//    // 2) 并行遍历 3-D 数据，slice 级别开任务
//    //------------------------------------------------------------
//    for (auto& slice : *anchors) {
//        // 按值捕获 slice、rowLocks、unique_anchors 等
//        pool.enqueue(
//            [slice,                                               // 值
//            rowLocks,                                            // 值 (shared_ptr)
//            unique_anchors, repeat_anchors,                      // 值 (shared_ptr)
//            &ref_fasta_manager, &query_fasta_manager,            // 引用
//            query_chr_cnt] () mutable
//            {
//                if (slice.empty()) return;
//
//                const Match& first = slice.front().front();
//                uint_t sIdx = (first.strand == REVERSE ? 1u : 0u);
//
//            uint_t qIdx = std::visit([&first](auto&& manager) {
//                using PtrType = std::decay_t<decltype(manager)>;
//                if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
//                    return manager->getSequenceId(first.query_region.chr_name);
//                } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//                    return manager->getSequenceId(first.query_region.chr_name);
//                } else {
//                    throw std::runtime_error("Unhandled manager type in variant.");
//                }
//            }, query_fasta_manager);
//
//            if(qIdx == SeqPro::SequenceIndex::INVALID_ID) return;
//
//            for (auto& vec : slice)
//            {
//                if (vec.empty()) continue;
//
//                uint_t rIdx = std::visit([&vec](auto&& manager) {
//                    using PtrType = std::decay_t<decltype(manager)>;
//                    if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
//                        return manager->getSequenceId(vec.front().ref_region.chr_name);
//                    } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
//                        return manager->getSequenceId(vec.front().ref_region.chr_name);
//                    } else {
//                        throw std::runtime_error("Unhandled manager type in variant.");
//                    }
//                }, ref_fasta_manager);
//
//                if(rIdx == SeqPro::SequenceIndex::INVALID_ID) continue;
//
//                // 计算锁下标并加锁
//                uint_t lockIdx = sIdx * query_chr_cnt + qIdx;
//                std::lock_guard<std::mutex> lk((*rowLocks)[lockIdx]);
//
//                MatchVec& tgt = (vec.size() == 1)
//                    ? (*unique_anchors)[sIdx][qIdx][rIdx]
//                    : (*repeat_anchors)[sIdx][qIdx][rIdx];
//
//                    tgt.insert(tgt.end(),
//                        std::make_move_iterator(vec.begin()),
//                        std::make_move_iterator(vec.end()));
//                    vec.clear();
//                    vec.shrink_to_fit();
//                }
//            });
//    }
//    
//}
// ---------------------------------------------------------------------------
// groupMatchByQueryRef – SERIAL, LOW‑MEMORY VERSION
// ---------------------------------------------------------------------------
// - 完全取消并行，不依赖 ThreadPool，也不需要锁。
// - 按需扩展 3‑D 目标结构，避免一次性为空槽分配 vector。
// - 每次把 slice 内部 vec 的元素 move 到目标后，立即 clear()+shrink_to_fit()
//   以释放多余 capacity，降低峰值 RSS。
// ---------------------------------------------------------------------------
// 依赖类型说明（保持与原工程一致）：
//   using Match               = ...;
//   using MatchVec            = std::vector<Match>;
//   using MatchVec2D          = std::vector<MatchVec>;                // [ref]
//   using MatchVec3D          = std::vector<MatchVec2D>;              // [query][ref]
//   using MatchVec3DPtr       = std::shared_ptr<MatchVec3D>;          // slice 别名
//   using MatchByStrandByQueryRef = std::vector<MatchVec3D>;          // [strand][query][ref]
//   using MatchByStrandByQueryRefPtr = std::shared_ptr<MatchByStrandByQueryRef>;
//   enum Strand { FORWARD, REVERSE };
// ---------------------------------------------------------------------------

void groupMatchByQueryRef(
    MatchVec3DPtr& anchors,
    MatchByStrandByQueryRefPtr unique_anchors,
    MatchByStrandByQueryRefPtr repeat_anchors,
    SeqPro::ManagerVariant& ref_fasta_manager,
    SeqPro::ManagerVariant& query_fasta_manager,
    ThreadPool & pool)
{
    constexpr uint_t STRAND_CNT = 2; // 0 = FWD, 1 = REV

    // ----------------- 0) 获取染色体数 -----------------
    const uint_t ref_chr_cnt = std::visit([](auto& m) { return m->getSequenceCount(); }, ref_fasta_manager);
    const uint_t qry_chr_cnt = std::visit([](auto& m) { return m->getSequenceCount(); }, query_fasta_manager);

    // ----------------- 1) 按需扩展工具 lambda -----------------
    auto ensure_slot = [&](MatchByStrandByQueryRefPtr& tgt,
        uint_t s, uint_t q, uint_t r) -> MatchVec&
        {
            //if (s >= tgt->size())            tgt->resize(STRAND_CNT);
            //if (q >= (*tgt)[s].size())       (*tgt)[s].resize(qry_chr_cnt);
            //if (r >= (*tgt)[s][q].size())    (*tgt)[s][q].resize(ref_chr_cnt);
            return (*tgt)[s][q][r];
        };

    auto initTarget = [&](MatchByStrandByQueryRefPtr& tgt) {
    tgt->resize(STRAND_CNT);
    for (uint_t s = 0; s < STRAND_CNT; ++s) {
        (*tgt)[s].resize(qry_chr_cnt);
        for (uint_t q = 0; q < qry_chr_cnt; ++q)
            (*tgt)[s][q].resize(ref_chr_cnt);
    }
    };
    initTarget(unique_anchors);
    initTarget(repeat_anchors);

    // ----------------- 2) 顺序遍历所有 slice -----------------
    for (auto& slice : *anchors) {
        if (slice.empty()) continue;

        // sIdx & qIdx 只需要从 slice 第一个元素即可得出
        const Match& first = slice.front().front();
        const uint_t sIdx = (first.strand() == REVERSE ? 1u : 0u);

        ChrIndex qIdx = first.qry_chr_index;
        if (qIdx == SeqPro::SequenceIndex::INVALID_ID) continue;

        // --- 遍历 slice 内的每个 MatchVec（同一 query & strand, 不同 ref）
        for (auto& vec : slice) {
            if (vec.empty()) continue;

            ChrIndex rIdx = vec.front().ref_chr_index;
            if (rIdx == SeqPro::SequenceIndex::INVALID_ID) continue;

            MatchVec& dest = (vec.size() == 1)
                ? ensure_slot(unique_anchors, sIdx, qIdx, rIdx)
                : ensure_slot(repeat_anchors, sIdx, qIdx, rIdx);

            if (dest.empty()) dest.reserve(vec.size()); // 减少后续扩容

            dest.insert(dest.end(),
                std::make_move_iterator(vec.begin()),
                std::make_move_iterator(vec.end()));

            // --- 立即回收 vec 占用容量 ---
            vec.clear();
            vec.shrink_to_fit();
        }
    }

    // ----------------- 3) anchors 自身可以释放 -----------------
    anchors->clear();
    anchors->shrink_to_fit();

    // 遍历unique_anchors，使用shrink_to_fit()
    auto shrink_matrix = [](MatchByStrandByQueryRefPtr& matrix) {
        for (auto& strand_v : *matrix) {
            for (auto& qry_v : strand_v) {
                for (auto& ref_v : qry_v) {
                    if (ref_v.capacity() > ref_v.size())
                        ref_v.shrink_to_fit();
                }
            }
        }
        };
    shrink_matrix(unique_anchors);
    shrink_matrix(repeat_anchors);
}



// 重载版本：支持 SharedManagerVariant，通过解引用调用原始版本
void groupMatchByQueryRef(MatchVec3DPtr& anchors,
    MatchByStrandByQueryRefPtr unique_anchors,
    MatchByStrandByQueryRefPtr repeat_anchors,
    SeqPro::SharedManagerVariant& ref_fasta_manager,
    SeqPro::SharedManagerVariant& query_fasta_manager,
    ThreadPool& pool)
{
    // 解引用 SharedManagerVariant 并调用原始函数
    groupMatchByQueryRef(anchors, unique_anchors, repeat_anchors,
                        *ref_fasta_manager, *query_fasta_manager, pool);
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
    constexpr size_t MIN_SIZE_FOR_PARALLEL = 100;  // 阈值：小于100个元素直接串行排序
    constexpr size_t BATCH_SIZE = 50;              // 批处理大小
    
    std::vector<std::future<void>> futures;
    std::vector<std::vector<MatchVec*>> batches;
    std::vector<MatchVec*> current_batch;
    
    // 收集所有需要排序的向量，并按大小分类
    for (auto& strand_layer : *anchors) {
        for (auto& query_row : strand_layer) {
            for (auto& vec : query_row) {
                if (vec.empty()) continue;
                
                if (vec.size() >= MIN_SIZE_FOR_PARALLEL) {
                    // 大向量：单独提交任务
                    futures.emplace_back(pool.enqueue([v = &vec]() {
                        std::sort(v->begin(), v->end(),
                            [](const Match& a, const Match& b) {
                                return a.qry_start < b.qry_start;
                            });
                    }));
                } else {
                    // 小向量：加入批处理
                    current_batch.push_back(&vec);
                    if (current_batch.size() >= BATCH_SIZE) {
                        batches.push_back(std::move(current_batch));
                        current_batch.clear();
                    }
                }
            }
        }
    }
    
    // 处理剩余的小向量
    if (!current_batch.empty()) {
        batches.push_back(std::move(current_batch));
    }
    
    // 批处理小向量
    for (auto& batch : batches) {
        futures.emplace_back(pool.enqueue([batch = std::move(batch)]() {
            for (auto* vec : batch) {
                std::sort(vec->begin(), vec->end(),
                    [](const Match& a, const Match& b) {
                        return a.qry_start < b.qry_start;
                    });
            }
        }));
    }
    
    // 等待所有任务完成
    for (auto& future : futures) {
        future.get();
    }
}

//
//AnchorPtrVec findNonOverlapAnchors(const AnchorVec& anchors)
//{
//    AnchorPtrVec uniques;
//    const size_t n = anchors.size();
//    if (n == 0) return uniques;
//
//    Coord_t max_end_so_far = 0;  // 记录到当前位置前的最大 end
//
//    for (size_t i = 0; i < n; ++i) {
//        const Anchor& curr = anchors[i];
//        Coord_t curr_start = curr.match.ref_region.start;
//        Coord_t curr_end = curr_start + curr.match.ref_region.length;
//
//        bool overlap_left = (i > 0 && max_end_so_far > curr_start);
//        bool overlap_right = (i + 1 < n &&
//            anchors[i + 1].match.ref_region.start < curr_end);
//
//        if (!overlap_left && !overlap_right)
//            uniques.emplace_back(std::make_shared<Anchor>(curr));
//
//        // 更新左侧最大 end
//        if (curr_end > max_end_so_far) max_end_so_far = curr_end;
//    }
//    return std::move(uniques);
//}

MatchClusterVecPtr
groupClustersToVec(const ClusterVecPtrByStrandByQueryRefPtr& src,
    ThreadPool& pool, uint_t thread_num)
{
    auto result = std::make_shared<MatchClusterVec>();
    
    if (!src || src->empty()) {
        return result;
    }
    
    // 预计算总的簇数量，用于预分配空间
    size_t total_clusters = 0;
    for (const auto& strand_layer : *src) {
        for (const auto& query_row : strand_layer) {
            for (const auto& cluster_vec_ptr : query_row) {
                if (cluster_vec_ptr && !cluster_vec_ptr->empty()) {
                    total_clusters += cluster_vec_ptr->size();
                }
            }
        }
    }
    
    result->reserve(total_clusters);
    
    // 使用线程安全的方式收集所有簇
    std::mutex result_mutex;
    std::vector<std::future<std::vector<MatchCluster>>> futures;
    
    // 将收集任务分配给多个线程
    size_t tasks_per_thread = std::max(size_t(1), total_clusters / thread_num);
    size_t current_task_count = 0;
    std::vector<MatchCluster> current_batch;
    
    for (const auto& strand_layer : *src) {
        for (const auto& query_row : strand_layer) {
            for (const auto& cluster_vec_ptr : query_row) {
                if (!cluster_vec_ptr || cluster_vec_ptr->empty()) continue;
                
                for (const auto& cluster : *cluster_vec_ptr) {
                    if (!cluster.empty()) {
                        current_batch.push_back(cluster);
                        current_task_count++;
                        
                        // 当批次足够大时，提交到线程池
                        if (current_task_count >= tasks_per_thread) {
                            auto batch_copy = current_batch;  // 复制用于线程安全
                            futures.emplace_back(pool.enqueue([batch_copy]() -> std::vector<MatchCluster> {
                                return batch_copy;
                            }));
                            
                            current_batch.clear();
                            current_task_count = 0;
                        }
                    }
                }
            }
        }
    }
    
    // 处理剩余的批次
    if (!current_batch.empty()) {
        futures.emplace_back(pool.enqueue([current_batch]() -> std::vector<MatchCluster> {
            return current_batch;
        }));
    }
    
    // 收集所有结果
    for (auto& future : futures) {
        auto batch_result = future.get();
        for (auto& cluster : batch_result) {
            result->emplace_back(std::move(cluster));
        }
    }
    
    return result;
}

// ------------------------------------------------------------------
// 把 3D: [strand][queryRef][ref] 的聚簇重新组织成 1D: [ref]
// ------------------------------------------------------------------
ClusterVecPtrByRefPtr groupClustersToRefVec(
    const ClusterVecPtrByStrandByQueryRefPtr& src,
    ThreadPool& pool, uint_t thread_num)
{
    // ---------- 1. 空输入 ----------
    if (!src || src->empty())
        return std::make_shared<ClusterVecPtrByRef>();

    // ---------- 2. 探测 ref 数 ----------
    size_t n_ref = 0;
    for (const auto& strandVec : *src) {
        for (const auto& queryVec : strandVec) {
            if (!queryVec.empty()) {
                n_ref = queryVec.size();
                break;
            }
        }
        if (n_ref) break;
    }
    if (n_ref == 0)
        return std::make_shared<ClusterVecPtrByRef>();

    // ---------- 3. 目标 ----------
    auto dst = std::make_shared<ClusterVecPtrByRef>(n_ref, nullptr);

    /* ******************************************************************
     * 4. 决定是否并行
     *    只能使用 pool 已有的信息，不能改 ThreadPool 接口。
     *    - 大部分线程池都会暴露  thread_count()/size()/parallelism() 之类的读接口；
     *      如果你的实现没有，可以在 ThreadPool 里加一个 const 方法，
     *      但这属于“实现细节”而不是“函数接口”变更，外部签名不动。
     ******************************************************************/
    const bool can_parallel =
        n_ref > 1 &&   // 工作量足够大才值得并行
        thread_num > 1;   // 剩余线程数至少 1 条

    /* =================== 顺序分支 =================== */
    if (!can_parallel) {
        for (size_t ref_id = 0; ref_id < n_ref; ++ref_id) {
            auto combined = std::make_shared<MatchClusterVec>();
            for (const auto& strandVec : *src)
                for (const auto& queryVec : strandVec)
                    if (ref_id < queryVec.size() && queryVec[ref_id])
                        combined->insert(combined->end(),
                            queryVec[ref_id]->begin(),
                            queryVec[ref_id]->end());
            (*dst)[ref_id] = std::move(combined);
        }
        return dst;                   // ✅ 单线程直接结束
    }

    /* =================== 并行分支（原逻辑） =================== */
    std::vector<std::future<void>> futures;
    futures.reserve(n_ref);

    for (size_t ref_id = 0; ref_id < n_ref; ++ref_id) {
        futures.emplace_back(
            pool.enqueue([&, ref_id] {
                auto combined = std::make_shared<MatchClusterVec>();
                for (const auto& strandVec : *src)
                    for (const auto& queryVec : strandVec)
                        if (ref_id < queryVec.size() && queryVec[ref_id])
                            combined->insert(combined->end(),
                                queryVec[ref_id]->begin(),
                                queryVec[ref_id]->end());
                (*dst)[ref_id] = std::move(combined);
                }));
    }
    for (auto& fut : futures) fut.get();
    return dst;
}

// ------------------------------------------------------------------
// 把 3D: [strand][queryRef][ref] 的聚簇重新组织成 1D: [ref]
// ------------------------------------------------------------------
ClusterVecPtrByRefPtr
groupClustersToRefVec(const ClusterVecPtrByStrandByQueryRefPtr& src,
    ThreadPool& pool)
{
    // ---------- 1. 空输入快速返回 ----------
    if (!src || src->empty()) {
        return std::make_shared<ClusterVecPtrByRef>();
    }

    // ---------- 2. 探测参考染色体 (Ref) 数 ----------
    size_t n_ref = 0;
    for (const auto& strandVec : *src) {
        for (const auto& queryVec : strandVec) {
            if (!queryVec.empty()) {
                n_ref = queryVec.size();   // 第一次遇到非空即足够
                break;
            }
        }
        if (n_ref) break;
    }
    if (n_ref == 0) {
        return std::make_shared<ClusterVecPtrByRef>();   // 没有 Ref
    }

    // ---------- 3. 准备目标结构 ----------
    auto dst = std::make_shared<ClusterVecPtrByRef>(n_ref, nullptr);

    // ---------- 4. 为 “每个 Ref ⇨ 一条任务” ----------
    std::vector<std::future<void>> futures;
    futures.reserve(n_ref);

    for (size_t ref_id = 0; ref_id < n_ref; ++ref_id) {
        futures.emplace_back(
            pool.enqueue([&, ref_id] {
                auto combined = std::make_shared<MatchClusterVec>();

                // 聚合：遍历所有 [strand][query]，把 ref_id 处的聚簇并入
                for (const auto& strandVec : *src) {
                    for (const auto& queryVec : strandVec) {
                        if (ref_id < queryVec.size() && queryVec[ref_id]) {
                            combined->insert(combined->end(),
                                queryVec[ref_id]->begin(),
                                queryVec[ref_id]->end());
                        }
                    }
                }
                // 无竞态：独占写入各自槽位
                (*dst)[ref_id] = std::move(combined);
                })
        );
    }

    // ---------- 5. 等待全部任务完成 ----------
    for (auto& fut : futures) fut.get();

    return dst;
}




//--------------------------------------------------------------------
// Anchor 验证功能实现
//--------------------------------------------------------------------

ValidationResult validateAnchorsCorrectness(
    const MatchVec3DPtr& anchors,
    const SeqPro::ManagerVariant& ref_manager,
    const SeqPro::ManagerVariant& query_manager
) {
    spdlog::info("开始验证 anchors 结果的正确性…");
    
    ValidationResult result;
    
    // 反向互补函数
    auto reverseComplement = [](const std::string &seq) -> std::string {
        std::string result;
        result.reserve(seq.length());
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            result.push_back(BASE_COMPLEMENT[static_cast<unsigned char>(*it)]);
        }
        return result;
    };

    // 计算总工作项数
    uint64_t total_work_items = 0;
    for (const auto &level1: *anchors) {
        for (const auto &level2: level1) {
            total_work_items += level2.size();
        }
    }

    if (total_work_items == 0) {
        spdlog::info("没有需要验证的匹配项。");
        return result;
    }

    spdlog::info("总计需要验证 {} 个匹配项。", total_work_items);

    // 进度报告相关变量
    std::atomic<uint64_t> processed_items = 0;
    const uint64_t report_interval = std::max(1ULL, total_work_items / 100ULL);
    std::atomic<bool> diagnostic_dump_done = false;

    // 用于OpenMP reduction的临时变量
    uint64_t total_matches = 0;
    uint64_t correct_matches = 0;
    uint64_t incorrect_matches = 0;

#pragma omp parallel reduction(+: total_matches, correct_matches, incorrect_matches)
    {
#pragma omp for schedule(dynamic) nowait
        for (size_t i = 0; i < anchors->size(); ++i) {
            for (size_t j = 0; j < (*anchors)[i].size(); ++j) {
                for (size_t k = 0; k < (*anchors)[i][j].size(); ++k) {
                    uint64_t current_processed = processed_items.fetch_add(1, std::memory_order_relaxed) + 1;

                    if (current_processed % report_interval == 0) {
                        spdlog::info("验证进度: {} / {} ({:.2f}%)",
                                     current_processed,
                                     total_work_items,
                                     (100.0 * current_processed / total_work_items));
                    }

                    const Match &match = (*anchors)[i][j][k];
                    ++total_matches;

                    try {
                        // 提取参考序列
                        std::string ref_seq = std::visit([&match](auto &&manager_ptr) -> std::string {
                            using PtrType = std::decay_t<decltype(manager_ptr)>;
                            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                                return manager_ptr->getSubSequence(match.ref_chr_index, match.ref_start,
                                                                   match.match_len());
                            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                                SeqPro::MaskedSequenceManager> >) {
                                return manager_ptr->getSubSequence(match.ref_chr_index, match.ref_start,
                                    match.match_len());
                            } else {
                                throw std::runtime_error("Unhandled manager type in variant.");
                            }
                        }, ref_manager);
                        
                        // 提取查询序列
                        std::string query_seq = std::visit([&match](auto &&manager_ptr) -> std::string {
                            using PtrType = std::decay_t<decltype(manager_ptr)>;
                            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                                return manager_ptr->getSubSequence(match.qry_chr_index,
                                                                   match.qry_start, match.match_len());
                            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                                SeqPro::MaskedSequenceManager> >) {
                                return manager_ptr->getSubSequence(match.qry_chr_index,
                                    match.qry_start, match.match_len());
                            } else {
                                throw std::runtime_error("Unhandled manager type in variant.");
                            }
                        }, query_manager);

                        // 如果是反向链，进行反向互补
                        if (match.strand() == REVERSE) {
                            query_seq = reverseComplement(query_seq);
                        }

                        // 比较序列
                        if (ref_seq == query_seq) {
                            ++correct_matches;
                        } else {
                            ++incorrect_matches;

                            // 线程安全的首错捕获逻辑
                            bool expected = false;
                            if (!diagnostic_dump_done.load(std::memory_order_relaxed) &&
                                diagnostic_dump_done.compare_exchange_strong(expected, true)) {
                                spdlog::warn(
                                    "序列不匹配: ref_chr={}, ref_start={}, query_chr={}, "
                                    "query_start={}, length={}, strand={}\n"
                                    "  Ref Seq:    {}\n"
                                    "  Query Seq{}: {}",
                                    match.ref_chr_index, match.ref_start, match.qry_chr_index,
                                    match.qry_start, match.match_len(),
                                    (match.strand() == FORWARD ? "FORWARD" : "REVERSE"), ref_seq,
                                    (match.strand() == REVERSE ? " (RC)" : ""), query_seq
                                );
                                spdlog::error("--- [CAPTURED FIRST MISMATCH] INITIATING DIAGNOSTIC DUMP ---");

                                // 打印错误匹配的详细信息
                                spdlog::error("Failing Match Details:");
                                spdlog::error("  - Reference: {}:{} (len:{})", match.ref_chr_index,
                                              match.ref_start, match.match_len());
                                spdlog::error("  - Query:     {}:{} (len:{})", match.qry_chr_index,
                                              match.qry_start, match.match_len());
                                spdlog::error("  - Strand:    {}", (match.strand() == FORWARD ? "FORWARD" : "REVERSE"));
                            }
                        }
                    } catch (const std::exception &e) {
                        ++incorrect_matches;
                        spdlog::warn("处理匹配项时发生异常: {}", e.what());
                    }
                }
            }
        } // omp for
    } // omp parallel

    // 将OpenMP reduction结果赋值给结果结构体
    result.total_matches = total_matches;
    result.correct_matches = correct_matches;
    result.incorrect_matches = incorrect_matches;

    // 确保最终进度是100%
    spdlog::info("验证进度: {} / {} (100.00%)", total_work_items, total_work_items);

    spdlog::info("验证完成: 总匹配数={}, 正确匹配数={}, 错误匹配数={}, 正确率={:.2f}%",
                 result.total_matches,
                 result.correct_matches,
                 result.incorrect_matches,
                 result.accuracy());

    return result;
}

// 重载版本：支持 SharedManagerVariant
ValidationResult validateAnchorsCorrectness(
    const MatchVec3DPtr& anchors,
    const SeqPro::SharedManagerVariant& ref_manager,
    const SeqPro::SharedManagerVariant& query_manager
) {
    // 解引用 SharedManagerVariant 并调用原始函数
    return validateAnchorsCorrectness(anchors, *ref_manager, *query_manager);
}


/*------------------------------------------------------------------*/
/*  WaveFront 配置：常规 affine + ends-free 版本                     */
/*------------------------------------------------------------------*/
inline wavefront_aligner_t* newDefaultWFA() {
    static wavefront_aligner_attr_t attr = [] {
        wavefront_aligner_attr_t a = wavefront_aligner_attr_default;
        a.distance_metric = gap_affine;
        a.affine_penalties.mismatch = 2;
        a.affine_penalties.gap_opening = 3;
        a.affine_penalties.gap_extension = 1;
        a.memory_mode = wavefront_memory_ultralow;
        a.heuristic.strategy = wf_heuristic_wfadaptive;
        a.heuristic.min_wavefront_length = 10;
        a.heuristic.max_distance_threshold = 50;
        a.heuristic.steps_between_cutoffs = 1;
        return a;
        }();
    return wavefront_aligner_new(&attr);
}

/* 可指定 ends-free 参数的 WFA（用于左右延伸）*/
inline wavefront_aligner_t* newEndsFreeWFA(
    int pat_begin_free, int pat_end_free,
    int txt_begin_free, int txt_end_free)
{
    wavefront_aligner_attr_t a = wavefront_aligner_attr_default;
    a.distance_metric = gap_affine;
    a.affine_penalties = { 3, 1, 2 };        // 与 newDefaultWFA 相同
    a.memory_mode = wavefront_memory_ultralow;
    a.alignment_form.span = alignment_endsfree;
    a.alignment_form.pattern_begin_free = pat_begin_free;
    a.alignment_form.pattern_end_free = pat_end_free;
    a.alignment_form.text_begin_free = txt_begin_free;
    a.alignment_form.text_end_free = txt_end_free;
    return wavefront_aligner_new(&a);
}

/*------------------------------------------------------------------*/
/*  子串提取工具（支持 MaskedManager）                               */
/*------------------------------------------------------------------*/
inline std::string subSeq(const SeqPro::ManagerVariant& mv,
    ChrIndex chr, Coord_t b, Coord_t l)
{
    return std::visit([&](auto& ptr) -> std::string {
        using P = std::decay_t<decltype(ptr)>;
        if constexpr (std::is_same_v<P,
            std::unique_ptr<SeqPro::SequenceManager>>)
            return ptr->getSubSequence(chr, b, l);
        else
            return ptr->getOriginalManager().getSubSequence(chr, b, l);
        }, mv);
}

/*------------------------------------------------------------------*/
/*     对单个 Anchor 进行左右 ends-free 延伸，更新坐标/CIGAR         */
/*------------------------------------------------------------------*/
static void extendAnchorEnds(Anchor& a,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr,
    bool   left_extend = true,
    bool   right_extend = true,
    Length_t extend_len = 10000)
{
    if (!left_extend && !right_extend) return;

    const bool fwd = (a.strand == FORWARD);
    auto& cigar = a.cigar;
    Length_t ref_chr_length = std::visit(
        [&](auto& ptr) -> SeqPro::Length {
            using P = std::decay_t<decltype(ptr)>;
            if constexpr (std::is_same_v<P,
                std::unique_ptr<SeqPro::SequenceManager>>) {
                // plain SequenceManager → original length
                return ptr->getSequenceLength(a.ref_chr_index);
            }
            else {
                // MaskedSequenceManager → forward to its original manager
                return ptr->getOriginalManager()
                    .getSequenceLength(a.ref_chr_index);
            }
        },
        ref_mgr
    );

    Length_t qry_chr_length = std::visit(
        [&](auto& ptr) -> SeqPro::Length {
            using P = std::decay_t<decltype(ptr)>;
            if constexpr (std::is_same_v<P,
                std::unique_ptr<SeqPro::SequenceManager>>) {
                // plain SequenceManager → original length
                return ptr->getSequenceLength(a.qry_chr_index);
            }
            else {
                // MaskedSequenceManager → forward to its original manager
                return ptr->getOriginalManager()
                    .getSequenceLength(a.qry_chr_index);
            }
        },
        qry_mgr
    );

    /* ---------- 左端延伸 ---------- */
/* ---------- 左端延伸 ---------- */
    if (left_extend) {
        while (true) {
            // 1) 计算这轮要对齐的 ref 区间 [ref_sub_start, ref_sub_start+ref_sub_len)
            bool to_end = false;
            Coord_t ref_sub_start = a.ref_start;
            Length_t ref_sub_len = extend_len;
			if (ref_sub_start < extend_len) {
				// 如果 ref_start 小于 extend_len，说明要从头开始
				ref_sub_len = ref_sub_start;  // 限制在染色体长度内
				ref_sub_start = 0;             // 到头了
                to_end = true;
			}
			else {
				ref_sub_start -= extend_len;   // 向左延伸
			}

            // 2) 计算这轮要对齐的 qry 区间 [qry_sub_start, qry_sub_start+qry_sub_len)
            Coord_t qry_sub_start;
            Length_t qry_sub_len = extend_len;
            
            if (fwd) {
                // + 链：向左延伸就是在 qry_start 之前取窗口
                if (a.qry_start > extend_len) {
                    qry_sub_start = a.qry_start - extend_len;
                }
                else {
                    qry_sub_start = 0;
                    qry_sub_len = a.qry_start;
                    to_end = true; // 到头了
                }
            }
            else {
                // - 链：向左延伸对应 query 上在 (qry_start+qry_len) 之后
                qry_sub_start = a.qry_start + a.qry_len;
                if (qry_sub_start + qry_sub_len > qry_chr_length) {
                    qry_sub_len = qry_chr_length - qry_sub_start;
                    to_end = true;
                }
            }
            if (!ref_sub_len || !qry_sub_len) {
                // 如果子序列长度为 0，直接跳过
                break;
            }
			Length_t sub_len = std::min(ref_sub_len, qry_sub_len);
            // 3) 抽出子序列，revcomp（如果需要）
            std::string ref_sub_seq = subSeq(ref_mgr, a.ref_chr_index,
                ref_sub_start, sub_len);
            std::string qry_sub_seq = subSeq(qry_mgr, a.qry_chr_index,
                qry_sub_start, sub_len);
            reverseSequence(ref_sub_seq);
            if (fwd) {
                reverseSequence(qry_sub_seq);
            }
            else {
                complementSequence(qry_sub_seq);
            }

            // 4) 用 KSW2 延伸
            
            int tmp_ref_len = 0, tmp_qry_len = 0;
            Cigar_t tmp_cigar = extendAlignKSW2(
                ref_sub_seq, qry_sub_seq);

            bool reach_break_len = false;
            Cigar_t final_cigar;
            for (auto& unit : tmp_cigar) {
                uint32_t len;
                char  op;   // 0=M,1=I,2=D,7='=',8='X'
                intToCigar(unit, op, len);

                if (op == 'I') {
                    if (len > 200) {
                        reach_break_len = true;
                        break;
                    }
                    tmp_qry_len += len;

                }
                else if (op == 'D') {
                    if (len > 200) {
                        reach_break_len = true;
                        break;
                    }
                    tmp_ref_len += len;
                }
                else {
                    tmp_ref_len += len;
                    tmp_qry_len += len;
                }
                final_cigar.push_back(unit);
            }

            // 5) 把本轮 tmp_cigar 插到总 cigar 前面
			prependCigar(cigar, final_cigar);

            // 6) 更新 Anchor 的 ref_start/ref_len
			a.ref_start -= tmp_ref_len;
			a.ref_len += tmp_ref_len;

            if (fwd) {
				a.qry_start -= tmp_qry_len;

            }
            a.qry_len += tmp_qry_len;
            // 8) 退出条件：z-drop 或者已经到染色体/查询序列边界
            bool ref_neq = tmp_ref_len != sub_len;
            bool qry_neq = tmp_qry_len != sub_len;
            if (ref_neq || qry_neq || to_end || reach_break_len) {
                break;
            }
        }
    }


    /* ---------- 右端延伸 ---------- */
    if (right_extend) {

        while (true) {

            Coord_t ref_start = a.ref_start + a.ref_len;
            Length_t ref_sub_len = extend_len;
            Coord_t qry_start;
            if (fwd) {
                qry_start = a.qry_start + a.qry_len;
            }
            else {
                qry_start = a.qry_start;
            }
            Length_t qry_sub_len = extend_len;
            bool to_end = false;
           
            if (ref_start + ref_sub_len > ref_chr_length) {
                ref_sub_len = ref_chr_length - ref_start;  // 限制在染色体长度内
				to_end = true;  // 到达染色体末尾
            }

            if (fwd) {
                if (qry_start + qry_sub_len > qry_chr_length) {
                    qry_sub_len = qry_chr_length - qry_start;  // 限制在染色体长度内
                    to_end = true;
                }
            }
            else {  
                if (qry_start < qry_sub_len) {
                    qry_sub_len = qry_start;  // 限制在染色体长度内
                    qry_start = 0;            // 反向链从末尾开始
                    to_end = true;
                }
                else {
                    qry_start -= qry_sub_len;
                }
            }
			if (!ref_sub_len || !qry_sub_len) {
				// 如果子序列长度为 0，直接跳过
				break;
			}
            Length_t sub_len = std::min(ref_sub_len, qry_sub_len);
            std::string ref_sub_seq = subSeq(ref_mgr, a.ref_chr_index, ref_start, sub_len);
            std::string qry_sub_seq = subSeq(qry_mgr, a.qry_chr_index, qry_start, sub_len);
            if (!fwd) {
                reverseComplement(qry_sub_seq);
            }
            
            int tmp_ref_len = 0, tmp_qry_len = 0;
            Cigar_t tmp_cigar = extendAlignKSW2(ref_sub_seq, qry_sub_seq);
            bool reach_break_len = false;
            Cigar_t final_cigar;
            for (auto& unit : tmp_cigar) {
                uint32_t len;
                char  op;   // 0=M,1=I,2=D,7='=',8='X'
                intToCigar(unit, op, len);

                if (op == 'I') {
                    if (len > 200) {
                        reach_break_len = true;
                        break;
                    }
                    tmp_qry_len += len;

                }
                else if (op == 'D') {
                    if (len > 200) {
                        reach_break_len = true;
                        break;
                    }
                    tmp_ref_len += len;
                }
                else {
                    tmp_ref_len += len;
                    tmp_qry_len += len;
                }
				final_cigar.push_back(unit);
            }

            appendCigar(cigar, final_cigar);
            a.ref_len += tmp_ref_len;

            if (!fwd) {
                a.qry_start -= tmp_qry_len;
            }

            a.qry_len += tmp_qry_len;

            bool ref_neq = tmp_ref_len != sub_len;
            bool qry_neq = tmp_qry_len != sub_len;
            if (ref_neq || qry_neq || to_end || reach_break_len) {
                break;
            }

            
        }
        
    }
}

/* ---------- 核心实现：把 MatchCluster 扩展成一组 Anchor ---------- */
static AnchorVec extendClusterCore(const MatchCluster& cluster,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& qry_mgr,
    bool extend_end = true,
    Length_t extend_len = 10000,
    bool   split_large_gap = false)
{
    AnchorVec anchors;
    if (cluster.empty()) return anchors;

    /* === 初始 anchor 状态 === */
    const Match& first = cluster.front();
    const Strand strand = first.strand();
    const bool   fwd = strand == FORWARD;
    const ChrIndex ref_chr = first.ref_chr_index;
    const ChrIndex qry_chr = first.qry_chr_index;

    Coord_t ref_beg = start1(first);
    Coord_t ref_end = ref_beg;
    Coord_t qry_beg = fwd ? start2(first)
        : start2(first) + len2(first);
    Coord_t qry_end = qry_beg;

    Cigar_t cigar; cigar.reserve(cluster.size() * 2);
    Length_t aln_len = 0;
    Length_t aln_match_len = 0;

    auto pushEq = [&](uint32_t len) {
        appendCigarOp(cigar, 'M', len);
        aln_len += len;
        aln_match_len += len;
        ref_end += len;
        fwd ? (qry_end += len) : (qry_end -= len);
        };
    auto flush = [&] {
        Anchor a;
        if (fwd)
            a = Anchor(ref_chr, ref_beg, ref_end - ref_beg,
                qry_chr, qry_beg, qry_end - qry_beg,
                strand, aln_len, aln_match_len, std::move(cigar));
        else
            a = Anchor(ref_chr, ref_beg, ref_end - ref_beg,
                qry_chr, qry_end, qry_beg - qry_end,
                strand, aln_len, aln_match_len, std::move(cigar));
            
        anchors.emplace_back(a);
        /* reset state for next segment */
        cigar.clear(); cigar.shrink_to_fit(); cigar.reserve(16);
        aln_len = 0;
        aln_match_len = 0;
        ref_beg = ref_end;
        qry_beg = qry_end;
        };

    /* === WaveFront 对齐器 === */
    std::unique_ptr<wavefront_aligner_t, decltype(&wavefront_aligner_delete)>
        wf(newDefaultWFA(), wavefront_aligner_delete);

    /* ========== 遍历 cluster ========== */
    for (size_t i = 0; i < cluster.size(); ++i) {
        const Match& m = cluster[i];
        pushEq(len1(m));                         // 精确 match

        if (i + 1 == cluster.size()) break;

        /* ---- 计算 gap ---- */
        const Match& nxt = cluster[i + 1];
        const Coord_t rgBeg = start1(m) + len1(m);
        const Coord_t rgEnd = start1(nxt);
        Coord_t qgBeg = fwd ? start2(m) + len2(m)
            : start2(nxt) + len2(nxt);
        Coord_t qgEnd = fwd ? start2(nxt)
            : start2(m);

        const uint32_t rgLen = rgEnd > rgBeg ? rgEnd - rgBeg : 0;
        const uint32_t qgLen = qgEnd > qgBeg ? qgEnd - qgBeg : 0;
        if (rgLen == 0 && qgLen == 0) continue;

        Cigar_t buf, gap;
        if (rgLen == 0 || qgLen == 0) {
            /* 纯 I / 纯 D */
            gap.push_back(cigarToInt(rgLen ? 'D' : 'I',
                rgLen ? rgLen : qgLen));
        }
        else {
            /* 需要对齐的 Indel 片段 */
            std::string ref_gap = subSeq(ref_mgr, ref_chr, rgBeg, rgLen);
            std::string qry_gap = subSeq(qry_mgr, qry_chr, qgBeg, qgLen);
            if (!fwd) reverseComplement(qry_gap);

            const auto Lt = ref_gap.size();
            const auto Lq = qry_gap.size();
            const auto rho = std::abs((int64_t)Lt - (int64_t)Lq) /
                double(std::min(Lt, Lq));

            if (rho <= 0.3 && Lt > 10 && Lq > 10) {   // 用 WFA
                int cigar_len; uint32_t* tmp;
                wavefront_align(wf.get(), ref_gap.c_str(), Lt,
                    qry_gap.c_str(), Lq);
                cigar_get_CIGAR(wf->cigar, false, &tmp, &cigar_len);
                gap.assign(tmp, tmp + cigar_len);
            }
            else {                                  // KSW2 全局对齐
                gap = globalAlignKSW2(ref_gap, qry_gap);
            }
        }

        /* ---- 扫描 gap-cigar ---- */
        buf.reserve(gap.size());
        for (auto unit : gap) {
            uint32_t len = unit >> 4;
            uint8_t  op = unit & 0xf;   // 0=M,1=I,2=D,7='=',8='X'
            bool big = split_large_gap &&
                (op == 1 || op == 2) && len > 50;

            if (big) {
                if (!buf.empty()) { appendCigar(cigar, buf); buf.clear(); }
                flush();

                /* 纯大 Gap：移动起点 */
                if (op == 1) fwd ? (qry_beg += len) : (qry_beg -= len);
                else         ref_beg += len;

                ref_end = ref_beg;
                qry_end = qry_beg;
            }
            else {
                buf.push_back(unit);
                /* 更新末端坐标 */
                switch (op) {
                case 1: fwd ? (qry_end += len) : (qry_end -= len); break;
                case 2: ref_end += len;                            break;
                default: ref_end += len;
                    fwd ? (qry_end += len) : (qry_end -= len);
                }
                if (op != 3) {
                    aln_len += len;
                    aln_match_len += len;
                }   // op==3 (N) 不出现
            }
        }
        if (!buf.empty()) appendCigar(cigar, buf);
    }
    flush();
    if (extend_end) {
        extendAnchorEnds(anchors.front(), ref_mgr, qry_mgr, true, false, extend_len);
        extendAnchorEnds(anchors.back(), ref_mgr, qry_mgr, false, true, extend_len);
    }
    return anchors;
}
/*------------------------------------------------------------------*/
/*  对外包装（老接口 + 新接口）                                      */
/*------------------------------------------------------------------*/
AnchorVec extendClusterToAnchorVec(const MatchCluster& c,
    const SeqPro::ManagerVariant& r,
    const SeqPro::ManagerVariant& q,
    bool extend_end,
    Length_t extend_len)
{
    return extendClusterCore(c, r, q,
        extend_end, extend_len, true);
}

Anchor extendClusterToAnchor(const MatchCluster& c,
    const SeqPro::ManagerVariant& r,
    const SeqPro::ManagerVariant& q,
    bool extend_end,
    Length_t extend_len)
{
    AnchorVec v = extendClusterCore(c, r, q,
        extend_end, extend_len, false);
    return v.empty() ? Anchor() : std::move(v.front());
}


/// 同时验证 ref/query，两者都通过才算成功
void validateClusters(const ClusterVecPtrByStrandByQueryRefPtr& cluster_vec_ptr)
{
    if (!cluster_vec_ptr) {
        spdlog::debug("validateClusters: cluster_vec_ptr is null, nothing to check.");
        return;
    }

    std::size_t total_clusters = 0;
    std::size_t failed_clusters = 0;

    bool reverse_cluster = false;

    for (std::size_t strand_i = 0; strand_i < cluster_vec_ptr->size(); ++strand_i) {
        const auto& by_query = (*cluster_vec_ptr)[strand_i];

        for (std::size_t q_i = 0; q_i < by_query.size(); ++q_i) {
            const auto& by_ref = by_query[q_i];

            for (std::size_t r_i = 0; r_i < by_ref.size(); ++r_i) {
                const auto& clusters_ptr = by_ref[r_i];
                if (!clusters_ptr) continue;

                for (std::size_t c_i = 0; c_i < clusters_ptr->size(); ++c_i) {
                    ++total_clusters;
                    const MatchCluster& cluster = (*clusters_ptr)[c_i];
 
                    Strand strand = cluster[0].strand();

					if (strand == REVERSE && cluster.size() > 5) {
						reverse_cluster = true;
					}   

                    for (std::size_t m_i = 0; m_i < cluster.size(); ++m_i) {
                        if (cluster[m_i].strand() != strand) {
                            spdlog::debug(
                                "❌ Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                "strand mismatch (expected {}, got {})",
                                strand_i, q_i, r_i, c_i, static_cast<int>(strand), static_cast<int>(cluster[m_i].strand()));
                        }
                    }


                    Coord_t ref_last_end = std::numeric_limits<Coord_t>::min();
                    Coord_t qry_last_pos;                 // 启动值依赖方向
					if (strand == FORWARD) {
						qry_last_pos = std::numeric_limits<Coord_t>::min();
					}
					else { // REVERSE
						qry_last_pos = std::numeric_limits<Coord_t>::max();
					}
                    

                    for (std::size_t m_i = 0; m_i < cluster.size(); ++m_i) {
                        const Match& m = cluster[m_i];

                        //-------------------//
                        // 1) 参考坐标检查
                        //-------------------//
                        Coord_t ref_start = m.ref_start;
                        Coord_t ref_end = ref_start + m.match_len(); // 右开

                        if (ref_start < ref_last_end) {
                            spdlog::debug(
                                "❌ Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                "ref_start < ref_last_end ({} < {})",
                                strand_i, q_i, r_i, c_i, ref_start, ref_last_end);
                        }
                            
                        ref_last_end =  ref_end;

                        //-------------------//
                        // 2) 查询坐标检查
                        //-------------------//
                        Coord_t qry_start = m.qry_start;
                        Coord_t qry_end = qry_start + m.match_len();


                        
                        if (m.strand() == FORWARD) {
                            // 正向：要求升序且不重叠
                            if (qry_start < qry_last_pos) {
                                spdlog::debug(
                                    "❌ Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                    "qry_start < qry_last_pos ({} < {})",
                                    strand_i, q_i, r_i, c_i, qry_start, qry_last_pos);
                            }
                            qry_last_pos = qry_end;
                        }
                        else { // REVERSE
                            // 反向：查询坐标应当递减，区间不重叠
                            if (qry_end > qry_last_pos) {
                                spdlog::debug(
                                    "❌ Cluster FAILED (strand={},query={},ref={},cluster={}): "
                                    "qry_start > qry_last_pos ({} > {})",
                                    strand_i, q_i, r_i, c_i, qry_start, qry_last_pos);
                            }
                            qry_last_pos = qry_start;
                        }
                        
                    } // end matches loop


                }
            }
        }
    }
    if (reverse_cluster == false) {
		spdlog::debug("reverse chain may fail");
		return;
    }

    spdlog::debug("validateClusters finished: {} clusters checked, {} failed.",
        total_clusters, failed_clusters);
}

