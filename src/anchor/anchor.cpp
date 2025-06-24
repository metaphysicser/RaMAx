#include "anchor.h"

#include <SeqPro.h>

#include "data_process.h"

//--------------------------------------------------------------------
// 适配 4‑D 结构的聚类分发
//--------------------------------------------------------------------
void groupMatchByQueryRef(MatchVec3DPtr& anchors,
    MatchByStrandByQueryRefPtr unique_anchors,
    MatchByStrandByQueryRefPtr repeat_anchors,
    SeqPro::ManagerVariant& ref_fasta_manager,
    SeqPro::ManagerVariant& query_fasta_manager,
	ThreadPool& pool)
{
    //------------------------------------------------------------
    // 0) 初始化输出矩阵 [strand][query][ref]
    //------------------------------------------------------------
    constexpr uint_t STRAND_CNT = 2;                         // 0=FWD,1=REV
    const uint_t ref_chr_cnt = std::visit([](auto&& manager) {
        using PtrType = std::decay_t<decltype(manager)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager->getSequenceCount();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager->getSequenceCount();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, ref_fasta_manager);

    const uint_t query_chr_cnt = std::visit([](auto&& manager) {
        using PtrType = std::decay_t<decltype(manager)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager->getSequenceCount();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager->getSequenceCount();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, query_fasta_manager);
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

            uint_t qIdx = std::visit([&first](auto&& manager) {
                using PtrType = std::decay_t<decltype(manager)>;
                if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                    return manager->getSequenceId(first.query_region.chr_name);
                } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    return manager->getSequenceId(first.query_region.chr_name);
                } else {
                    throw std::runtime_error("Unhandled manager type in variant.");
                }
            }, query_fasta_manager);

            if(qIdx == SeqPro::SequenceIndex::INVALID_ID) return;

            for (auto& vec : slice)
            {
                if (vec.empty()) continue;

                uint_t rIdx = std::visit([&vec](auto&& manager) {
                    using PtrType = std::decay_t<decltype(manager)>;
                    if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                        return manager->getSequenceId(vec.front().ref_region.chr_name);
                    } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                        return manager->getSequenceId(vec.front().ref_region.chr_name);
                    } else {
                        throw std::runtime_error("Unhandled manager type in variant.");
                    }
                }, ref_fasta_manager);

                if(rIdx == SeqPro::SequenceIndex::INVALID_ID) continue;

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
                                return a.query_region.start < b.query_region.start;
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
                        return a.query_region.start < b.query_region.start;
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







