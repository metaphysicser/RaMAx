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
                                return manager_ptr->getSubSequence(match.ref_region.chr_name, match.ref_region.start,
                                                                   match.ref_region.length);
                            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                                SeqPro::MaskedSequenceManager> >) {
                                return manager_ptr->getSubSequence(match.ref_region.chr_name, match.ref_region.start,
                                                                   match.ref_region.length);
                            } else {
                                throw std::runtime_error("Unhandled manager type in variant.");
                            }
                        }, ref_manager);
                        
                        // 提取查询序列
                        std::string query_seq = std::visit([&match](auto &&manager_ptr) -> std::string {
                            using PtrType = std::decay_t<decltype(manager_ptr)>;
                            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
                                return manager_ptr->getSubSequence(match.query_region.chr_name,
                                                                   match.query_region.start, match.query_region.length);
                            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
                                SeqPro::MaskedSequenceManager> >) {
                                return manager_ptr->getSubSequence(match.query_region.chr_name,
                                                                   match.query_region.start, match.query_region.length);
                            } else {
                                throw std::runtime_error("Unhandled manager type in variant.");
                            }
                        }, query_manager);

                        // 如果是反向链，进行反向互补
                        if (match.strand == REVERSE) {
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
                                    match.ref_region.chr_name, match.ref_region.start, match.query_region.chr_name,
                                    match.query_region.start, match.ref_region.length,
                                    (match.strand == FORWARD ? "FORWARD" : "REVERSE"), ref_seq,
                                    (match.strand == REVERSE ? " (RC)" : ""), query_seq
                                );
                                spdlog::error("--- [CAPTURED FIRST MISMATCH] INITIATING DIAGNOSTIC DUMP ---");

                                // 打印错误匹配的详细信息
                                spdlog::error("Failing Match Details:");
                                spdlog::error("  - Reference: {}:{} (len:{})", match.ref_region.chr_name,
                                              match.ref_region.start, match.ref_region.length);
                                spdlog::error("  - Query:     {}:{} (len:{})", match.query_region.chr_name,
                                              match.query_region.start, match.query_region.length);
                                spdlog::error("  - Strand:    {}", (match.strand == FORWARD ? "FORWARD" : "REVERSE"));
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


AnchorVec extendClusterToAnchor(const MatchCluster& cluster,
    const SeqPro::ManagerVariant& ref_mgr,
    const SeqPro::ManagerVariant& query_mgr,
    const KSW2AlignConfig& cfg)
{
    AnchorVec anchors;
    if (cluster.empty()) return anchors;

    // -- 快速 slice 提取：visit 一次，避免重复 λ 创建 --
    auto subSeq = [&](const SeqPro::ManagerVariant& mv,
        const ChrName& chr, Coord_t b, Coord_t l)->std::string {
            return std::visit([&](auto& p) { return p->getSubSequence(chr, b, l);}, mv);
        };

    /* ===== 初始 anchor 状态 ===== */
    const Match& first = cluster.front();
    Strand strand = first.strand;
    ChrName ref_chr = first.ref_region.chr_name;
    ChrName qry_chr = first.query_region.chr_name;

    Coord_t ref_beg = start1(first);
    Coord_t qry_beg = start2(first);
    Coord_t ref_end = ref_beg;
    Coord_t qry_end = qry_beg;

    Cigar_t cigar; cigar.reserve(cluster.size() * 2);  // 预估
    Coord_t aln_len = 0;

    auto pushEq = [&](uint32_t len) {
        appendCigarOp(cigar, '=', len);
        aln_len += len;  ref_end += len;  qry_end += len;
        };
    auto flush = [&] {
        Region rR{ ref_chr, ref_beg, ref_end - ref_beg };
        Region qR{ qry_chr, qry_beg, qry_end - qry_beg };
        anchors.emplace_back(Match{ rR,qR,strand }, aln_len, std::move(cigar));
        cigar.clear(); cigar.shrink_to_fit(); cigar.reserve(16);
        aln_len = 0;
        };

    /* ==================== 遍历 cluster ==================== */
    for (size_t i = 0;i < cluster.size();++i) {
        const Match& m = cluster[i];
        pushEq(len1(m));                                  // 精确 match

        if (i + 1 == cluster.size()) break;

        const Match& nxt = cluster[i + 1];
        Coord_t rgBeg = start1(m) + len1(m), rgEnd = start1(nxt);
        Coord_t qgBeg = start2(m) + len2(m), qgEnd = start2(nxt);

        // 无 gap
        if (rgBeg >= rgEnd && qgBeg >= qgEnd) continue;

        // 2) 获取 gap 片段
        std::string ref_gap = subSeq(ref_mgr, ref_chr, rgBeg, rgEnd - rgBeg);
        std::string qry_gap = subSeq(query_mgr, qry_chr, qgBeg, qgEnd - qgBeg);
        if (strand == REVERSE) reverseComplement(qry_gap);

        Cigar_t gap = globalAlignKSW2(ref_gap, qry_gap, cfg);   // 调一次 KSW2

        /* ---- 扫描 gap-cigar，遇 >50bp I/D 即分段 ---- */
        Cigar_t buf; buf.reserve(gap.size());
        for (auto unit : gap) {
            uint32_t len = unit >> 4;
            uint8_t  op = unit & 0xf;          // 0=M,1=I,2=D,7='=',8='X'
            bool big = ((op == 1 || op == 2) && len > 50);

            if (big) {
                // 先把已有片段 merge
                if (!buf.empty()) { appendCigar(cigar, buf); buf.clear(); }
                flush();                        // 输出 anchor
                // 移动起点：I 影响 query，D 影响 ref
                if (op == 1) qry_beg += len;
                else       ref_beg += len;
                ref_end = ref_beg; qry_end = qry_beg;
            }
            else {
                buf.push_back(unit);
                // 更新末端坐标
                if (op == 1)            qry_end += len;
                else if (op == 2)       ref_end += len;
                else { ref_end += len; qry_end += len; }
                if (op != 3) aln_len += len;     // 3(N)不会出现
            }
        }
        if (!buf.empty()) appendCigar(cigar, buf);
    }
    flush();                // 收尾
    return anchors;
}
