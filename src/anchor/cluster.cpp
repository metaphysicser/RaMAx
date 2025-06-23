#include "anchor.h"
#include "data_process.h"

// ────────────────────────────────────────────
//  Constructors / reset
// ────────────────────────────────────────────
UnionFind::UnionFind(std::size_t n)
{
    reset(n);
}

void UnionFind::reset(std::size_t n)
{
    parent_.assign(n, -1);                   // -1 表示单元素集合
    component_cnt_ = static_cast<int_t>(n);    // 初始有 n 个独立集合
}

// ────────────────────────────────────────────
//  Query
// ────────────────────────────────────────────
int_t UnionFind::find(int_t x)
{
    // 路径折半（迭代版）
    while (parent_[x] >= 0 && parent_[parent_[x]] >= 0)
    {
        parent_[x] = parent_[parent_[x]];
        x = parent_[x];
    }
    return (parent_[x] < 0) ? x : parent_[x];
}

int_t UnionFind::set_size(int_t x)
{
    return -parent_[find(x)];
}

bool UnionFind::same(int_t a, int_t b)
{
    return find(a) == find(b);
}

// ────────────────────────────────────────────
//  Modification
// ────────────────────────────────────────────
bool UnionFind::unite(int_t a, int_t b)
{
    a = find(a);
    b = find(b);
    if (a == b) return false;

    // 按大小合并：把小树挂到大树
    if (parent_[a] > parent_[b]) std::swap(a, b);

    parent_[a] += parent_[b];  // 更新根 a 的大小（负值）
    parent_[b] = a;           // b 挂到 a
    --component_cnt_;
    return true;
}


// 1. 根据 max_gap / diagdiff / diagfactor 把 unique_match 聚成若干簇
MatchClusterVec buildClusters(MatchVec& unique_match,
    int_t      max_gap,
    int_t      diagdiff,
    double     diagfactor)
{
    const uint_t N = static_cast<uint_t>(unique_match.size());
    MatchClusterVec clusters;
    if (N < 2) return clusters;

    UnionFind uf(N);
    for (uint_t i = 0; i < N; ++i) {
        uint_t i_end = start2(unique_match[i]) + len2(unique_match[i]);
        int_t  i_diag = diag(unique_match[i]);
        for (uint_t j = i + 1; j < N; ++j) {
            int_t sep = start2(unique_match[j]) - i_end;
            if (sep > static_cast<int_t>(max_gap)) break;
            int_t diag_diff = std::abs(diag(unique_match[j]) - i_diag);
            int_t th = std::max(diagdiff, static_cast<int_t>(diagfactor * sep));
            if (diag_diff <= th) uf.unite(i, j);
        }
    }

    // 根节点映射 → cluster_id
    std::vector<int_t> root_map(N, -1);
    for (uint_t idx = 0; idx < N; ++idx) {
        int_t root = uf.find(idx);
        int_t& cid = root_map[root];
        if (cid == -1) {
            cid = static_cast<int_t>(clusters.size());
            clusters.emplace_back();
        }
        clusters[cid].push_back(std::move(unique_match[idx]));
        // clusters[cid].push_back(unique_match[idx]);
    }
    return clusters;
}

ClusterVecPtrByStrandByQueryRefPtr
clusterAllChrMatch(const MatchByStrandByQueryRefPtr& unique_anchors,
    const MatchByStrandByQueryRefPtr& repeat_anchors,
    ThreadPool& pool)
{
    // ---------- 0. 判空 ----------
    if (!unique_anchors || !repeat_anchors) {
        spdlog::warn("[clusterAllChrAnchors] empty anchor ptr");
        return std::make_shared<ClusterVecPtrByStrandByQueryRef>();   // 返回空指针容器
    }

    // ---------- 1. 构建结果桶 ----------
    auto cluster_ptr = std::make_shared<ClusterVecPtrByStrandByQueryRef>();
    cluster_ptr->resize(2);  // strand 维度

    for (auto& query_layer : *cluster_ptr) {
        query_layer.resize(unique_anchors->front().size());           // query-chr
        for (auto& ref_row : query_layer)
            ref_row.resize(unique_anchors->front().front().size());   // ref-chr
    }

    // ---------- 2. 提交任务 ----------
    using ClusterFuture = std::future<std::shared_ptr<MatchClusterVec>>;
    std::vector<ClusterFuture> futures;

    for (uint_t k = 0; k < 2; ++k) {
        for (uint_t i = 0; i < (*unique_anchors)[k].size(); ++i) {
            for (uint_t j = 0; j < (*unique_anchors)[k][i].size(); ++j) {
                
                MatchVec& uniq = (*unique_anchors)[k][i][j];
                MatchVec& rept = (*repeat_anchors)[k][i][j];

                futures.emplace_back(
    pool.enqueue(
        [p_uniq = &uniq, p_rept = &rept]() -> std::shared_ptr<MatchClusterVec> {
            return clusterChrMatch(*p_uniq, *p_rept);
        }
    )
);
            }
        }
    }

    // ---------- 3. 收集结果 ----------
    size_t idx = 0;
    for (uint_t k = 0; k < 2; ++k) {
        for (uint_t i = 0; i < (*unique_anchors)[k].size(); ++i) {
            for (uint_t j = 0; j < (*unique_anchors)[k][i].size(); ++j) {
                (*cluster_ptr)[k][i][j] = futures[idx++].get();
            }
        }
    }

    return cluster_ptr;   // 返回 shared_ptr
}

/*!
 * \brief  把 src（三维：strand × query × ref）按 ref 合并
 *         每个 ref 得到一个合并后的 MatchClusterVecPtr
 *
 * \param  src  输入三维结构（shared_ptr 外围保持不变）
 * \return      shared_ptr< vector< MatchClusterVecPtr > >
 */
ClusterVecPtrByRefPtr
groupClustersByRef(const ClusterVecPtrByStrandByQueryRefPtr& src)
{
    /* ---------- 0. 早退 ---------- */
    auto by_ref = std::make_shared<ClusterVecPtrByRef>();   // 最终返回对象
    if (!src || src->empty() || (*src)[0].empty()) return by_ref;

    /* ---------- 1. 基础信息 ---------- */
    const size_t ref_n = (*src)[0][0].size();               // 所有行列数一致

    /* ---------- 2. 统计每个 ref 的簇总量 & 记录首矢量指针 ---------- */
    std::vector<size_t> cluster_counts(ref_n, 0);
    std::vector<MatchClusterVec*> first_src_vec(ref_n, nullptr);

    for (auto& query_layer : *src)          // strand
        for (auto& ref_row : query_layer)   // query
            for (size_t r = 0; r < ref_n; ++r) {
                auto& vec = *ref_row[r];
                cluster_counts[r] += vec.size();
                if (!first_src_vec[r] && !vec.empty())
                    first_src_vec[r] = &vec;             // 只记录第一次出现
            }

    /* ---------- 3. 构造目标数组，首矢量直接 swap ---------- */
    by_ref->resize(ref_n);
    std::vector<bool> swapped(ref_n, false);

    for (size_t r = 0; r < ref_n; ++r) {
        auto dst = std::make_shared<MatchClusterVec>();
        if (first_src_vec[r]) {                          // 有首矢量就直接搬进来
            dst->swap(*first_src_vec[r]);                // 零拷贝
            swapped[r] = true;
        }
        dst->reserve(cluster_counts[r]);                 // 再保证容量充足
        (*by_ref)[r] = std::move(dst);
    }

    /* ---------- 4. 把其余簇 append 进去 ---------- */
    for (auto& query_layer : *src)
        for (auto& ref_row : query_layer)
            for (size_t r = 0; r < ref_n; ++r) {
                auto& src_vec = *ref_row[r];
                if (src_vec.empty()) continue;

                auto& dst = *(*by_ref)[r];
                if (swapped[r]) {
                    dst.insert(dst.end(),
                        std::make_move_iterator(src_vec.begin()),
                        std::make_move_iterator(src_vec.end()));
                }
                else {                                 // 之前全空，直接 swap
                    dst.swap(src_vec);
                    swapped[r] = true;
                }
            }

    return by_ref;                                       // shared_ptr 返回
}

MatchClusterVecPtr
groupClustersToVec(const ClusterVecPtrByStrandByQueryRefPtr& src,
    ThreadPool& pool, uint_t thread_num)
{
    /* ---------- 0. 早退 ---------- */
    auto result = std::make_shared<MatchClusterVec>();
    if (!src || src->empty() || (*src)[0].empty()) return result;

    /* ---------- 1. 收集每行信息并统计总量 ---------- */
    struct RowInfo { size_t si, qi, count; };
    std::vector<RowInfo> rows;
    uint_t total_clusters = 0;

    const uint_t strand_n = src->size();
    for (uint_t si = 0; si < strand_n; ++si) {
        auto& strand_layer = (*src)[si];
        const uint_t query_n = strand_layer.size();
        for (uint_t qi = 0; qi < query_n; ++qi) {
            uint_t cnt = 0;
            for (auto& ref_vec_ptr : strand_layer[qi])
                if (ref_vec_ptr) cnt += ref_vec_ptr->size();
            if (cnt) {
                rows.push_back({ si, qi, cnt });
                total_clusters += cnt;
            }
        }
    }
    if (total_clusters == 0) return result;

    /* ---------- 2. 划分任务（顺序装箱） ---------- */

    const uint_t avg = (total_clusters + thread_num - 1) / thread_num;   // 每桶理想负载
    std::vector<std::vector<RowInfo>> buckets(thread_num);

    size_t cur_bucket = 0, bucket_load = 0;
    for (const auto& r : rows) {
        if (bucket_load >= avg && cur_bucket + 1 < thread_num) {
            ++cur_bucket;
            bucket_load = 0;
        }
        buckets[cur_bucket].push_back(r);
        bucket_load += r.count;
    }

    /* ---------- 3. 预分配结果向量 & 原子写指针 ---------- */
    result->resize(total_clusters);
    std::atomic<size_t> pos{ 0 };

    /* ---------- 4. 并行执行 ---------- */
    for (const auto& todo : buckets) {
        if (todo.empty()) continue;
        pool.enqueue([src, result, &pos, todo]() {
            for (const auto& row : todo) {
                auto& query_row = (*src)[row.si][row.qi];
                for (auto& ref_vec_ptr : query_row) {
                    if (!ref_vec_ptr || ref_vec_ptr->empty()) continue;
                    auto& src_vec = *ref_vec_ptr;            // MatchClusterVec
                    for (auto& cluster : src_vec) {
                        size_t idx = pos.fetch_add(1, std::memory_order_relaxed);
                        (*result)[idx] = std::move(cluster); // 顺序无关
                    }
                    MatchClusterVec().swap(src_vec);         // 释放（若需保留删掉此行）
                }
            }
            });
    }

    /* ---------- 5. 等待全部完成 ---------- */
    pool.waitAllTasksDone();
    // assert(pos == total_clusters);

    return result;
}



// 2. 给定一个 cluster，返回其最佳非交叉链（DP O(N^2)）
MatchVec bestChainDP(MatchVec& cluster, double diagfactor)
{
    if (cluster.empty()) return {};
    if (cluster.size() == 1) return std::move(cluster);

    std::sort(cluster.begin(), cluster.end(),
        [](const Match& a, const Match& b) { return start2(a) < start2(b); });

    const uint_t N = static_cast<uint_t>(cluster.size());
    std::vector<int_t> score(N), pred(N, -1);
    uint_t best_idx = 0;

    for (uint_t i = 0; i < N; ++i) {
        score[i] = len2(cluster[i]);
        for (uint_t j = 0; j < i; ++j) {
            if (start1(cluster[i]) <= start1(cluster[j])) continue;
            int_t prev_end2 = start2(cluster[j]) + len2(cluster[j]);
            if (start2(cluster[i]) <= prev_end2) continue;
            int_t sep = start2(cluster[i]) - prev_end2;
            int_t d = std::abs(diag(cluster[i]) - diag(cluster[j]));
            int_t cand = score[j] + len2(cluster[i]) - (sep + static_cast<int_t>(diagfactor * d));
            if (cand > score[i]) { score[i] = cand; pred[i] = static_cast<int_t>(j); }
        }
        if (score[i] > score[best_idx]) best_idx = i;
    }

    MatchVec chain;
    for (int_t k = static_cast<int_t>(best_idx); k != -1; k = pred[k])
        chain.emplace_back(std::move(cluster[k]));
    std::reverse(chain.begin(), chain.end());
    return chain;
}

/* ───────────────────────────────────────────────────────── *
 * 对外主函数                                              *
 * ───────────────────────────────────────────────────────── */
MatchClusterVecPtr clusterChrMatch(MatchVec& unique_match, MatchVec& repeat_match, int_t max_gap, int_t diagdiff, double diagfactor, int_t min_cluster_length)
{

    auto best_chain_clusters = std::make_shared<MatchClusterVec>();
    if (unique_match.size() == 0) {
        return best_chain_clusters;
    }
    // 1. 聚簇
    MatchClusterVec clusters = buildClusters(unique_match, max_gap, diagdiff, diagfactor);

    // 2. 每簇提链 + 长度过滤
    
    best_chain_clusters->reserve(clusters.size());

    for (auto& cluster : clusters) {
        if (cluster.empty()) continue;

        MatchVec best_chain = bestChainDP(cluster, diagfactor);
        if (best_chain.empty()) { releaseCluster(cluster); continue; }

        int_t span = start1(best_chain.back()) + len1(best_chain.back()) - start1(best_chain.front());
        if (span >= min_cluster_length) {
            best_chain_clusters->emplace_back(std::move(best_chain));
        }
        else {
            std::move(best_chain.begin(), best_chain.end(), std::back_inserter(repeat_match));
        }
        releaseCluster(cluster);     // 回收 cluster 剩余元素
    }

    return best_chain_clusters;
}

 /*!
  * \brief  把 cl 依照冲突矩形切成 ≤2 段
  *         · 仅当 ref_hit 为 true 时才检查 [bad_r_beg,bad_r_end)
  *         · 仅当 query_hit 为 true 时才检查 [bad_q_beg,bad_q_end)
  *
  * \return 0 段：整簇都冲突；1 段：完全无冲突；2 段：被矩形切成前/后两段
  */
MatchClusterVec
splitCluster(const MatchCluster& cl,
    bool  ref_hit,
    int_t bad_r_beg, int_t bad_r_end,
    bool  query_hit,
    int_t bad_q_beg, int_t bad_q_end)
{
    MatchClusterVec parts;
    if (cl.empty()) return parts;

    /* ---------- 1. λ 判断“命中” ---------- */
    auto hit = [&](const Match& m) -> bool
        {
            bool r_overlap = false;   
            bool q_overlap = false;

            if (ref_hit) {
                int_t rb = start1(m), re = rb + len1(m);
                r_overlap = !(re <= bad_r_beg || rb >= bad_r_end);
            }
            if (query_hit) {
                int_t qb = start2(m), qe = qb + len2(m);
                q_overlap = !(qe <= bad_q_beg || qb >= bad_q_end);
            }

            /* 两个维度都要满足要检查的那部分 */
            return r_overlap || q_overlap;
        };

    /* ---------- 2. 找第一次 / 最后一次命中 ---------- */
    size_t first_hit = cl.size();   // 未命中
    size_t last_hit = 0;

    for (size_t i = 0; i < cl.size(); ++i)
        if (hit(cl[i])) {
            first_hit = std::min(first_hit, i);
            last_hit = i;
        }

    /* ---------- 3. 根据命中情况切分 ---------- */
    if (first_hit == cl.size()) {              // 完全不冲突
        parts.emplace_back(cl);
        return parts;
    }

    if (first_hit > 0)                         // 前段
        parts.emplace_back(cl.begin(), cl.begin() + first_hit);

    if (last_hit + 1 < cl.size())              // 后段
        parts.emplace_back(cl.begin() + last_hit + 1, cl.end());

    return parts;                              // 0、1 或 2 段
}

/*------------------------------------------------------------------*
 *  by_ref        : shared_ptr<vector<MatchClusterVecPtr>>  (ref 维)
 *  返回值 by_refQ : shared_ptr<vector< vector<MatchClusterVecPtr> >>
 *                   └──ref──┘└────query────┘
 *------------------------------------------------------------------*/
ClusterVecPtrByRefQueryPtr
groupClustersByRefQuery(const ClusterVecPtrByRefPtr& by_ref,
    SeqPro::ManagerVariant& query_fasta_manager,
    ThreadPool& pool)
{
    /* ---------- 0. 提前创建返回对象 ---------- */
    auto by_ref_query = std::make_shared<ClusterVecPtrByRefQuery>();
    if (!by_ref || by_ref->empty()) return by_ref_query;

    /* ---------- 1. 查询染色体总数 ---------- */
    const uint32_t query_cnt = std::visit([](auto&& mgr) {
        return static_cast<uint32_t>(mgr->getSequenceCount());
        }, query_fasta_manager);

    /* ---------- 2. 预分配 [ref][query] 矩阵 ---------- */
    const uint32_t ref_cnt = static_cast<uint32_t>(by_ref->size());
    by_ref_query->resize(ref_cnt);
    for (uint32_t r = 0; r < ref_cnt; ++r)
        (*by_ref_query)[r].resize(query_cnt);          // 初值 = 空 vector

    /* ---------- 3. 每个 ref 开一个线程任务 ---------- */
    for (uint32_t r = 0; r < ref_cnt; ++r) {
        auto src_ptr = (*by_ref)[r];                   // 可能为空
        if (!src_ptr || src_ptr->empty()) continue;

        /* 捕获：src_ptr 按值、by_ref_query 按 shared_ptr 值 */
        pool.enqueue([src_ptr,
            by_ref_query,
            &query_fasta_manager,
            r, query_cnt]
            {
                /* 指向目标行，避免 lambda 内多级索引 */
                ClusterVecPtrByRef* dst_row = &(*by_ref_query)[r];

                for (auto& cluster : *src_ptr) {
                    if (cluster.empty()) continue;

                    /* ---- 3-1 计算 query 染色体 ID ---- */
                    const auto& m = cluster.front();
                    uint32_t qIdx = std::visit([&](auto&& mgr) {
                        return static_cast<uint32_t>(
                            mgr->getSequenceId(m.query_region.chr_name));
                        }, query_fasta_manager);

                    if (qIdx >= query_cnt) continue;      // 无效或未收录

                    /* ---- 3-2 若目标指针为空就创建 ---- */
                    if (!dst_row->at(qIdx))
                        dst_row->at(qIdx) = std::make_shared<MatchClusterVec>();

                    /* ---- 3-3 把簇移动到目标桶 ---- */
                    dst_row->at(qIdx)->emplace_back(std::move(cluster));
                }
                /* src_ptr 中元素已 move，src_ptr 本身随后被析构 */
            });
    }

    /* ---------- 4. 同步，保证数据就绪 ---------- */
    return by_ref_query;
}
