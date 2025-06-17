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
        query_layer.resize(unique_anchors->size());           // query-chr
        for (auto& ref_row : query_layer)
            ref_row.resize(unique_anchors->front().size());   // ref-chr
    }

    // ---------- 2. 提交任务 ----------
    using ClusterFuture = std::future<std::shared_ptr<MatchClusterVec>>;
    std::vector<ClusterFuture> futures;

    for (uint_t k = 0; k < 2; ++k) {
        for (uint_t i = 0; i < unique_anchors->size(); ++i) {
            for (uint_t j = 0; j < (*unique_anchors)[i].size(); ++j) {

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
        for (uint_t i = 0; i < unique_anchors->size(); ++i) {
            for (uint_t j = 0; j < (*unique_anchors)[i].size(); ++j) {
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
inline ClusterVecPtrByRefPtr
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
    // 1. 聚簇
    MatchClusterVec clusters = buildClusters(unique_match, max_gap, diagdiff, diagfactor);

    // 2. 每簇提链 + 长度过滤
    auto best_chain_clusters = std::make_shared<MatchClusterVec>();
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
inline MatchClusterVec
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



/*!
 * \brief  在原 vector 上直接执行“两级 map + 贪婪拆分”过滤
 * \param  clusters_ptr   ref 上的所有簇
 * \param  MIN_SPAN       最小跨度阈值
 * \param  MIN_MATCH      最小匹配数（若不想改，可写死 50）
 */
void keepWithSplitGreedy(MatchClusterVecPtr clusters_ptr,
    int_t              MIN_SPAN)
{
    if (!clusters_ptr || clusters_ptr->empty()) return;

    /* ---------- 1. 最大堆 ---------- */
    struct Node { MatchCluster cl; int_t sc; };
    auto cmp = [](const Node& a, const Node& b) { return a.sc < b.sc; };

    std::vector<Node> heap;  heap.reserve(clusters_ptr->size());
    for (auto& cl : *clusters_ptr) {
        if (cl.empty()) continue;

        int_t sc = clusterSpan(cl);                // ① 先安全读取
        heap.push_back({ std::move(cl), sc });     // ② 再移动进堆
    }
    clusters_ptr->clear();                              // **清空原容器**
    std::make_heap(heap.begin(), heap.end(), cmp);

    /* ---------- 2. 两级 interval map ---------- */
    IntervalMap refMap;
    std::unordered_map<std::string, IntervalMap> qMaps;

    /* ---------- 3. 贪婪循环 ---------- */
    while (!heap.empty())
    {
        std::pop_heap(heap.begin(), heap.end(), cmp);
        Node cur = std::move(heap.back()); heap.pop_back();
        if (cur.sc < MIN_SPAN) continue;

        uint_t rb = start1(cur.cl.front());
        uint_t re = start1(cur.cl.back()) + len1(cur.cl.back());
        uint_t qb = start2(cur.cl.front());
        uint_t qe = start2(cur.cl.back()) + len2(cur.cl.back());
        const ChrName& qChr = cur.cl.front().query_region.chr_name;

        int_t RL = 0, RR = 0, QL = 0, QR = 0;
        bool ref_hit = overlap1D(refMap, rb, re, RL, RR);
        bool query_hit = overlap1D(qMaps[qChr], qb, qe, QL, QR);

        if (!ref_hit && !query_hit) {
            insertInterval(refMap, rb, re);
            insertInterval(qMaps[qChr], qb, qe);
            clusters_ptr->emplace_back(std::move(cur.cl));   // **写回原 vector**
            continue;
        }


        for (auto& part : splitCluster(cur.cl, ref_hit, RL, RR, query_hit, QL, QR)) {
            int_t sc = clusterSpan(part);
            if (sc >= MIN_SPAN) {
                heap.push_back({ std::move(part), sc });
                std::push_heap(heap.begin(), heap.end(), cmp);
            }

        }
        /* clusters_ptr 现在就是过滤后的结果 */
    }
}

/* ============================================================= *
 *  把三维 clusters  ->  按 ref 的一维 clusters
 *  并行对每个 ref 走 keepWithSplitGreedy，再写回三维结构
 * ============================================================= */
void filterClustersByGreedy(ClusterVecPtrByStrandByQueryRefPtr& cluster_ptr,
    ThreadPool& pool,
    int_t                                min_span)
{
    if (!cluster_ptr || cluster_ptr->empty()) return;

    /* 0. 先按 ref 汇总（保持 shared_ptr 共享语义） */
    ClusterVecPtrByRefPtr by_ref = groupClustersByRef(cluster_ptr);
    if (!by_ref || by_ref->empty()) return;

    for (size_t r = 0; r < by_ref->size(); ++r) {
        pool.enqueue([&, r] {
            keepWithSplitGreedy((*by_ref)[r], min_span);
            });
    }
	pool.waitAllTasksDone();

    return;

}

