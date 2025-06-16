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

ClusterVecPtrByRefPtrVec
groupClustersByRef(ClusterVecPtrByStrandByQueryRefPtr& src)
{
    ClusterVecPtrByRefPtrVec by_ref;
    if (!src || src->empty()) return by_ref;

    const uint_t strand_n = static_cast<uint_t>(src->size());

    /* ---------- ① 预扫描：确定 ref_n、query_max ---------- */
    uint_t ref_n = 0, query_max = 0;
    for (const auto& q_layer_by_strand : *src) {
        query_max = std::max<uint_t>(query_max,
            static_cast<uint_t>(q_layer_by_strand.size()));
        for (const auto& r_row_by_query : q_layer_by_strand)
            ref_n = std::max<uint_t>(ref_n,
                static_cast<uint_t>(r_row_by_query.size()));
    }
    if (ref_n == 0) return by_ref;

    /* ---------- ② 初始化 by_ref 并提前 reserve ---------- */
    by_ref.resize(ref_n);
    const size_t expect_per_ref = static_cast<size_t>(strand_n) * query_max; // 上界
    for (auto& vec_ptr : by_ref) {
        vec_ptr = std::make_shared<ClusterVecPtrByRef>();
        vec_ptr->reserve(expect_per_ref);     // 可显著减少 push_back 的 realloc 次数
    }

    /* ---------- ③ 单趟扫描直接归位 ---------- */
    for (uint_t s = 0; s < strand_n; ++s) {
        const auto& q_layer = (*src)[s];
        for (const auto& r_row : q_layer) {
            for (uint_t r = 0; r < r_row.size(); ++r)
                by_ref[r]->push_back(r_row[r]);   // 仅拷贝 shared_ptr（轻量）
        }
    }

    return by_ref;
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
MatchClusterVecPtr clusterChrMatch(MatchVec& unique_match,
    MatchVec& repeat_match,
    int_t      max_gap,
    int_t      diagdiff,
    double     diagfactor,
    int_t      min_cluster_length)
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

/* ---------- 辅助：簇按 query/ref 连续性拆分 ---------- */
std::vector<MatchCluster> splitCluster(const MatchCluster& cl,
    int_t bad_q_beg, int_t bad_q_end,
    int_t bad_r_beg, int_t bad_r_end)
{
    std::vector<MatchCluster> parts;
    MatchCluster current;

    for (const auto& m : cl)
    {
        int_t qb = start1(m), qe = qb + len1(m);
        int_t rb = start2(m), re = rb + len2(m);

        bool hit = !(qe <= bad_q_beg || qb >= bad_q_end ||
            re <= bad_r_beg || rb >= bad_r_end);
        if (hit) {
            if (!current.empty()) parts.emplace_back(std::move(current));
            current.clear();
        }
        else {
            current.emplace_back(m);
        }
    }
    if (!current.empty()) parts.emplace_back(std::move(current));
    return parts;
}

/* ===================================================================== *
 * 主函数：无需 R-tree，仅用两级 map
 * ===================================================================== */
std::shared_ptr<MatchClusterVec>
keepWithSplitGreedy(std::shared_ptr<MatchClusterVec> clusters_ptr,
    int_t MIN_SPAN,
    int_t MIN_MATCH)
{
    if (!clusters_ptr || clusters_ptr->empty()) return clusters_ptr;

    /* ---------- 1. 最大堆（一次 make_heap） ---------- */
    struct Node { MatchCluster cl; int_t sc; };
    auto cmp = [](const Node& a, const Node& b) { return a.sc < b.sc; };

    std::vector<Node> heap;  heap.reserve(clusters_ptr->size());
    for (auto& cl : *clusters_ptr)
        if (!cl.empty())
            heap.push_back({ std::move(cl), clusterSpan(cl) });
    clusters_ptr->clear();
    std::make_heap(heap.begin(), heap.end(), cmp);

    /* ---------- 2. 两级 interval map ---------- */
    IntervalMap refMap;                                         // 已选 ref 区间
    std::unordered_map<std::string, IntervalMap> qMaps;         // 每条 query-chr

    auto keep_ptr = std::make_shared<MatchClusterVec>();
    keep_ptr->reserve(heap.size());

    /* ---------- 3. 贪心循环 ---------- */
    while (!heap.empty())
    {
        std::pop_heap(heap.begin(), heap.end(), cmp);
        Node cur = std::move(heap.back()); heap.pop_back();
        if (cur.sc < MIN_SPAN || cur.cl.size() < (size_t)MIN_MATCH) continue;

        /* 当前簇的端点 & 染色体名 */
        int_t qb = start1(cur.cl.front());
        int_t qe = start1(cur.cl.back()) + len1(cur.cl.back());
        int_t rb = start2(cur.cl.front());
        int_t re = start2(cur.cl.back()) + len2(cur.cl.back());
        const ChrName& qChr = cur.cl.front().query_region.chr_name;

        /* ---- 3-1. 判重 ---- */
        int_t RL = 0, RR = 0, QL = 0, QR = 0;
        bool refHit = overlap1D(refMap, rb, re, RL, RR);
        bool queryHit = overlap1D(qMaps[qChr], qb, qe, QL, QR);

        if (!refHit && !queryHit) {
            /* ---- 3-2. 无冲突：选中并更新两张表 ---- */
            insertInterval(refMap, rb, re);
            insertInterval(qMaps[qChr], qb, qe);
            keep_ptr->emplace_back(std::move(cur.cl));
            continue;
        }

        /* ---- 3-3. 有冲突：按触发的维度拆分 ---- */
        if (refHit) {                           // 先按 ref 轴拆
            for (auto& part : splitCluster(cur.cl, RL, RR, qb, qe)) {
                int_t sc = clusterSpan(part);
                if (sc >= MIN_SPAN && part.size() >= (size_t)MIN_MATCH) {
                    heap.push_back({ std::move(part), sc });
                    std::push_heap(heap.begin(), heap.end(), cmp);
                }
            }
        }
        else {                                // 只在 query 轴冲突
            for (auto& part : splitCluster(cur.cl, qb, qe, QL, QR)) {
                int_t sc = clusterSpan(part);
                if (sc >= MIN_SPAN && part.size() >= (size_t)MIN_MATCH) {
                    heap.push_back({ std::move(part), sc });
                    std::push_heap(heap.begin(), heap.end(), cmp);
                }
            }
        }
    }
    return keep_ptr;
}


