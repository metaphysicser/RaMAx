#include "anchor.h"

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
    }
    return clusters;
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
