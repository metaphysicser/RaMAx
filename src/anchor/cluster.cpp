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

inline uint_t start1(const Match& m) { return static_cast<uint_t>(m.ref_region.start); }
inline uint_t start2(const Match& m) { return static_cast<uint_t>(m.query_region.start); }
inline uint_t len2(const Match& m) { return static_cast<uint_t>(m.query_region.length); }
inline int_t diag(const Match& m) {
    return start2(m) - start1(m);
}

void clusterChrMatch(MatchVec& unique_match, MatchVec& repeat_match, int_t max_gap, int_t diagdiff, double diagfactor)
{
	uint_t match_count = unique_match.size();
    if (match_count < 2) return;

    UnionFind uf(match_count);
    
    for (uint_t i = 0; i < match_count; ++i)                             // i: 0 … N-1
    {
        uint_t i_end = start2(unique_match[i]) + len2(unique_match[i]);
        int_t i_diag = diag(unique_match[i]);


        for (uint_t j = i + 1; j < match_count; ++j)                     // j: i+1 … N-1
        {
            int_t sep = start2(unique_match[j]) - i_end;
            if (sep > static_cast<int_t>(max_gap)) break;               // 太远直接退出内层

            int_t diag_diff = std::abs(diag(unique_match[j]) - i_diag);
            int_t threshold = std::max(diagdiff,
                static_cast<int_t>(diagfactor * sep));

            if (diag_diff <= threshold)
                uf.unite(i, j);                             // 同簇
        }
    }

    /*─────────────────────────────────────────────────────────────*
    * 2. 按根节点把所有元素聚到 cluster                          *
    *─────────────────────────────────────────────────────────────*/
    std::vector<int_t> root_map(match_count, -1);          // 根 → cluster_id
    MatchClusterVec clusters;
    for (uint_t idx = 0; idx < match_count; ++idx)
    {
        int_t root = uf.find(idx);
        int_t cid = root_map[root];
        if (cid == -1) {                       // 新簇
            cid = static_cast<int_t>(clusters.size());
            root_map[root] = cid;
            clusters.emplace_back();           // 创建一个 MatchCluster
        }
        // 把当前 Match 放到 clusters[cid] 中的第一个 MatchVec
        clusters[cid].push_back(std::move(unique_match[idx]));
    }

	return;
}

