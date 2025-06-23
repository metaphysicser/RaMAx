#ifndef ANCHOR_H
#define ANCHOR_H

// ------------------------------------------------------------------
// 引入必要的头文件
// ------------------------------------------------------------------
#include "config.hpp"              // 包含通用配置（如类型定义）
#include "align.h"
#include <memory>                 // 用于智能指针（如 std::shared_ptr）
#include <vector>                 // 用于 std::vector 容器
#include <algorithm>              // 提供算法（如 sort 等）
#include <list>                   // 用于 std::list 容器
#include <SeqPro.h>
#include <queue>


#define ANCHOR_EXTENSION "anchor"  // Anchor 文件保存使用的默认扩展名

// ------------------------------------------------------------------
// 类型别名定义
// ------------------------------------------------------------------
// 定义坐标类型（可在 config.hpp 中设为 uint32_t 或 uint64_t）
using Coord_t = uint_t;

/* 区间与 map 类型 */
using IntervalMap = std::map<int_t, int_t>;          // key = beg, value = end   (右开)

/* 插入非重叠区间到 map（调用处已保证不重叠，只需 emplace 即可） */
inline void insertInterval(IntervalMap& m, int_t l, int_t r) { m.emplace(l, r); }

/* 判断 [l,r) 是否与 map 中已有区间重叠；若重叠，把那条区间的 [L,R) 返到 outBox */
inline bool overlap1D(const IntervalMap& m,
    int_t l, int_t r,
    int_t& L, int_t& R)
{
    if (m.empty()) return false;

    auto it = m.lower_bound(l);          // 第一个 beg ≥ l
    if (it != m.end() && it->first < r) {  // 与右邻区间重叠
        L = it->first; R = it->second;
        return true;
    }
    if (it != m.begin()) {                 // 与左邻区间重叠？
        --it;
        if (it->second > l) {
            L = it->first; R = it->second;
            return true;
        }
    }
    return false;
}


// ------------------------------------------------------------------
// 区域 Region 与 比对 Match 的结构定义
// ------------------------------------------------------------------

// 表示基因组上的一个区域
struct Region {
    ChrName chr_name{};     // 染色体名称
    Coord_t start{ 0 };     // 起始坐标
    Coord_t length{ 0 };    // 区域长度

    Region() = default;
    Region(ChrName chr_name, Coord_t s, Coord_t len)
        : chr_name{ chr_name }, start{ s }, length{ len } {
    }
};

using RegionVec = std::vector<Region>;  // 区域集合

// 链接方向：正向或反向
enum Strand { FORWARD, REVERSE };

// 表示参考与查询序列之间的一段匹配区域
struct Match {
    Region ref_region;      // 参考基因组上的区域
    Region query_region;    // 查询基因组上的区域
    Strand strand;          // 链接方向

    // 构造函数：使用详细参数创建一个 Match
    Match(ChrName r_chr, Coord_t r_start, Coord_t r_len,
        ChrName q_chr, Coord_t q_start, Coord_t q_len,
        Strand sd = FORWARD)
        : ref_region(r_chr, r_start, r_len),
        query_region(q_chr, q_start, q_len),
        strand(sd) {
    }

    // 构造函数：使用已构建的 Region
    Match(Region ref_region, Region query_region, Strand sd = FORWARD)
        : ref_region(ref_region),
        query_region(query_region),
        strand(sd) {
    }

    Match() = default;
};

using MatchVec = std::vector<Match>;
using MatchPtr = std::shared_ptr<Match>; // 智能指针类型
using MatchPtrVec = std::vector<MatchPtr>; // 智能指针集合

using MatchVec2D = std::vector<MatchVec>;
using MatchVec2DPtr = std::shared_ptr<MatchVec2D>;
using MatchVec3D = std::vector<std::vector<MatchVec>>;
using MatchVec3DPtr = std::shared_ptr<MatchVec3D>;

using SpeciesMatchVec3DPtrMap = std::unordered_map<std::string, MatchVec3DPtr>;
using SpeciesMatchVec3DPtrMapPtr = std::shared_ptr<SpeciesMatchVec3DPtrMap>;

using MatchByRef = std::vector<MatchVec>;
using MatchByQueryRef = std::vector<MatchByRef>;
using MatchByStrandByQueryRef = std::vector<MatchByQueryRef>;
using MatchByStrandByQueryRefPtr = std::shared_ptr<MatchByStrandByQueryRef>;

using SpeciesMatchByStrandByQueryRefPtrMap = std::unordered_map<SpeciesName, MatchByStrandByQueryRefPtr>;

using MatchCluster = MatchVec; // 匹配簇，包含多个匹配向量
using MatchClusterVec = std::vector<MatchCluster>;
using MatchClusterVecPtr = std::shared_ptr<MatchClusterVec>;

using ClusterSynteny = std::vector<MatchCluster>;
using ClusterSyntenyVec = std::vector<ClusterSynteny>;
using ClusterSyntenyVecPtr = std::shared_ptr<ClusterSyntenyVec>;
using ClusterSyntenyVecPtrByRef = std::vector<ClusterSyntenyVecPtr>;

using ClusterVecPtrByRef = std::vector<MatchClusterVecPtr>;
using ClusterVecPtrByQueryRef = std::vector<ClusterVecPtrByRef>;
using ClusterVecPtrByStrandByQueryRef = std::vector<ClusterVecPtrByQueryRef>;

using ClusterVecPtrByRefPtr = std::shared_ptr<ClusterVecPtrByRef>;
using ClusterVecPtrByRefQuery = std::vector<ClusterVecPtrByRef>;
using ClusterVecPtrByRefQueryPtr = std::shared_ptr<ClusterVecPtrByRefQuery>;
using ClusterVecPtrByStrandByQueryRefPtr = std::shared_ptr<ClusterVecPtrByStrandByQueryRef>;

using SpeciesClusterMap =
std::unordered_map<SpeciesName, ClusterVecPtrByStrandByQueryRefPtr>;
using SpeciesClusterMapPtr = std::shared_ptr<SpeciesClusterMap>;

inline uint_t start1(const Match& m) { return static_cast<uint_t>(m.ref_region.start); }
inline uint_t start2(const Match& m) { return static_cast<uint_t>(m.query_region.start); }
inline uint_t len1(const Match& m) { return static_cast<uint_t>(m.ref_region.length); }
inline uint_t len2(const Match& m) { return static_cast<uint_t>(m.query_region.length); }
inline int_t diag(const Match& m) {
    return start2(m) - start1(m);
}

// 3. 彻底释放 cluster 内剩余元素的内存（swap 技巧）
inline void releaseCluster(MatchVec& cluster) {
    MatchVec tmp; tmp.swap(cluster);
}

/* ---------- 辅助：计算簇跨度 ---------- */
inline int_t clusterSpan(const MatchCluster& cl)
{
    return start1(cl.back()) + len1(cl.back())
        - start1(cl.front());
}

void groupMatchByQueryRef(MatchVec3DPtr& anchors,
    MatchByStrandByQueryRefPtr unique_anchors,
    MatchByStrandByQueryRefPtr repeat_anchors,
    SeqPro::ManagerVariant& ref_fasta_manager,
    SeqPro::ManagerVariant& query_fasta_manager,
    ThreadPool& pool);

MatchClusterVecPtr
groupClustersToVec(const ClusterVecPtrByStrandByQueryRefPtr& src,
    ThreadPool& pool, uint_t thread_num);


ClusterVecPtrByRefPtr
groupClustersByRef(const ClusterVecPtrByStrandByQueryRefPtr& src);
// void sortMatchByRefStart(MatchByQueryRefPtr& anchors, ThreadPool& pool);

ClusterVecPtrByRefQueryPtr
groupClustersByRefQuery(const ClusterVecPtrByRefPtr& by_ref,
    SeqPro::ManagerVariant& query_fasta_manager,
    ThreadPool& pool);

void sortMatchByQueryStart(MatchByStrandByQueryRefPtr& anchors, ThreadPool& pool);

MatchClusterVecPtr clusterChrMatch(MatchVec& unique_match, MatchVec& repeat_match, int_t max_gap = 90, int_t diagdiff = 5, double diagfactor = 0.12, int_t min_cluster_length = 50);

MatchVec bestChainDP(MatchVec& cluster, double diagfactor);

MatchClusterVec buildClusters(MatchVec& unique_match,
    int_t      max_gap,
    int_t      diagdiff,
    double     diagfactor);

ClusterVecPtrByStrandByQueryRefPtr
clusterAllChrMatch(const MatchByStrandByQueryRefPtr& unique_anchors,
    const MatchByStrandByQueryRefPtr& repeat_anchors,
    ThreadPool& pool);

MatchClusterVec
splitCluster(const MatchCluster& cl,
    bool  ref_hit,
    int_t bad_r_beg, int_t bad_r_end,
    bool  query_hit,
    int_t bad_q_beg, int_t bad_q_end);

// 一个比对锚点（Anchor）表示一对匹配区域之间的精确比对信息
struct Anchor {
    Match match;                    // 匹配信息
    Coord_t alignment_length{ 0 }; // 对齐长度
    Cigar_t cigar;                 // 对齐的 CIGAR 字符串
    Score_t alignment_score{ 0 };  // 对齐得分

    Anchor() = default;

    Anchor(const Match m, Coord_t aln_len,
        const Cigar_t c, Score_t score)
        : match(m), alignment_length(aln_len),
        cigar(c), alignment_score(score) {
    }
};

using AnchorVec = std::vector<Anchor>;           // 一组 Anchor
using AnchorPtr = std::shared_ptr<Anchor>;       // Anchor 的智能指针
using AnchorPtrVec = std::vector<AnchorPtr>;    // Anchor 智能指针的集合
using AnchorPtrList = std::list<AnchorPtr>;      // Anchor 指针链表
using AnchorPtrListVec = std::vector<AnchorPtrList>; // Anchor 指针链表的集合（如按染色体分组）
using AnchorVec2D = std::vector<AnchorVec>;
using AnchorVec2DPtr = std::shared_ptr<AnchorVec2D>;
using AnchorVec3D = std::vector<std::vector<AnchorVec>>;
using AnchorVec3DPtr = std::shared_ptr<AnchorVec3D>;

using AnchorsByRef = std::vector<AnchorVec>;
using AnchorsByQueryRef = std::vector<AnchorsByRef>;


// ------------------------------------------------------------------
// 智能序列分块函数
// ------------------------------------------------------------------

/**
 * @brief 智能预分配序列分块，支持自动分块策略选择
 * 
 * 该函数会根据数据量大小智能选择分块策略：
 * - 当序列数量较多时（>= max_sequences_threshold），按序列划分而不是按长度切分
 * - 当序列数量较少但总长度很大时，按指定大小切分序列并保留重叠区域
 * 
 * @param seq_manager SeqPro::Manager对象，用于获取序列信息
 * @param chunk_size 每个chunk的目标大小（碱基数）
 * @param overlap_size chunk之间的重叠大小（碱基数）
 * @param max_sequences_threshold 序列数量阈值，超过此值时按序列划分
 * @param min_chunk_size 最小chunk大小，小于此值的序列不会被进一步切分
 * @return RegionVec 分块结果，每个Region代表一个chunk
 */
RegionVec preAllocateChunks(const SeqPro::ManagerVariant& seq_manager,
                           uint_t chunk_size,
                           uint_t overlap_size = 0,
                           size_t max_sequences_threshold = 1000,
                           uint_t min_chunk_size = 10000);

/**
 * @brief 按序列划分的简单分块策略
 * 
 * 当序列数量较多时使用，每个序列作为一个独立的chunk
 * 
 * @param seq_manager SeqPro::Manager对象
 * @return RegionVec 分块结果，每个序列一个chunk
 */
RegionVec preAllocateChunksBySequence(const SeqPro::ManagerVariant& seq_manager);

/**
 * @brief 按大小划分的分块策略（原始FastaManager逻辑）
 * 
 * 当序列数量较少但总长度很大时使用，会对长序列进行切分
 * 
 * @param seq_manager SeqPro::Manager对象
 * @param chunk_size 每个chunk的目标大小
 * @param overlap_size chunk之间的重叠大小
 * @param min_chunk_size 最小chunk大小，小于此值的序列不会被切分
 * @return RegionVec 分块结果
 */
RegionVec preAllocateChunksBySize(const SeqPro::ManagerVariant& seq_manager,
                                 uint_t chunk_size,
                                 uint_t overlap_size = 0,
                                 uint_t min_chunk_size = 10000);

// ------------------------------------------------------------------
// Cereal 序列化支持：用于将结构写入文件或从文件读取
// ------------------------------------------------------------------
namespace cereal {

    // Region 的序列化
    template <class Archive>
    void serialize(Archive& ar, Region& r) {
        ar(r.chr_name, r.start, r.length);
    }

    // Match 的序列化
    template <class Archive>
    void serialize(Archive& ar, Match& m) {
        ar(m.ref_region, m.query_region, m.strand);
    }

    // Anchor 的序列化
    template <class Archive>
    void serialize(Archive& ar, Anchor& a) {
        ar(a.match,
            a.alignment_length,
            a.cigar,
            a.alignment_score);
    }

} // namespace cereal

// ------------------------------------------------------------------
// Anchor 文件的读写接口
// ------------------------------------------------------------------

// 保存 anchors 到文件中（返回 true 表示成功）
bool saveAnchors(const std::string& filename, const AnchorPtrListVec& anchors);

// 从文件中读取 anchors（返回 true 表示成功）
bool loadAnchors(const std::string& filename, AnchorPtrListVec& anchors);

// 保存：all_sets -> filename
bool saveAnchorsSets(const std::string& filename,
    const std::vector<AnchorPtrListVec>& all_sets);

// 读取：filename -> all_sets
bool loadAnchorsSets(const std::string& filename,
    std::vector<AnchorPtrListVec>& all_sets);

bool loadAnchorVec3D(const std::string& filename, AnchorVec3DPtr& data);
bool saveAnchorVec3D(const std::string& filename, const AnchorVec3DPtr& data);

bool saveMatchVec3D(const std::string& filename, const MatchVec3DPtr& data);
bool loadMatchVec3D(const std::string& filename, MatchVec3DPtr& data);

bool loadSpeciesMatchMap(const std::string& filename, SpeciesMatchVec3DPtrMapPtr& map_ptr);
bool saveSpeciesMatchMap(const std::string& filename, const SpeciesMatchVec3DPtrMapPtr& map_ptr);

class UnionFind
{
public:
    /// 构造并初始化为 n 个独立元素
    explicit UnionFind(std::size_t n = 0);

    /// 清空并重建大小为 n 的并查集
    void reset(std::size_t n);

    /// 查找 x 所在集合的代表（根）
    int_t find(int_t x);

    /// 返回 x 所在集合的元素个数
    int_t set_size(int_t x);

    /// 判断 a 与 b 是否位于同一集合
    bool same(int_t a, int_t b);

    /// 当前独立集合个数
    int_t components() const { return component_cnt_; }

    /**
     * @brief 合并 a、b 所在集合
     * @return true  如果二者原本不连通并成功合并
     *         false 如果二者已在同一集合
     */
    bool unite(int_t a, int_t b);

private:
    std::vector<int_t> parent_;   // 负数：根 & -size；非负：父索引
    int_t component_cnt_{ 0 };      // 当前集合数
};

#endif // ANCHOR_H
