#ifndef ANCHOR_H
#define ANCHOR_H

// ------------------------------------------------------------------
// 引入必要的头文件
// ------------------------------------------------------------------
#include "config.hpp"              // 包含通用配置（如类型定义）
#include "cigar.h"                 // CIGAR 字符串表示比对信息
#include <memory>                 // 用于智能指针（如 std::shared_ptr）
#include <vector>                 // 用于 std::vector 容器
#include <algorithm>              // 提供算法（如 sort 等）
#include <cereal/types/vector.hpp> // cereal 序列化支持 vector
#include <cereal/types/list.hpp>   // cereal 序列化支持 list
#include <cereal/types/memory.hpp> // cereal 序列化支持智能指针

// boost 空间索引相关头文件（当前注释掉）
// #include <boost/geometry.hpp>
// #include <boost/geometry/geometries/point.hpp>
// #include <boost/geometry/geometries/box.hpp>
// #include <boost/geometry/index/rtree.hpp>

#define ANCHOR_EXTENSION "anchor"  // Anchor 文件保存使用的默认扩展名

// 前向声明
class FastaManager;
// ------------------------------------------------------------------
// 类型别名定义
// ------------------------------------------------------------------
// 定义坐标类型（可在 config.hpp 中设为 uint32_t 或 uint64_t）
using Coord_t = uint_t;

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

using MatchByRef = std::vector<MatchVec>;
using MatchByQueryRef = std::vector<MatchByRef>;
using MatchByStrandByQueryRef = std::vector<MatchByQueryRef>;

using MatchCluster = MatchVec; // 匹配簇，包含多个匹配向量
using MatchClusterVec = std::vector<MatchCluster>;
using MatchClusterVecPtr = std::shared_ptr<MatchClusterVec>;
using ClusterVecPtrByRefByQuery = std::vector<std::vector<MatchClusterVecPtr>>;

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

void groupMatchByQueryRef(MatchVec3DPtr& anchors,
    MatchByStrandByQueryRef& unique_anchors,
    MatchByStrandByQueryRef& repeat_anchors,
    FastaManager& ref_fasta_manager,
    FastaManager& query_fasta_manager,
    uint_t thread_num);

void sortMatchByRefStart(MatchByQueryRef& anchors, uint_t thread_num);
void sortMatchByQueryStart(MatchByQueryRef& anchors, uint_t thread_num);

AnchorPtrVec findNonOverlapAnchors(const AnchorVec& anchors);

MatchClusterVecPtr clusterChrMatch(MatchVec& unique_match, MatchVec& repeat_match, int_t max_gap = 90, int_t diagdiff = 5, double diagfactor = 0.12, int_t min_cluster_length = 50);

MatchVec bestChainDP(MatchVec& cluster, double diagfactor);

MatchClusterVec buildClusters(MatchVec& unique_match,
    int_t      max_gap,
    int_t      diagdiff,
    double     diagfactor);

// ------------------------------------------------------------------
// 空间索引（注释掉的部分，若启用 Boost RTree 可用于高效的空间查询）
// ------------------------------------------------------------------
// using AnchorPoint = bg::model::point<Coord_t, 2, bg::cs::cartesian>;
// using AnchorBox = bg::model::box<AnchorPoint>;
// using AnchorEntry = std::pair<AnchorBox, AnchorPtr>;
// using AnchorRTree = bgi::rtree<AnchorEntry, bgi::quadratic<64>>;
//
// // 构建一个 Anchor 的空间包围盒，用于空间索引
// inline AnchorBox make_box(const Anchor& a) {
//     const auto& r = a.match.ref_region;
//     const auto& q = a.match.query_region;
//     AnchorPoint min_pt{ r.start, q.start };
//     AnchorPoint max_pt{ r.start + r.length, q.start + q.length };
//     return AnchorBox{ min_pt, max_pt };
// }

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
// ------------------------------------------------------------------
// （注释掉）双基因组比对结构体：用于管理两个物种之间的 Anchor 信息
// ------------------------------------------------------------------
// class PairGenomeAnchor {
// public:
//     PairGenomeAnchor() = default;
//     PairGenomeAnchor(SpeciesName ref, SpeciesName query,
//         AnchorPtrListVec init_anchors);
//
//     void rebuild();                            // 重新构建 R-Tree
//     void insertAnchorRtree(AnchorPtr ap);      // 插入 anchor 到空间索引
//     void removeAnchorRtree(AnchorPtr ap);      // 从空间索引中移除 anchor
//     std::vector<AnchorPtr> query(const AnchorBox& region) const; // 查询 anchor
//
// private:
//     SpeciesName       ref_species{};           // 参考物种名称
//     SpeciesName       query_species{};         // 查询物种名称
//     AnchorPtrListVec  anchors{};               // anchor 列表
//     AnchorRTree   anchor_rtree{};              // 空间索引树
// };

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
