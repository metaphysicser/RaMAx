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

void groupAnchorsByQueryRef(AnchorVec3DPtr& anchors,
    AnchorsByQueryRef& unique_anchors,
    AnchorsByQueryRef& repeat_anchors,
    FastaManager& ref_fasta_manager,
    FastaManager& query_fasta_manager,
    uint_t thread_num);

void sortAnchorsByRefStart(AnchorsByQueryRef& anchors, uint_t thread_num);

AnchorPtrVec findNonOverlapAnchors(const AnchorVec& anchors);

void clusterChrAnchors(const AnchorVec& anchors);
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

#endif // ANCHOR_H
