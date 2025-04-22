#ifndef ANCHOR_H
#define ANCHOR_H

#include "config.hpp"
#include <memory>
#include <vector>
#include <algorithm>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// ------------------------------------------------------------------
// 基础别名
// ------------------------------------------------------------------
using Coord_t = uint_t;
using Score_t = int_t;
using CigarUnit = uint32_t;              // 一个 CIGAR 操作
using Cigar_t = std::vector<CigarUnit>; // 整条 CIGAR

// ------------------------------------------------------------------
// 区 域 & 匹 配
// ------------------------------------------------------------------
struct Region {
    Chr chr_name{};
    Coord_t start{ 0 };
    Coord_t length{ 0 };
    Region() = default;
    Region(Chr chr_name, Coord_t s, Coord_t len) : chr_name{chr_name}, start{ s }, length{ len } {}
};
using RegionVec = std::vector<Region>;

enum Strand { FORWARD, REVERSE };

struct Match {
    Region ref_region;
    Region query_region;
    Strand strand;

    Match(Chr r_chr, Coord_t r_start, Coord_t r_len,
        Chr q_chr, Coord_t q_start, Coord_t q_len,
        Strand sd = FORWARD)
        :  ref_region(r_chr, r_start, r_len),
        query_region(q_chr, q_start, q_len),
        strand(sd) {
    }

	Match(Region ref_region, Region query_region, Strand sd = FORWARD) :
		ref_region(ref_region),
		query_region(query_region),
		strand(sd) {
	}
    Match() = default;
};

struct Anchor {
    Match   match;
    Coord_t alignment_length{ 0 };
    Cigar_t cigar;
    Score_t alignment_score{ 0 };
    Anchor() = default;
    Anchor(const Match& m, Coord_t aln_len,
        const Cigar_t& c, Score_t score)
        : match(m), alignment_length(aln_len), cigar(c), alignment_score(score) {
    }
};
using AnchorVec = std::vector<Anchor>;
// Alias for shared pointers to Anchor
using AnchorPtr = std::shared_ptr<Anchor>;
using AnchorPtrVec = std::vector<AnchorPtr>;

// 空间索引类型
using AnchorPoint = bg::model::point<Coord_t, 2, bg::cs::cartesian>;
using AnchorBox = bg::model::box<AnchorPoint>;
using AnchorEntry = std::pair<AnchorBox, AnchorPtr>;
using AnchorRTree = bgi::rtree<AnchorEntry, bgi::quadratic<64>>;

inline AnchorBox make_box(const Anchor& a) {
    const auto& r = a.match.ref_region;
    const auto& q = a.match.query_region;
    AnchorPoint min_pt{ r.start, q.start };
    AnchorPoint max_pt{ r.start + r.length, q.start + q.length };
    return AnchorBox{ min_pt, max_pt };
}

// ------------------------------------------------------------------
// 类：PairGenomeAnchor
// ------------------------------------------------------------------
class PairGenomeAnchor {
public:
    PairGenomeAnchor() = default;
    PairGenomeAnchor(Species ref, Species query,
        AnchorPtrVec init_anchors);

    void rebuild();
    void insert(AnchorPtr ap);
    void remove(AnchorPtr ap);
    std::vector<AnchorPtr> query(const AnchorBox& region) const;

private:
    Species       ref_species{};
    Species       query_species{};
    AnchorPtrVec  anchors{};
    AnchorRTree   anchor_rtree{};
};
#endif // ANCHOR_H
