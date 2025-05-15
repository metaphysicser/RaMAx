#include "anchor.h"

PairGenomeAnchor::PairGenomeAnchor(SpeciesName ref,
    SpeciesName query,
    AnchorPtrListVec init_anchors)
    : ref_species(ref)
    , query_species(query)
    , anchors(std::move(init_anchors))
{
    rebuild();
}

void PairGenomeAnchor::rebuild() {
    // 1) 清空旧的树
    anchor_rtree.clear();

    // 2) 扁平化地把所有 list 中的 AnchorPtr 插入 R-tree
    for (auto const& list_group : anchors) {
        for (auto const& ap : list_group) {
            anchor_rtree.insert({ make_box(*ap), ap });
        }
    }
}

void PairGenomeAnchor::insertAnchorRtree(AnchorPtr ap) {
    anchor_rtree.insert({ make_box(*ap), ap });
}

void PairGenomeAnchor::removeAnchorRtree(AnchorPtr ap) {
    // 1) 从 R-tree 中移除
    anchor_rtree.remove({ make_box(*ap), ap });

}

std::vector<AnchorPtr> PairGenomeAnchor::query(const AnchorBox& region) const {
    // 在 R-tree 中查找所有相交的 AnchorEntry
    std::vector<AnchorEntry> hits;
    anchor_rtree.query(bgi::intersects(region),
        std::back_inserter(hits));

    // 抽取出所有的 AnchorPtr
    std::vector<AnchorPtr> out;
    out.reserve(hits.size());
    for (auto const& e : hits) {
        out.push_back(e.second);
    }
    return out;
}

// 返回 true/false 表示保存是否成功
bool saveAnchors(const std::string& filename, const AnchorPtrListVec& anchors) {
    std::ofstream os(filename, std::ios::binary);
    if (!os) return false;
    cereal::BinaryOutputArchive oar(os);
    oar(anchors);
    return static_cast<bool>(os);
}

// 将文件里的结果读到 anchors，返回是否成功
bool loadAnchors(const std::string& filename, AnchorPtrListVec& anchors) {
    std::ifstream is(filename, std::ios::binary);
    if (!is) return false;
    cereal::BinaryInputArchive iar(is);
    iar(anchors);
    return static_cast<bool>(is);
}
