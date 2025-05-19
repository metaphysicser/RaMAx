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

// ------------------------------------------------------------------
// 高速保存：先写外层数量，再写每个 list 数量，随后连续写 Anchor
// ------------------------------------------------------------------
bool saveAnchors(const std::string& filename,
    const AnchorPtrListVec& anchors)
{
    std::ofstream os(filename, std::ios::binary);
    if (!os) return false;

    cereal::BinaryOutputArchive oar(os);

    // 1) 写外层 vector 大小
    uint64_t outer = anchors.size();
    oar(outer);

    // 2) 逐 list 写入其大小和全部 Anchor
    for (const auto& lst : anchors) {
        uint64_t inner = lst.size();
        oar(inner);
        for (const auto& ap : lst) {
            oar(*ap);            // 只序列化 Anchor 对象本身
        }
    }
    return static_cast<bool>(os);
}

// ------------------------------------------------------------------
// 高速加载：按同样格式读回，现场 new Anchor 并建 shared_ptr
// ------------------------------------------------------------------
bool loadAnchors(const std::string& filename,
    AnchorPtrListVec& anchors)
{
    std::ifstream is(filename, std::ios::binary);
    if (!is) return false;

    cereal::BinaryInputArchive iar(is);

    uint64_t outer;
    iar(outer);

    anchors.clear();
    anchors.reserve(static_cast<size_t>(outer));

    for (uint64_t i = 0; i < outer; ++i) {
        uint64_t inner;
        iar(inner);

        AnchorPtrList lst;
        for (uint64_t j = 0; j < inner; ++j) {
            auto ap = std::make_shared<Anchor>();
            iar(*ap);             // 直接把内容读到 *ap
            lst.push_back(std::move(ap));
        }
        anchors.emplace_back(std::move(lst));
    }
    return static_cast<bool>(is);
}

