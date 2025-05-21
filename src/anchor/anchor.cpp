#include "anchor.h"

//PairGenomeAnchor::PairGenomeAnchor(SpeciesName ref,
//    SpeciesName query,
//    AnchorPtrListVec init_anchors)
//    : ref_species(ref)
//    , query_species(query)
//    , anchors(std::move(init_anchors))
//{
//    rebuild();
//}
//
//void PairGenomeAnchor::rebuild() {
//    // 1) ��վɵ���
//    anchor_rtree.clear();
//
//    // 2) ��ƽ���ذ����� list �е� AnchorPtr ���� R-tree
//    for (auto const& list_group : anchors) {
//        for (auto const& ap : list_group) {
//            anchor_rtree.insert({ make_box(*ap), ap });
//        }
//    }
//}
//
//void PairGenomeAnchor::insertAnchorRtree(AnchorPtr ap) {
//    anchor_rtree.insert({ make_box(*ap), ap });
//}
//
//void PairGenomeAnchor::removeAnchorRtree(AnchorPtr ap) {
//    // 1) �� R-tree ���Ƴ�
//    anchor_rtree.remove({ make_box(*ap), ap });
//
//}
//
//std::vector<AnchorPtr> PairGenomeAnchor::query(const AnchorBox& region) const {
//    // �� R-tree �в��������ཻ�� AnchorEntry
//    std::vector<AnchorEntry> hits;
//    anchor_rtree.query(bgi::intersects(region),
//        std::back_inserter(hits));
//
//    // ��ȡ�����е� AnchorPtr
//    std::vector<AnchorPtr> out;
//    out.reserve(hits.size());
//    for (auto const& e : hits) {
//        out.push_back(e.second);
//    }
//    return out;
//}

// ------------------------------------------------------------------
// ���ٱ��棺��д�����������дÿ�� list �������������д Anchor
// ------------------------------------------------------------------
bool saveAnchors(const std::string& filename,
    const AnchorPtrListVec& anchors)
{
    std::ofstream os(filename, std::ios::binary);
    if (!os) return false;

    cereal::BinaryOutputArchive oar(os);

    // 1) д��� vector ��С
    uint64_t outer = anchors.size();
    oar(outer);

    // 2) �� list д�����С��ȫ�� Anchor
    for (const auto& lst : anchors) {
        uint64_t inner = lst.size();
        oar(inner);
        for (const auto& ap : lst) {
            oar(*ap);            // ֻ���л� Anchor ������
        }
    }
    return static_cast<bool>(os);
}

// ------------------------------------------------------------------
// ���ټ��أ���ͬ����ʽ���أ��ֳ� new Anchor ���� shared_ptr
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
            iar(*ap);             // ֱ�Ӱ����ݶ��� *ap
            lst.push_back(std::move(ap));
        }
        anchors.emplace_back(std::move(lst));
    }
    return static_cast<bool>(is);
}

