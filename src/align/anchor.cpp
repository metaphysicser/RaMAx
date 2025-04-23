#include "anchor.h"

PairGenomeAnchor::PairGenomeAnchor(SpeciesName ref, SpeciesName query,
    AnchorPtrVec init_anchors)
    : ref_species(ref)
    , query_species(query)
    , anchors(std::move(init_anchors))
{
    rebuild();
}

void PairGenomeAnchor::rebuild() {
    std::vector<AnchorEntry> entries;
    entries.reserve(anchors.size());
    for (auto const& ap : anchors) {
        entries.emplace_back(make_box(*ap), ap);
    }
    anchor_rtree = AnchorRTree(entries.begin(), entries.end());
}

void PairGenomeAnchor::insert(AnchorPtr ap) {
    anchors.push_back(ap);
    anchor_rtree.insert({ make_box(*ap), ap });
}

void PairGenomeAnchor::remove(AnchorPtr ap) {
    anchor_rtree.remove({ make_box(*ap), ap });
    anchors.erase(std::remove(anchors.begin(), anchors.end(), ap),
        anchors.end());
}

std::vector<AnchorPtr> PairGenomeAnchor::query(const AnchorBox& region) const {
    std::vector<AnchorEntry> hits;
    anchor_rtree.query(bgi::intersects(region),
        std::back_inserter(hits));
    std::vector<AnchorPtr> out;
    out.reserve(hits.size());
    for (auto const& e : hits) out.push_back(e.second);
    return out;
}