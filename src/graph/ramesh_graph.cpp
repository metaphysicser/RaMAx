#include "ramesh.h"

namespace RaMesh {


    GenomeEnd::GenomeEnd()
        : head(HeadSegment::create()),
        tail(TailSegment::create())
    {
        // nothing else
    }

    /* ② 带参：一次性设置 */
    GenomeEnd::GenomeEnd(const std::shared_ptr<HeadSegment>& h,
        const std::shared_ptr<TailSegment>& t)
        : head(h), tail(t) {
    }


    /* ② 只有物种名 */
    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName& sp)
        : species_name(sp)
    {
        /* nothing else; chr2end 空，按需 later add */
    }

    /* ③ 物种名 + 染色体列表 */
    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName& sp,
        const std::vector<ChrName>& chrs)
        : species_name(sp)
    {
        chr2end.reserve(chrs.size());
        for (const auto& c : chrs) {
            chr2end.emplace(c, GenomeEnd{});   // 先放空哨兵，后续可补 head/tail
        }
    }

    RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(
        std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map)
    {
        for (auto& [species_name, mgr_var] : seqpro_map) {

            std::vector<ChrName> chr_names = std::visit(
                [](auto& up)->std::vector<ChrName> {
                    return up ? up->getSequenceNames() : std::vector<ChrName>{};
                }, mgr_var);

            species_graphs.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(species_name),
                std::forward_as_tuple(species_name, chr_names));
        }
    }


} // namespace RaMesh