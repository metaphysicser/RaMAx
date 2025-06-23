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

    void RaMeshMultiGenomeGraph::insertClusterIntoGraph(SpeciesName ref_name, SpeciesName query_name, const MatchCluster& cluster)
    {
		if (cluster.empty()) return; // 空簇不处理   

		/* —— 1. 获取 ref_region 的物种名 —— */
        SpeciesName ref_chr = cluster.front().ref_region.chr_name;
		SpeciesName query_chr = cluster.front().query_region.chr_name; 

		std::shared_ptr<HeadSegment> ref_head = species_graphs[ref_name].chr2end[ref_chr].head;
		std::shared_ptr<HeadSegment> query_head = species_graphs[query_name].chr2end[query_chr].head;

        for (const auto& m : cluster) {
            auto blk = RaMesh::Block::create_from_match(m);  
        }
    }


} // namespace RaMesh