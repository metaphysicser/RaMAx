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

    /* ① 默认：占位构造 */
    RaMeshGenomeGraph::RaMeshGenomeGraph() = default;

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

} // namespace RaMesh