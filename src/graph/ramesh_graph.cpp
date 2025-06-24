/*******************************************************
*  ramesh.cpp – 适配 ramesh.h (v0.4-rev, 单一 Segment)
*  2025-06-24
*******************************************************/
#include "ramesh.h"
#include <iomanip>
#include <shared_mutex>

namespace RaMesh {

    /* ------------------------------------------------------------------ */
    /* 1. RaMeshGenomeGraph                                               */
    /* ------------------------------------------------------------------ */
    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName& sp)
        : species_name(sp) {
    }

    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName& sp,
        const std::vector<ChrName>& chrs)
        : species_name(sp)
    {
        chr2end.reserve(chrs.size());
        for (const auto& c : chrs) chr2end.emplace(c, GenomeEnd{});
    }

    void RaMeshGenomeGraph::debug_print(std::ostream& os,
        bool show_secondary) const
    {
        std::shared_lock lg(rw);                 // 读锁

        os << "=== GenomeGraph <" << species_name << "> ===\n";
        for (const auto& [chr, end] : chr2end)
        {
            os << "\n[Chromosome " << chr << "]\n";
            os << std::left << std::setw(6) << "Idx"
                << std::setw(12) << "Start"
                << std::setw(10) << "Len"
                << std::setw(4) << "Str"
                << std::setw(6) << "Role";
            if (show_secondary) os << "  |  Mate(ptr)";
            os << '\n';

            Segment* headSeg = end.head->load();
            Segment* tailSeg = end.tail->load();

            SegAtomPtr currAtom = headSeg->primary_path.next;
            std::size_t idx = 0;

            while (currAtom && currAtom->load() != tailSeg)
            {
                const Segment* seg = currAtom->load();

                os << std::setw(6) << idx++
                    << std::setw(12) << seg->start
                    << std::setw(10) << seg->length
                    << std::setw(4) << (seg->strand == Strand::FORWARD ? "+" : "-")
                    << std::setw(6) << (seg->is_primary() ? "Pri" : "Sec");

                if (show_secondary) {
                    Segment* mate = seg->secondary_path.next
                        ? seg->secondary_path.next->load()
                        : nullptr;
                    os << "  |  " << mate;
                }
                os << '\n';

                currAtom = seg->primary_path.next;
            }
            os << "Total segments: " << idx << '\n';
        }
        os << "=== End of Graph ===\n";
    }

    /* ------------------------------------------------------------------ */
    /* 2. RaMeshMultiGenomeGraph – ctor                                   */
    /* ------------------------------------------------------------------ */
    //RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(
    //    std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map)
    //{
    //    for (auto& [sp, mgr_var] : seqpro_map) {
    //        std::vector<ChrName> chr_names = std::visit(
    //            [](auto& up)->std::vector<ChrName> {
    //                return up ? up->getSequenceNames()
    //                    : std::vector<ChrName>{};
    //            }, mgr_var);

    //        species_graphs.emplace(
    //            std::piecewise_construct,
    //            std::forward_as_tuple(sp),
    //            std::forward_as_tuple(sp, chr_names));
    //    }
    //}

    /* ------------------------------------------------------------------ */
/* 2. RaMeshMultiGenomeGraph – ctor                                   */
/* ------------------------------------------------------------------ */
    RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(
       std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map)
    {
        for (const auto& [sp, mgr_var] : seqpro_map)
        {
            /* 取出染色体名称列表 */
            std::vector<ChrName> chr_names = std::visit(
                [](const auto& up) -> std::vector<ChrName>
                {
                    return up ? up->getSequenceNames()
                        : std::vector<ChrName>{};
                },
                mgr_var);

            /* 原地构造 <key, RaMeshGenomeGraph> */
            species_graphs.try_emplace(sp,           // map key
                sp,           // RaMeshGenomeGraph::species_name
                chr_names);   // vector<ChrName>
        }
    }


    /* ------------------------------------------------------------------ */
    /* 3. insertClusterIntoGraph                                          */
    /* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* insertClusterIntoGraph – 节点级写锁版本                            */
/* ------------------------------------------------------------------ */
    void RaMeshMultiGenomeGraph::insertClusterIntoGraph(
        SpeciesName         ref_name,
        SpeciesName         qry_name,
        const MatchCluster& cluster)
    {
        if (cluster.empty()) return;

        const ChrName& ref_chr = cluster.front().ref_region.chr_name;
        const ChrName& qry_chr = cluster.front().query_region.chr_name;

        auto& refEnd = species_graphs[ref_name].chr2end[ref_chr];
        auto& qryEnd = species_graphs[qry_name].chr2end[qry_chr];

        std::scoped_lock chromWrite(refEnd.rw, qryEnd.rw);

        /* 1. 找包围区间（只读，无锁即可） */
        uint_t ref_beg = cluster.front().ref_region.start;
        uint_t ref_end = cluster.back().ref_region.start
            + cluster.back().ref_region.length;
        uint_t qry_beg = cluster.front().query_region.start;
        uint_t qry_end = cluster.back().query_region.start
            + cluster.back().query_region.length;

        auto [ref_prevAtom, ref_nextAtom] =
            refEnd.findSurroundingSegmentRange(ref_beg, ref_end);
        auto [qry_prevAtom, qry_nextAtom] =
            qryEnd.findSurroundingSegmentRange(qry_beg, qry_end);

        SegAtomPtr lastRefAtom = ref_prevAtom;
        SegAtomPtr lastQryAtom = qry_prevAtom;

        /* 2. 逐 match 插入（对局部节点加写锁） */
        for (const auto& m : cluster)
        {
            BlockPtr blk = Block::make(2);
            blk->ref_chr = ref_chr;

            /* --- 构建两条新 Segment --- */
            SegPtr refSeg = Segment::create_from_region(
                const_cast<Region&>(m.ref_region), m.strand,
                Cigar_t{}, AlignRole::PRIMARY,
                SegmentRole::SEGMENT, blk);
            SegAtomPtr refAtom = make_atom(refSeg);

            SegPtr qrySeg = Segment::create_from_region(
                const_cast<Region&>(m.query_region), m.strand,
                Cigar_t{}, AlignRole::PRIMARY,
                SegmentRole::SEGMENT, blk);
            SegAtomPtr qryAtom = make_atom(qrySeg);

            blk->anchors[{ref_name, ref_chr}] = refAtom;
            blk->anchors[{qry_name, qry_chr}] = qryAtom;

            /* ===== 插入 ref 链 ===== */
            {
                refSeg->primary_path.prev = lastRefAtom;
                refSeg->primary_path.next = ref_nextAtom;

                lastRefAtom->load()->primary_path.next = refAtom;
                if (ref_nextAtom)
                    ref_nextAtom->load()->primary_path.prev = refAtom;
            }
            lastRefAtom = refAtom;

            /* ===== 插入 qry 链 ===== */
            {
                qrySeg->primary_path.prev = lastQryAtom;
                qrySeg->primary_path.next = qry_nextAtom;

                lastQryAtom->load()->primary_path.next = qryAtom;
                if (qry_nextAtom)
                    qry_nextAtom->load()->primary_path.prev = qryAtom;
            }
            lastQryAtom = qryAtom;

            /* ---- 记录 block ---- */
            {
                std::unique_lock poolLock(rw);
                blocks.emplace_back(WeakBlock(blk));
            }
        }
    }


    /* ------------------------------------------------------------------ */
    /* 4. debug_print (multi-genome)                                       */
    /* ------------------------------------------------------------------ */
    void RaMeshMultiGenomeGraph::debug_print(std::ostream& os,
        bool show_sec) const
    {
        std::shared_lock gLock(rw);

        os << "\n********  Multi-Genome Graph  ********\n";
        for (const auto& [sp, g] : species_graphs)
            g.debug_print(os, show_sec);
        os << "********  End of Graphs  ********\n";
    }

} // namespace RaMesh
