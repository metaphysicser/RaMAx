// ramesh.cpp  – 适配 ramesh.h (v0.4)
#include "ramesh.h"

namespace RaMesh {
    /* ------------------------------------------------------------------ */
    /* RaMeshGenomeGraph / MultiGenomeGraph                               */
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

    RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(
        std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map)
    {
        for (auto& [sp, mgr_var] : seqpro_map) {
            std::vector<ChrName> chr_names = std::visit(
                [](auto& up)->std::vector<ChrName> {
                    return up ? up->getSequenceNames() : std::vector<ChrName>{};
                }, mgr_var);
            species_graphs.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(sp),
                std::forward_as_tuple(sp, chr_names));
        }
    }

    /* ------------------------------------------------------------------ */
    /* insertClusterIntoGraph                                             */
    /* ------------------------------------------------------------------ */
    void RaMeshMultiGenomeGraph::insertClusterIntoGraph(
        SpeciesName         ref_name,
        SpeciesName         query_name,
        const MatchCluster& cluster)
    {
        if (cluster.empty()) return;

        const ChrName& ref_chr = cluster.front().ref_region.chr_name;
        const ChrName& qry_chr = cluster.front().query_region.chr_name;

        uint_t ref_start = cluster.front().ref_region.start;
        uint_t ref_end = cluster.back().ref_region.start +
            cluster.back().ref_region.length;

        uint_t qry_start = cluster.front().query_region.start;
        uint_t qry_end = cluster.back().query_region.start +
            cluster.back().query_region.length;

        auto& refEnd = species_graphs[ref_name].chr2end[ref_chr];
        auto& qryEnd = species_graphs[query_name].chr2end[qry_chr];

        auto [ref_prev, ref_next] =
            refEnd.findSurroundingSegmentRange(ref_start, ref_end);
        auto [qry_prev, qry_next] =
            qryEnd.findSurroundingSegmentRange(qry_start, qry_end);

        /* 把 prev 节点的 SegAtomPtr 取出备用 */
        auto getAtomPtr = [](SegPtr node)->SegAtomPtr {
            if (!node) return nullptr;
            auto* base = reinterpret_cast<RaMeshNode*>(node);
            // head/tail/segment 均有 primary_path.next，可取其 prev/next 来获得 SegAtomPtr。
            // 若 prev/next 为 nullptr，则需要在插入时重新创建新的 SegAtomPtr。
            return base->primary_path.next; // 只为示例；真正使用时按需要取 prev / next
            };

        SegAtomPtr lastRefAtom =
            (ref_prev == reinterpret_cast<SegPtr>(refEnd.head.get()))
            ? refEnd.head->primary_path.next
            : std::make_shared<SegAtom>(ref_prev);

        SegAtomPtr lastQryAtom =
            (qry_prev == reinterpret_cast<SegPtr>(qryEnd.head.get()))
            ? qryEnd.head->primary_path.next
            : std::make_shared<SegAtom>(qry_prev);

        /* ------ 逐 Match 插入，保持顺序 ------ */
        for (const auto& m : cluster)
        {
            BlockPtr blk = Block::make(2);
            blk->ref_chr = ref_chr;

            /* —— ref side —— */
            SegPtr refSeg = Segment::create_from_region(
                const_cast<Region&>(m.ref_region),
                m.strand,
                Cigar_t{},
                AlignRole::PRIMARY,
                blk);
            SegAtomPtr refAtom = std::make_shared<SegAtom>(refSeg);

            /* —— query side —— */
            SegPtr qrySeg = Segment::create_from_region(
                const_cast<Region&>(m.query_region),
                m.strand,
                Cigar_t{},
                AlignRole::SECONDARY,
                blk);
            SegAtomPtr qryAtom = std::make_shared<SegAtom>(qrySeg);

            // 双向 secondary_link
            refSeg->secondary_path.next = qryAtom;
            qrySeg->secondary_path.prev = refAtom;

            /* ------ 插入 ref 链尾 ------ */
            auto* prevRefNode = reinterpret_cast<RaMeshNode*>(lastRefAtom->load());
            prevRefNode->primary_path.next = refAtom;
            refSeg->primary_path.prev = lastRefAtom;
            lastRefAtom = refAtom;                       // 向后推进

            /* ------ 插入 query 链尾 ------ */
            auto* prevQryNode = reinterpret_cast<RaMeshNode*>(lastQryAtom->load());
            prevQryNode->primary_path.next = qryAtom;
            qrySeg->primary_path.prev = lastQryAtom;
            lastQryAtom = qryAtom;                       // 向后推进

            /* ------ 登记到 Block anchors ------ */
            blk->anchors[{ref_name, ref_chr }] = refAtom;
            blk->anchors[{query_name, qry_chr }] = qryAtom;

            /* ------ 全局 Block 弱引用 ------ */
            {
                std::lock_guard<std::mutex> gLock(rw);
                blocks.emplace_back(WeakBlock(blk));
            }
        }

        /* ------ 接尾巴与 next 哨兵/节点 ------ */
        if (ref_next) {
            auto* last = reinterpret_cast<RaMeshNode*>(lastRefAtom->load());
            last->primary_path.next = std::make_shared<SegAtom>(ref_next);
            reinterpret_cast<RaMeshNode*>(ref_next)->primary_path.prev =
                lastRefAtom;
        }
        if (qry_next) {
            auto* last = reinterpret_cast<RaMeshNode*>(lastQryAtom->load());
            last->primary_path.next = std::make_shared<SegAtom>(qry_next);
            reinterpret_cast<RaMeshNode*>(qry_next)->primary_path.prev =
                lastQryAtom;
        }
    }

} // namespace RaMesh
