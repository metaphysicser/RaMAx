// =============================================================
//  File: ramesh_graph.cpp –  High‑level graph ops (v0.6‑alpha)
// =============================================================
#include "ramesh.h"
#include <iomanip>
#include <shared_mutex>

namespace RaMesh {

    /* =============================================================
     * 1.  RaMeshGenomeGraph
     * ===========================================================*/
    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName& sp)
        : species_name(sp) {
    }

    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName& sp,
        const std::vector<ChrName>& chrs)
        : species_name(sp)
    {
        chr2end.reserve(chrs.size());
        for (const auto& c : chrs) chr2end.try_emplace(c);
    }

    size_t RaMeshGenomeGraph::debugPrint(bool show_detail) const
    {
        std::shared_lock lg(rw);
        size_t total_idx = 0;

        std::cout << "=== GenomeGraph <" << species_name << "> ===\n";

        for (const auto& [chr, end] : chr2end) {
            SegPtr first = end.head->primary_path.next.load(std::memory_order_acquire);

            /*----------- 统计这一条染色体的段数 -----------*/
            size_t chr_count = 0;
            for (SegPtr p = first; p && !p->isTail();
                p = p->primary_path.next.load(std::memory_order_acquire))
                ++chr_count;

            std::cout << "\n[Chromosome " << chr << "]\n";

            if (show_detail) {
                /*----------- 打印表头 -----------*/
                std::cout << std::left
                    << std::setw(6) << "Idx"
                    << std::setw(12) << "Start"
                    << std::setw(12) << "End"
                    << std::setw(10) << "Len"
                    << std::setw(4) << "Str"
                    << std::setw(6) << "Role"
                    << std::setw(14) << "Ptr"
                    << '\n';

                /*----------- 打印明细 -----------*/
                size_t idx = 0;
                for (SegPtr p = first; p && !p->isTail();
                    p = p->primary_path.next.load(std::memory_order_acquire)) {
                    std::cout << std::setw(6) << idx++
                        << std::setw(12) << p->start
                        << std::setw(12) << (p->start + p->length - 1)
                        << std::setw(10) << p->length
                        << std::setw(4) << (p->strand == Strand::FORWARD ? "+" : "-")
                        << std::setw(6) << (p->isPrimary() ? "Pri" : "Sec")
                        << '\n';
                }
            }

            std::cout << "Total segments: " << chr_count << '\n';
            total_idx += chr_count;
        }

        std::cout << "=== End of Graph ===\n";
        return total_idx;
    }


    /* =============================================================
     * 2.  RaMeshMultiGenomeGraph – ctor
     * ===========================================================*/
    RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map)
    {
        for (const auto& [sp, mgr_var] : seqpro_map) {
            std::vector<ChrName> chr_names = std::visit(
                [](const auto& up) -> std::vector<ChrName> {
                    return up ? up->getSequenceNames() : std::vector<ChrName>{};
                }, mgr_var);
            species_graphs.try_emplace(sp, sp, chr_names);
        }
    }

    /* =============================================================
     * 3.  Cluster insertion (public API)
     * ===========================================================*/
    void RaMeshMultiGenomeGraph::insertClusterIntoGraph(SpeciesName ref_name,
        SpeciesName qry_name,
        const MatchCluster& cluster)
    {
        if (cluster.empty()) return;

        // 1. Locate ends for reference & query chromosomes
        const ChrName& ref_chr = cluster.front().ref_region.chr_name;
        const ChrName& qry_chr = cluster.front().query_region.chr_name;

        auto& ref_end = species_graphs[ref_name].chr2end[ref_chr];
        auto& qry_end = species_graphs[qry_name].chr2end[qry_chr];

        // 2. Build blocks & segments once (no shared‑state yet)
        std::vector<SegPtr> ref_segs; ref_segs.reserve(cluster.size());
        std::vector<SegPtr> qry_segs; qry_segs.reserve(cluster.size());

        for (const auto& m : cluster) {
            BlockPtr blk = Block::create(2);
            blk->ref_chr = ref_chr;

            auto [r_seg, q_seg] = Block::createSegmentPair(m, ref_name, qry_name, ref_chr, qry_chr, blk);

            ref_segs.emplace_back(r_seg);
            qry_segs.emplace_back(q_seg);

            // register to global pool
            {
                std::unique_lock pool_lock(rw);
                blocks.emplace_back(WeakBlock(blk));
            }
        }

        // 3. Link internal chains locally (single‑threaded)
        Segment::linkChain(ref_segs);
        Segment::linkChain(qry_segs);

        // 4. Atomically splice into genome graph
        uint_t ref_beg = cluster.front().ref_region.start;
        uint_t ref_end_pos = cluster.back().ref_region.start + cluster.back().ref_region.length;
        uint_t qry_beg = cluster.front().query_region.start;
        uint_t qry_end_pos = cluster.back().query_region.start + cluster.back().query_region.length;

        ref_end.spliceSegmentChain(ref_segs, ref_beg, ref_end_pos);
        qry_end.spliceSegmentChain(qry_segs, qry_beg, qry_end_pos);
    }

    /* =============================================================
 * 3.  Anchor insertion (public API)
 * ===========================================================*/
    void RaMeshMultiGenomeGraph::insertAnchorIntoGraph(SpeciesName ref_name, SpeciesName qry_name,
        const AnchorVec& anchor_vec)
    {
        if (anchor_vec.empty()) return;

        // 1. Locate ends for reference & query chromosomes
        const ChrName& ref_chr = anchor_vec.front().match.ref_region.chr_name;
        const ChrName& qry_chr = anchor_vec.front().match.query_region.chr_name;

        auto& ref_end = species_graphs[ref_name].chr2end[ref_chr];
        auto& qry_end = species_graphs[qry_name].chr2end[qry_chr];

        // 2. Build blocks & segments once (no shared‑state yet)
        std::vector<SegPtr> ref_segs; ref_segs.reserve(anchor_vec.size());
        std::vector<SegPtr> qry_segs; qry_segs.reserve(anchor_vec.size());

        for (const Anchor& m : anchor_vec) {
            BlockPtr blk = Block::make(2);
            blk->ref_chr = ref_chr;

            auto [r_seg, q_seg] = Block::createSegmentPair(m, ref_name, qry_name, ref_chr, qry_chr, blk);

            ref_segs.emplace_back(r_seg);
            qry_segs.emplace_back(q_seg);

            // register to global pool
            {
                std::unique_lock pool_lock(rw);
                blocks.emplace_back(WeakBlock(blk));
            }
        }

        // 3. Link internal chains locally (single‑threaded)
        Segment::linkChain(ref_segs);
        Segment::linkChain(qry_segs);

        // 4. Atomically splice into genome graph
        uint_t ref_beg = anchor_vec.front().match.ref_region.start;
        uint_t ref_end_pos = anchor_vec.back().match.ref_region.start + anchor_vec.back().match.ref_region.length;
        uint_t qry_beg = anchor_vec.front().match.query_region.start;
        uint_t qry_end_pos = anchor_vec.back().match.query_region.start + anchor_vec.back().match.query_region.length;

        ref_end.spliceSegmentChain(ref_segs, ref_beg, ref_end_pos);
        qry_end.spliceSegmentChain(qry_segs, qry_beg, qry_end_pos);
    }

    /* ==============================================================
 * 4.  debugPrint (multi-genome)  -- 新版，参数改为 show_detail
 * ==============================================================*/
    void RaMeshMultiGenomeGraph::debugPrint(bool show_detail) const
    {
        std::shared_lock gLock(rw);

        /*------------- 页眉 -------------*/
        std::cout << "\n********  Multi-Genome Graph  ********\n";

        /*------------- 逐物种打印 + 计数 -------------*/
        std::vector<std::pair<std::string, size_t>> per_species;  // {species, seg_cnt}
        size_t grand_total = 0;

        for (const auto& [sp, g] : species_graphs) {
            // 假设 g.debugPrint 有重载：size_t debugPrint(std::ostream&, bool show_detail) const
            size_t seg_cnt = g.debugPrint(show_detail);
            per_species.emplace_back(sp, seg_cnt);
            grand_total += seg_cnt;
        }

        /*------------- 汇总区 -------------*/
        std::cout << "\n----------  Summary  ----------\n";
        for (const auto& [sp, cnt] : per_species) {
            std::cout << std::left << std::setw(15) << sp << ": "
                << cnt << " segments\n";
        }
        std::cout << "--------------------------------\n";
        std::cout << "Grand total  : " << grand_total
            << " segments in " << per_species.size() << " genome(s)\n";

        std::cout << "********  End of Graphs  ********\n";
    }

} // namespace RaMesh
