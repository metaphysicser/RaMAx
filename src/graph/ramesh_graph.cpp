// ******************************************************
//  ramesh.cpp – 适配 ramesh.h (v0.5-beta)
//  ✧  无锁并行插入版本  ✧
//  2025-06-24
// ******************************************************
#include "ramesh.h"
#include <iomanip>
#include <shared_mutex>

namespace RaMesh {


    /* ======================================================
     * 1. RaMeshGenomeGraph
     * ====================================================*/
    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName& sp)
        : species_name(sp) {
    }

    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName& sp,
        const std::vector<ChrName>& chrs)
        : species_name(sp)
    {
        chr2end.reserve(chrs.size());
        for (const auto& c : chrs)
            chr2end.try_emplace(c);
    }

    size_t RaMeshGenomeGraph::debug_print(std::ostream& os, bool /*show_secondary*/) const
    {
        std::shared_lock lg(rw);
        size_t total_idx = 0;
        os << "=== GenomeGraph <" << species_name << "> ===\n";
        for (const auto& [chr, end] : chr2end)
        {
            os << "\n[Chromosome " << chr << "]\n";
            os << std::left << std::setw(6) << "Idx"
                << std::setw(12) << "Start"
                << std::setw(10) << "Len"
                << std::setw(4) << "Str"
                << std::setw(6) << "Role" << '\n';

            SegPtr curr = end.head->primary_path.next.load(std::memory_order_acquire);
            size_t idx = 0;
            while (curr && !curr->is_tail())
            {
                os << std::setw(6) << idx++
                    << std::setw(12) << curr->start
                    << std::setw(10) << curr->length
                    << std::setw(4) << (curr->strand == Strand::FORWARD ? "+" : "-")
                    << std::setw(6) << (curr->is_primary() ? "Pri" : "Sec")
                    << '\n';
                curr = curr->primary_path.next.load(std::memory_order_acquire);
            }
            os << "Total segments: " << idx << '\n';
            total_idx += idx;
        }
        os << "=== End of Graph ===\n";
        return total_idx;
    }

    /* ======================================================
     * 2. RaMeshMultiGenomeGraph – ctor
     * ====================================================*/
    RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(
        std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map)
    {
        for (const auto& [sp, mgr_var] : seqpro_map)
        {
            std::vector<ChrName> chr_names = std::visit(
                [](const auto& up)->std::vector<ChrName> {
                    return up ? up->getSequenceNames() : std::vector<ChrName>{};
                }, mgr_var);
            species_graphs.try_emplace(sp, sp, chr_names);
        }
    }
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

        /* --------------------------------------------------
         * 逐 match 构建 Block 和 Segment, 并使用无锁插入
         * ------------------------------------------------*/
        for (const auto& m : cluster)
        {
            /* 创建 block, 预留两个 anchor */
            BlockPtr blk = Block::make(2);
            blk->ref_chr = ref_chr;

            /* —— ref Segment —— */
            SegPtr refSeg = Segment::create_from_region(
                const_cast<Region&>(m.ref_region), m.strand,
                Cigar_t{}, AlignRole::PRIMARY, SegmentRole::SEGMENT, blk);

            /* —— qry Segment —— */
            SegPtr qrySeg = Segment::create_from_region(
                const_cast<Region&>(m.query_region), m.strand,
                Cigar_t{}, AlignRole::PRIMARY, SegmentRole::SEGMENT, blk);

            /* 记录锚点到 block */
            {
                std::unique_lock lk(blk->rw);
                blk->anchors[{ref_name, ref_chr}] = refSeg;
                blk->anchors[{qry_name, qry_chr}] = qrySeg;
            }

            /* —— 插入主链（无锁 CAS）—— */
            concurrent_insert_segment(refEnd.head, refSeg);
            concurrent_insert_segment(qryEnd.head, qrySeg);

            /* —— 记录到全局 block 池 —— */
            {
                std::unique_lock poolLock(rw);
                blocks.emplace_back(WeakBlock(blk));
            }
        }
    }


    /* ======================================================
     * 4. debug_print (multi-genome)
     * ====================================================*/
    void RaMeshMultiGenomeGraph::debug_print(std::ostream& os, bool show_sec) const
    {
        std::shared_lock gLock(rw);

        std::vector<size_t> seg_num_vec;
        os << "\n********  Multi-Genome Graph  ********\n";
        for (const auto& [sp, g] : species_graphs) {
            size_t seg_num = g.debug_print(os, show_sec);
            seg_num_vec.push_back(seg_num);
        }

        size_t count = 0;
        for (const auto& [sp, g] : species_graphs) {
            // 打印每个图的segment数量
            os << seg_num_vec[count] << std::endl;
            count++;
        }
            
       
        os << "********  End of Graphs  ********\n";
    }

} // namespace RaMesh
