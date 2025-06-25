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
    /* ---------- 把 [first,last] 原子挂到 (prev,next) 之间 ---------- */
    static bool splice_range(SegPtr prev, SegPtr next,
        SegPtr first, SegPtr last)
    {
        /* Step-0: 先写内部指针 */
        first->primary_path.prev.store(prev, std::memory_order_relaxed);
        last->primary_path.next.store(next, std::memory_order_relaxed);
        std::atomic_thread_fence(std::memory_order_release);

        /* Step-1: CAS 修改 prev->next */
        SegPtr expected = next;
        if (!prev->primary_path.next.compare_exchange_strong(
            expected, first,
            std::memory_order_acq_rel,
            std::memory_order_acquire))
        {
            // 失败，caller 需要重试
            return false;
        }

        /* Step-2: 修补 next->prev */
        next->primary_path.prev.store(last, std::memory_order_release);
        return true;
    }


// 用于创建和链接 Block 和 Segment
std::tuple<SegPtr, SegPtr> create_segment_from_match(const Match& m, const SpeciesName& ref_name, const SpeciesName& qry_name, const ChrName& ref_chr, const ChrName& qry_chr, BlockPtr& blk)
{
    // 创建 ref Segment
    SegPtr refSeg = Segment::create_from_region(
        const_cast<Region&>(m.ref_region), m.strand,
        Cigar_t{}, AlignRole::PRIMARY, SegmentRole::SEGMENT, blk);

    // 创建 qry Segment
    SegPtr qrySeg = Segment::create_from_region(
        const_cast<Region&>(m.query_region), m.strand,
        Cigar_t{}, AlignRole::PRIMARY, SegmentRole::SEGMENT, blk);

    // 记录锚点到 block
    {
        std::unique_lock lk(blk->rw);
        blk->anchors[{ref_name, ref_chr}] = refSeg;
        blk->anchors[{qry_name, qry_chr}] = qrySeg;
    }

    return { refSeg, qrySeg };
}

// 用于链接一条链表（更新前一个和下一个元素的双向指针）
void link_chain(std::vector<SegPtr>& vec)
{
    for (size_t i = 0; i + 1 < vec.size(); ++i) {
        vec[i]->primary_path.next.store(vec[i + 1], std::memory_order_relaxed);
        vec[i + 1]->primary_path.prev.store(vec[i], std::memory_order_relaxed);
    }
}

/* ---------- 原子地把整段 splice 进去 ---------- */
static void splice_cluster_into_graph(GenomeEnd& end,
    std::vector<SegPtr>& vec,
    uint_t beg, uint_t ed)
{
    for (;;)
    {
        auto [prev, next] = end.find_surrounding(beg, ed);

        if (splice_range(prev, next, vec.front(), vec.back()))
            break;          // 真正成功才退出
        /* 否则继续 while 循环，重新定位 prev/next 后重试 */
    }
}

// 主要插入函数
void RaMeshMultiGenomeGraph::insertClusterIntoGraph(
    SpeciesName ref_name, SpeciesName qry_name,
    const MatchCluster& cluster)
{
    if (cluster.empty()) return;

    const ChrName& ref_chr = cluster.front().ref_region.chr_name;
    const ChrName& qry_chr = cluster.front().query_region.chr_name;

    auto& refEnd = species_graphs[ref_name].chr2end[ref_chr];
    auto& qryEnd = species_graphs[qry_name].chr2end[qry_chr];

    /* ---------- 1. 先一次性构建所有 Block & Segment ----------- */
    std::vector<SegPtr> refSegs; refSegs.reserve(cluster.size());
    std::vector<SegPtr> qrySegs; qrySegs.reserve(cluster.size());

    for (const auto& m : cluster)
    {
        BlockPtr blk = Block::make(2);
        blk->ref_chr = ref_chr;

        // 创建并获取 refSeg 和 qrySeg
        auto [rs, qs] = create_segment_from_match(m, ref_name, qry_name, ref_chr, qry_chr, blk);

        refSegs.emplace_back(rs);
        qrySegs.emplace_back(qs);

        // 记录到全局 block 池
        {
            std::unique_lock poolLock(rw);
            blocks.emplace_back(WeakBlock(blk));
        }
    }

    /* ---------- 2. 把每条链内部先串好 ---------- */
    link_chain(refSegs);
    link_chain(qrySegs);

    /* ---------- 3. 原子地把整段 splice 进去 ---------- */
    uint_t ref_beg = cluster.front().ref_region.start;
    uint_t ref_end = cluster.back().ref_region.start + cluster.back().ref_region.length;
    uint_t qry_beg = cluster.front().query_region.start;
    uint_t qry_end = cluster.back().query_region.start + cluster.back().query_region.length;

    splice_cluster_into_graph(refEnd, refSegs, ref_beg, ref_end);
    splice_cluster_into_graph(qryEnd, qrySegs, qry_beg, qry_end);
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
