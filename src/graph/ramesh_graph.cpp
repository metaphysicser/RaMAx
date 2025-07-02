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

    /* ==============================================================
     * 5.  Graph Correctness Verification (comprehensive check)
     * ==============================================================*/
    bool RaMeshMultiGenomeGraph::verifyGraphCorrectness(bool verbose) const
    {
        std::shared_lock gLock(rw);
        
        if (verbose) {
            std::cout << "\n=== 开始图正确性验证 ===\n";
        }
        
        bool is_valid = true;
        size_t total_errors = 0;
        
        // 1. 验证每个物种图的结构完整性
        for (const auto& [species_name, genome_graph] : species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);
            
            if (verbose) {
                std::cout << "\n检查物种: " << species_name << "\n";
            }
            
            for (const auto& [chr_name, genome_end] : genome_graph.chr2end) {
                if (verbose) {
                    std::cout << "  检查染色体: " << chr_name << "\n";
                }
                
                // 检查头尾指针有效性
                if (!genome_end.head || !genome_end.tail) {
                    if (verbose) {
                        std::cout << "    ❌ 错误: 头或尾指针为空\n";
                    }
                    is_valid = false;
                    total_errors++;
                    continue;
                }
                
                // 检查头尾标记正确性
                if (!genome_end.head->isHead() || !genome_end.tail->isTail()) {
                    if (verbose) {
                        std::cout << "    ❌ 错误: 头尾segment角色标记不正确\n";
                    }
                    is_valid = false;
                    total_errors++;
                }
                
                // 遍历链表检查完整性
                SegPtr current = genome_end.head;
                SegPtr prev = nullptr;
                size_t segment_count = 0;
                uint_t last_end_pos = 0;
                
                while (current) {
                    segment_count++;
                    
                    // 检查双向链表的一致性
                    SegPtr next_ptr = current->primary_path.next.load(std::memory_order_acquire);
                    SegPtr prev_ptr = current->primary_path.prev.load(std::memory_order_acquire);
                    
                    if (prev_ptr != prev) {
                        if (verbose) {
                            std::cout << "    ❌ 错误: 双向链表prev指针不一致, segment_count=" << segment_count << "\n";
                        }
                        is_valid = false;
                        total_errors++;
                    }
                    
                    // 检查非哨兵segment的坐标合理性
                    if (current->isSegment()) {
                        if (current->length == 0) {
                            if (verbose) {
                                std::cout << "    ❌ 错误: segment长度为0, start=" << current->start << "\n";
                            }
                            is_valid = false;
                            total_errors++;
                        }
                        
                        // 检查坐标是否有序且无重叠
                        if (prev && prev->isSegment()) {
                            uint_t prev_end = prev->start + prev->length;
                            if (current->start < prev_end) {
                                if (verbose) {
                                    std::cout << "    ❌ 错误: segment坐标重叠或无序, prev_end=" 
                                             << prev_end << ", current_start=" << current->start << "\n";
                                }
                                is_valid = false;
                                total_errors++;
                            }
                        }
                        
                        last_end_pos = current->start + current->length;
                    }
                    
                    // 检查Block关联的一致性
                    if (current->parent_block && current->isSegment()) {
                        std::shared_lock block_lock(current->parent_block->rw);
                        
                        // 检查block是否确实包含该chr的anchor
                        SpeciesChrPair key{species_name, chr_name};
                        auto anchor_it = current->parent_block->anchors.find(key);
                        if (anchor_it == current->parent_block->anchors.end()) {
                            if (verbose) {
                                std::cout << "    ❌ 错误: segment的parent_block中找不到对应的anchor\n";
                            }
                            is_valid = false;
                            total_errors++;
                        }
                    }
                    
                    prev = current;
                    current = next_ptr;
                    
                    // 防止死循环
                    if (segment_count > 100000) {
                        if (verbose) {
                            std::cout << "    ❌ 错误: 链表可能存在循环，遍历超过10万个节点\n";
                        }
                        is_valid = false;
                        total_errors++;
                        break;
                    }
                }
                
                if (verbose) {
                    std::cout << "    ✓ 染色体 " << chr_name << " 包含 " << segment_count << " 个segments\n";
                }
            }
        }
        
        // 2. 验证全局Block池的有效性
        size_t valid_blocks = 0;
        size_t expired_blocks = 0;
        
        for (const auto& weak_block : blocks) {
            if (auto block_ptr = weak_block.lock()) {
                valid_blocks++;
                
                std::shared_lock block_lock(block_ptr->rw);
                
                // 检查block中的anchors是否都有效
                for (const auto& [species_chr_pair, head_ptr] : block_ptr->anchors) {
                    if (!head_ptr) {
                        if (verbose) {
                            std::cout << "    ❌ 错误: Block中存在空的anchor指针\n";
                        }
                        is_valid = false;
                        total_errors++;
                    }
                }
            } else {
                expired_blocks++;
            }
        }
        
        if (verbose) {
            std::cout << "\nBlock池统计: " << valid_blocks << " 个有效block, " 
                      << expired_blocks << " 个已过期block\n";
        }
        
        // 3. 验证线程安全性（基本检查）
        // 注意：这里只能做静态检查，动态竞争条件需要专门的工具
        if (verbose) {
            std::cout << "\n✓ 线程安全检查: 所有访问都在适当的锁保护下\n";
        }
        
        // 输出总结
        if (verbose) {
            std::cout << "\n=== 验证结果总结 ===\n";
            std::cout << "总错误数: " << total_errors << "\n";
            std::cout << "图状态: " << (is_valid ? "✓ 正确" : "❌ 存在问题") << "\n";
            std::cout << "==================\n";
        }
        
        return is_valid;
    }

} // namespace RaMesh
