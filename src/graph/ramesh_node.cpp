// =============================================================
//  File: ramesh.cpp   –  Low‑level primitives (v0.6‑alpha)
//  ✧  Segment / Block / GenomeEnd implementation  ✧
// =============================================================
#include "ramesh.h"
#include <iomanip>
#include <shared_mutex>

namespace RaMesh {

    /* =============================================================
     * 0.  Segment factories & utilities
     * ===========================================================*/
    SegPtr Segment::create(uint_t start, uint_t len, Strand sd,
        Cigar_t cg, AlignRole rl, SegmentRole sl,
        const BlockPtr& bp)
    {
        auto* s = new Segment();
        s->start = start;
        s->length = len;
        s->strand = sd;
        s->cigar = std::move(cg);
        s->align_role = rl;
        s->seg_role = sl;
        if (bp) s->parent_block = bp;

        // list pointers – initialise nullptr
        s->primary_path.next.store(nullptr, std::memory_order_relaxed);
        s->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return s;
    }

    SegPtr Segment::createFromRegion(Region& region, Strand sd,
        Cigar_t cg, AlignRole rl, SegmentRole sl,
        const BlockPtr& bp)
    {
        return create(region.start, region.length, sd, std::move(cg), rl, sl, bp);
    }

    SegPtr Segment::createHead()
    {
        auto* h = new Segment();
        h->seg_role = SegmentRole::HEAD;
        h->align_role = AlignRole::PRIMARY;
        h->primary_path.next.store(nullptr, std::memory_order_relaxed);
        h->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return h;
    }

    SegPtr Segment::createTail()
    {
        auto* t = new Segment();
        t->seg_role = SegmentRole::TAIL;
        t->align_role = AlignRole::PRIMARY;
        t->primary_path.next.store(nullptr, std::memory_order_relaxed);
        t->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return t;
    }

    void Segment::linkChain(const std::vector<SegPtr>& segs)
    {
        for (size_t i = 0; i + 1 < segs.size(); ++i) {
            segs[i]->primary_path.next.store(segs[i + 1], std::memory_order_relaxed);
            segs[i + 1]->primary_path.prev.store(segs[i], std::memory_order_release);
        }
    }
    
    void Segment::unlinkSegment(SegPtr segment) {
        if (!segment || segment->isHead() || segment->isTail()) return;
        
        // 获取前驱和后继
        SegPtr prev = segment->primary_path.prev.load(std::memory_order_acquire);
        SegPtr next = segment->primary_path.next.load(std::memory_order_acquire);
        
        if (prev && next) {
            // 原子地更新链表指针
            prev->primary_path.next.store(next, std::memory_order_release);
            next->primary_path.prev.store(prev, std::memory_order_release);
            
            // 清空被删除segment的指针
            segment->primary_path.next.store(nullptr, std::memory_order_relaxed);
            segment->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        }
    }
    
    void Segment::deleteSegment(SegPtr segment) {
        if (!segment) return;
        
        // 先从链表中解除链接
        unlinkSegment(segment);
        
        // 清理parent_block引用
        if (segment->parent_block) {
            std::unique_lock block_lock(segment->parent_block->rw);
            
            // 从block的anchors中移除此segment
            for (auto it = segment->parent_block->anchors.begin(); 
                 it != segment->parent_block->anchors.end(); ++it) {
                if (it->second == segment) {
                    segment->parent_block->anchors.erase(it);
                    break;
                }
            }
            segment->parent_block.reset();
        }
        
        // 释放内存 - segment是原生指针，需要手动删除
        delete segment;
    }
    
    void Segment::deleteBatch(const std::vector<SegPtr>& segments) {
        if (segments.empty()) return;
        
        // 按block分组以减少锁竞争
        std::unordered_map<BlockPtr, std::vector<SegPtr>> block_groups;
        std::vector<SegPtr> orphaned_segments;
        
        for (SegPtr seg : segments) {
            if (!seg || seg->isHead() || seg->isTail()) continue;
            
            if (seg->parent_block) {
                block_groups[seg->parent_block].emplace_back(seg);
            } else {
                orphaned_segments.emplace_back(seg);
            }
        }
        
        // 批量处理每个block
        for (auto& [block, segs] : block_groups) {
            std::unique_lock block_lock(block->rw);
            
            for (SegPtr seg : segs) {
                // 从链表中解除链接
                unlinkSegment(seg);
                
                // 从anchors中移除
                for (auto it = block->anchors.begin(); it != block->anchors.end(); ++it) {
                    if (it->second == seg) {
                        block->anchors.erase(it);
                        break;
                    }
                }
                seg->parent_block.reset();
            }
        }
        
        // 处理孤立的segments
        for (SegPtr seg : orphaned_segments) {
            unlinkSegment(seg);
        }
        
        // 释放内存
        for (SegPtr seg : segments) {
            if (seg && !seg->isHead() && !seg->isTail()) {
                delete seg;
            }
        }
    }

    /* =============================================================
     * 1. GenomeEnd helpers
     * ===========================================================*/
    GenomeEnd::GenomeEnd() {
        head_holder.reset(Segment::createHead());
        tail_holder.reset(Segment::createTail());
        head = head_holder.get();
        tail = tail_holder.get();
        head->primary_path.next.store(tail, std::memory_order_relaxed);
        tail->primary_path.prev.store(head, std::memory_order_release);

        sample_vec.resize(1, head);   // slot 0 永远指向 head
    }

    /* ---------- 采样表维护 ---------- */
    void GenomeEnd::ensureSampleSize(uint_t pos) {
        std::size_t need = pos / kSampleStep + 1;
        if (need > sample_vec.size()) sample_vec.resize(need, nullptr);
    }

    void GenomeEnd::updateSampling(const std::vector<SegPtr>& segs) {
        if (segs.empty()) return;
        std::unique_lock lk(rw);      // ××× 写锁 —— 修改 sample_vec
        for (SegPtr s : segs) {
            std::size_t idx = s->start / kSampleStep;
            ensureSampleSize(s->start);
            if (!sample_vec[idx] || s->start < sample_vec[idx]->start)
                sample_vec[idx] = s;   // 只保留区间内最左端 segment
        }
    }

    std::pair<SegPtr, SegPtr>
        GenomeEnd::findSurrounding(uint_t range_start, uint_t range_end) {
        // 1) 读取采样表得到“最近前驱”的 hint
        std::shared_lock lk(rw);                 // 读锁即可
        std::size_t slot = range_start / kSampleStep;
        SegPtr hint = (slot < sample_vec.size() && sample_vec[slot])
            ? sample_vec[slot]
            : head;
            lk.unlock();                             // 之后只读链表，不再访问 sample_vec

            // 2) 保证 hint 在目标区间左侧
            while (!hint->isHead() && hint->start > range_start)
                hint = hint->primary_path.prev.load(std::memory_order_acquire);

            SegPtr prev = hint;
            SegPtr curr = hint->primary_path.next.load(std::memory_order_acquire);

            // 3) 向右遍历，直到越过 range_start
            while (curr && !curr->isTail() && curr->start < range_start) {
                prev = curr;
                curr = curr->primary_path.next.load(std::memory_order_acquire);
            }
            return { prev, curr ? curr : tail };
    }

    bool GenomeEnd::spliceRange(SegPtr prev, SegPtr next,
        SegPtr first, SegPtr last)
    {
        // Step‑0: stitch inside chain
        // 使用 release 语义，保证链表完整性对其他线程可见
        first->primary_path.prev.store(prev, std::memory_order_release);
        last->primary_path.next.store(next, std::memory_order_release);
        std::atomic_thread_fence(std::memory_order_release);

        // Step‑1: CAS prev‑>next
        SegPtr expected = next;
        if (!prev->primary_path.next.compare_exchange_strong(
            expected, first,
            std::memory_order_acq_rel,
            std::memory_order_acquire)) {
            return false; // caller must retry
        }

        // Step‑2: fix next‑>prev
        next->primary_path.prev.store(last, std::memory_order_release);
        return true;
    }

    void GenomeEnd::spliceSegmentChain(const std::vector<SegPtr>& segs,
        uint_t beg, uint_t end) {
        if (segs.empty()) return;
        for (;;) {
            auto [prev, next] = findSurrounding(beg, end);
            if (spliceRange(prev, next, segs.front(), segs.back()))
                break;        // 成功才退出
        }
        updateSampling(segs); // 插入成功后修补采样表
    }

    void GenomeEnd::clearAllSegments() {
        std::unique_lock lk(rw);

        // 1) 复位双向链表
        head->primary_path.next.store(tail, std::memory_order_relaxed);
        tail->primary_path.prev.store(head, std::memory_order_release);

        // 2) 清空采样表，只保留 slot0 指向 head
        sample_vec.clear();
        sample_vec.resize(1, head);
    }
    
    bool GenomeEnd::removeSegment(SegPtr segment) {
        if (!segment || segment->isHead() || segment->isTail()) return false;
        
        std::unique_lock lk(rw);
        
        // 验证segment确实在这个链表中
        SegPtr current = head->primary_path.next.load(std::memory_order_acquire);
        bool found = false;
        while (current && !current->isTail()) {
            if (current == segment) {
                found = true;
                break;
            }
            current = current->primary_path.next.load(std::memory_order_acquire);
        }
        
        if (!found) return false;
        
        // 从链表中移除
        Segment::unlinkSegment(segment);
        
        // 更新采样表
        invalidateSampling(segment->start, segment->start + segment->length);
        
        return true;
    }
    
    bool GenomeEnd::removeSegmentRange(uint_t range_start, uint_t range_end) {
        std::unique_lock lk(rw);
        
        std::vector<SegPtr> to_remove;
        SegPtr current = head->primary_path.next.load(std::memory_order_acquire);
        
        // 收集需要删除的segments
        while (current && !current->isTail()) {
            if (current->start >= range_start && 
                current->start + current->length <= range_end) {
                to_remove.emplace_back(current);
            }
            current = current->primary_path.next.load(std::memory_order_acquire);
        }
        
        // 批量删除
        for (SegPtr seg : to_remove) {
            Segment::unlinkSegment(seg);
        }
        
        if (!to_remove.empty()) {
            invalidateSampling(range_start, range_end);
        }
        
        return !to_remove.empty();
    }
    
    void GenomeEnd::removeBatch(const std::vector<SegPtr>& segments) {
        if (segments.empty()) return;
        
        std::unique_lock lk(rw);
        
        uint_t min_pos = UINT32_MAX;
        uint_t max_pos = 0;
        
        // 批量删除segments并记录范围
        for (SegPtr seg : segments) {
            if (seg && !seg->isHead() && !seg->isTail()) {
                Segment::unlinkSegment(seg);
                min_pos = std::min(min_pos, seg->start);
                max_pos = std::max(max_pos, seg->start + seg->length);
            }
        }
        
        // 更新采样表
        if (min_pos != UINT32_MAX) {
            invalidateSampling(min_pos, max_pos);
        }
    }
    
    void GenomeEnd::invalidateSampling(uint_t start, uint_t end) {
        // 计算受影响的采样区间
        std::size_t start_idx = start / kSampleStep;
        std::size_t end_idx = end / kSampleStep + 1;
        
        // 重建受影响区间的采样
        for (std::size_t idx = start_idx; idx <= end_idx && idx < sample_vec.size(); ++idx) {
            sample_vec[idx] = nullptr;
            
            // 重新扫描该区间找到最左的segment
            uint_t region_start = idx * kSampleStep;
            uint_t region_end = (idx + 1) * kSampleStep;
            
            SegPtr current = head->primary_path.next.load(std::memory_order_acquire);
            while (current && !current->isTail()) {
                if (current->start >= region_start && current->start < region_end) {
                    if (!sample_vec[idx] || current->start < sample_vec[idx]->start) {
                        sample_vec[idx] = current;
                    }
                }
                current = current->primary_path.next.load(std::memory_order_acquire);
            }
        }
    }

    /* -------------------------------------------------------------
     *  移除同一染色体链表中的重叠 Segment
     *  规则：两条 Segment 区间有任何交叠时，删除 length 较小者
     * -------------------------------------------------------------*/
    void RaMesh::GenomeEnd::removeOverlap()
    {
        std::unique_lock<std::shared_mutex> lk(rw);      // 独占写

        uint_t min_pos = UINT32_MAX;
        uint_t max_pos = 0;

        if (!head || !tail) return;

        SegPtr prev = head->primary_path.next.load(std::memory_order_acquire);
        if (!prev || prev == tail) return;

        SegPtr cur = prev->primary_path.next.load(std::memory_order_acquire);

        while (cur && cur != tail)
        {
            uint_t prev_end = prev->start + prev->length;
            uint_t cur_start = cur->start;

            /* -------- 检测交叠 -------- */
            if (prev_end > cur_start)
            {
                /* 选出要删除的较短段 */
                SegPtr victim = (prev->length <= cur->length) ? prev : cur;
                SegPtr keeper = (victim == prev) ? cur : prev;

                uint_t v_beg = victim->start;
                uint_t v_end = victim->start + victim->length;

                min_pos = std::min(min_pos, v_beg);
                max_pos = std::max(max_pos, v_end);
                if(victim->parent_block)
				    victim->parent_block->removeAllSegments(); // 从 block 中移除

                cur = keeper->primary_path.next.load(std::memory_order_acquire);
				prev = keeper;                // 继续比较 keeper 与新 cur
                continue;                           // 重新比较 keeper 与新 cur
            }

            /* -------- 无交叠，正常前进 -------- */
            prev = cur;
            cur = cur->primary_path.next.load(std::memory_order_acquire);
        }

        // 更新采样表
        if (min_pos != UINT32_MAX) {
            invalidateSampling(min_pos, max_pos);
        }
    }


    /* =============================================================
     * 2. Block factories
     * ===========================================================*/
    BlockPtr Block::create(std::size_t hint)
    {
        auto bp = std::make_shared<Block>();
        bp->anchors.reserve(hint);
        return bp;
    }

    BlockPtr Block::createEmpty(const ChrName& chr, std::size_t hint)
    {
        auto bp = Block::create(hint);
        bp->ref_chr = chr;
        return bp;
    }

    std::pair<SegPtr, SegPtr> Block::createSegmentPair(const Match& match,
        const SpeciesName& ref_name,
        const SpeciesName& qry_name,
        const ChrName& ref_chr,
        const ChrName& qry_chr,
        const BlockPtr& blk)
    {
        // Create ref segment
        SegPtr ref_seg = Segment::createFromRegion(const_cast<Region&>(match.ref_region),
            Strand::FORWARD, Cigar_t{ cigarToInt('M', match.ref_region.length) }, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Create qry segment
        SegPtr qry_seg = Segment::createFromRegion(const_cast<Region&>(match.query_region),
            match.strand, Cigar_t{ cigarToInt('M', match.query_region.length) }, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Register anchors
        {
            std::unique_lock lk(blk->rw);
            blk->anchors[{ ref_name, ref_chr }] = ref_seg;
            blk->anchors[{ qry_name, qry_chr }] = qry_seg;
        }
        return { ref_seg, qry_seg };
    }

    std::pair<SegPtr, SegPtr> Block::createSegmentPair(const Anchor& anchor,
        const SpeciesName& ref_name,
        const SpeciesName& qry_name,
        const ChrName& ref_chr,
        const ChrName& qry_chr,
        const BlockPtr& blk)
    {
        // Create ref segment
        SegPtr ref_seg = Segment::createFromRegion(const_cast<Region&>(anchor.match.ref_region),
            Strand::FORWARD, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Create qry segment
        SegPtr qry_seg = Segment::createFromRegion(const_cast<Region&>(anchor.match.query_region),
            anchor.match.strand, anchor.cigar, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Register anchors
        {
            std::unique_lock lk(blk->rw);
            blk->anchors[{ ref_name, ref_chr }] = ref_seg;
            blk->anchors[{ qry_name, qry_chr }] = qry_seg;
        }
        return { ref_seg, qry_seg };
    }
    
    void Block::removeAllSegments() {
        std::unique_lock lk(rw);
        
        // 从所有链表中解除链接并清理anchors
        for (auto& [species_chr, segment] : anchors) {
            if (segment) {
                Segment::unlinkSegment(segment);
                segment->parent_block.reset();
                delete segment;
            }
        }
        anchors.clear();
    }
    
    void Block::removeSegmentsBySpecies(const SpeciesName& species) {
        std::unique_lock lk(rw);
        
        std::vector<SpeciesChrPair> to_remove;
        
        // 收集需要删除的entries
        for (const auto& [species_chr, segment] : anchors) {
            if (species_chr.first == species) {
                to_remove.emplace_back(species_chr);
                if (segment) {
                    Segment::unlinkSegment(segment);
                    segment->parent_block.reset();
                    delete segment;
                }
            }
        }
        
        // 从anchors中移除
        for (const auto& key : to_remove) {
            anchors.erase(key);
        }
    }
    
    bool Block::removeSegment(const SpeciesName& species, const ChrName& chr) {
        std::unique_lock lk(rw);
        
        SpeciesChrPair key{species, chr};
        auto it = anchors.find(key);
        
        if (it != anchors.end()) {
            if (it->second) {
                Segment::unlinkSegment(it->second);
                it->second->parent_block.reset();
                delete it->second;
            }
            anchors.erase(it);
            return true;
        }
        return false;
    }

} // namespace RaMesh