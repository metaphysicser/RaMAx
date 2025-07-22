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
        return std::shared_ptr<Segment>(s);
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
        return std::shared_ptr<Segment>(h);
    }

    SegPtr Segment::createTail()
    {
        auto* t = new Segment();
        t->seg_role = SegmentRole::TAIL;
        t->align_role = AlignRole::PRIMARY;
        t->primary_path.next.store(nullptr, std::memory_order_relaxed);
        t->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return std::shared_ptr<Segment>(t);
    }

    //void Segment::linkChain(const std::vector<SegPtr>& segs)
    //{
    //    for (size_t i = 0; i + 1 < segs.size(); ++i) {
    //        segs[i]->primary_path.next.store(segs[i + 1], std::memory_order_relaxed);
    //        segs[i + 1]->primary_path.prev.store(segs[i], std::memory_order_release);
    //    }
    //}
    void Segment::linkChain(const std::vector<SegPtr>& segs)
    {
        for (size_t i = 0; i + 1 < segs.size(); ++i) {
            segs[i]->primary_path.next.store(segs[i + 1], std::memory_order_release);
            segs[i + 1]->primary_path.prev.store(segs[i], std::memory_order_release);
        }
    }

    
    //void Segment::unlinkSegment(SegPtr segment) {
    //    if (!segment || segment->isHead() || segment->isTail()) return;
    //    
    //    // 获取前驱和后继
    //    SegPtr prev = segment->primary_path.prev.load(std::memory_order_acquire);
    //    SegPtr next = segment->primary_path.next.load(std::memory_order_acquire);
    //    
    //    if (prev && next) {
    //        // 原子地更新链表指针
    //        prev->primary_path.next.store(next, std::memory_order_release);
    //        next->primary_path.prev.store(prev, std::memory_order_release);
    //        
    //        // 清空被删除segment的指针
    //        segment->primary_path.next.store(nullptr, std::memory_order_relaxed);
    //        segment->primary_path.prev.store(nullptr, std::memory_order_relaxed);
    //    }
    //}
// 把  | prev ⇆ seg ⇆ next |  改成  | prev ⇆ next |
    void Segment::unlinkSegment(SegPtr seg)
    {
        if (!seg || seg->isHead() || seg->isTail()) return;

        while (true) {
            SegPtr prev = seg->primary_path.prev.load(std::memory_order_acquire);
            SegPtr next = seg->primary_path.next.load(std::memory_order_acquire);

            // 已经被摘过
            if (!prev || !next) return;

            /* -- CAS  prev->next  -- */
            SegPtr exp = seg;
            if (!prev->primary_path.next.compare_exchange_weak(
                exp, next,
                std::memory_order_acq_rel,
                std::memory_order_acquire))
                continue;                       // 失败：再读一遍重试

            /* -- CAS  next->prev  -- */
            exp = seg;
            while (!next->primary_path.prev.compare_exchange_weak(
                exp, prev,
                std::memory_order_acq_rel,
                std::memory_order_acquire))
            {
                if (exp != seg) break;          // 别人已修好
            }

            //seg->primary_path.next.store(nullptr, std::memory_order_release);
            //seg->primary_path.prev.store(nullptr, std::memory_order_release);
            return;                             // 至此链表闭合
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
    }

    /* =============================================================
     * 1. GenomeEnd helpers
     * ===========================================================*/
    GenomeEnd::GenomeEnd() {
        head = Segment::createHead();
        tail = Segment::createTail();
        head->primary_path.next.store(tail, std::memory_order_relaxed);
        tail->primary_path.prev.store(head, std::memory_order_release);

        sample_vec.resize(1, head);   // slot 0 永远指向 head
    }

    /* ---------- 采样表维护 ---------- */
    void GenomeEnd::ensureSampleSize(uint_t pos) {
        std::size_t need = pos / kSampleStep + 1;
        if (need > sample_vec.size()) sample_vec.resize(need, nullptr);
    }

    void GenomeEnd::setToSampling(SegPtr cur) {
        // std::unique_lock lk(rw);
        std::size_t idx = cur->start / kSampleStep;
        
        if (idx == 0) {
            return;
        }
        if (idx + 1 > sample_vec.size()) sample_vec.resize(idx + 1, nullptr);
        //if (!sample_vec[need] || !sample_vec[need]->parent_block) {
        //    sample_vec[need] = cur; 
        //}else if (cur->start > sample_vec[need]->start) {
        //    sample_vec[need] = cur;
        //}
        if (!sample_vec[idx] || !sample_vec[idx]->parent_block || cur->start > sample_vec[idx]->start)
        {
            sample_vec[idx] = cur;

            SegPtr cur2 = head;
            bool t = true;
            while (cur2 != tail) {
                if (sample_vec[idx] == cur2) {
                    t = false;
                    break;
                }
                if (!cur2->primary_path.next.load()) {
                    std::cout << "";
                }
                cur2 = cur2->primary_path.next.load();
            }
            if (t) {
                std::cout << "";
            }
            //auto ptr = sample_vec[idx];
            //if (ptr) {
            //    std::cout
            //        << "sample_vec[" << idx << "] @ " << ptr.get()
            //        << "  start=" << ptr->start
            //        << "  length=" << ptr->length
            //        << "  strand=" << (ptr->strand == Strand::FORWARD ? '+' : '-')
            //        << '\n';
            //}
            //else {
            //    std::cout << "sample_vec[" << idx << "] is null\n";
            //}
        }
		
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

    SegPtr GenomeEnd::findSurrounding(uint_t range_start) {
        // 1) 读取采样表得到“最近前驱”的 hint
        //std::shared_lock lk(rw);                 // 读锁即可
        std::size_t slot = std::max((size_t)(range_start / kSampleStep) - 1, (size_t)0);
        //std::size_t slot = 1;
        //SegPtr hint = (slot < sample_vec.size() && sample_vec[slot])
        //? sample_vec[slot] : head;
        //SegPtr hint = head;
        //lk.unlock();                             // 之后只读链表，不再访问 sample_vec
        SegPtr hint = (slot < sample_vec.size() && sample_vec[slot])
            ? sample_vec[slot] : head;

        bool head_could_go = false;
        SegPtr a = head;
        while (a != tail) {
            if (a == hint) {
				head_could_go = true; // 说明 hint 在链表中
				break;
            }
			a = a->primary_path.next.load(std::memory_order_acquire);
        }
		if (!head_could_go) {
            std::cout << "";
		}
        SegPtr hint2 = head;
        // 2) 保证 hint 在目标区间左侧
        while (!hint->isHead() && hint->start > range_start) {
            SegPtr nxt = hint->primary_path.prev.load(std::memory_order_acquire);
            if (!nxt) {                    // 保险：有人在并发删除
                hint = head;               // 回到链表起点重新来
                break;
            }
            hint = nxt;
        }

        while (!hint2->isHead() && hint2->start > range_start) {
            SegPtr nxt = hint2->primary_path.prev.load(std::memory_order_acquire);
            if (!nxt) {                    // 保险：有人在并发删除
                hint2 = head;               // 回到链表起点重新来
                break;
            }
            hint2 = nxt;
        }


        SegPtr prev = hint;
        SegPtr curr = hint->primary_path.next.load(std::memory_order_acquire);

        // 3) 向右遍历，直到越过 range_start
        while (curr && !curr->isTail() && curr->start <= range_start) {
            prev = curr;
            curr = curr->primary_path.next.load(std::memory_order_acquire);
        }

		if (curr->primary_path.prev.load(std::memory_order_acquire) != prev) {
			std::cout << "";
		}

        SegPtr prev2 = hint2;
        SegPtr curr2 = hint2->primary_path.next.load(std::memory_order_acquire);

        bool t = true;
        // 3) 向右遍历，直到越过 range_start
        while (curr2 && !curr2->isTail() && curr2->start <= range_start) {
            prev2 = curr2;
            curr2 = curr2->primary_path.next.load(std::memory_order_acquire);
            if (prev2 == hint) {
                t == false;
            }
        }
        if (t) {
            std::cout << "";
        }

        if (prev2 != prev) {
            std::cout << "";
        }

        return prev2;
    }




    void GenomeEnd::insertSegment(const SegPtr seg)
    {
        if (!seg) return;

        uint_t beg = seg->start;
        if (beg == 1461897 || beg == 4014548) {
            std::cout << "";
        };
        // 1) 找到目标区间的前驱/后继（只读操作）
        SegPtr prev = findSurrounding(beg);
        
        SegPtr next = prev->primary_path.next.load();
        if (beg == 612415 || prev->start == 612415 || next->start == 612415) {
            std::cout << "";
        }

        if (prev->start == 1420786 || next->start == 1420786 || seg->start == 1420786) {
            std::cout << "";
        }

        if (prev->start == 1286480 || next->start == 1286480 || seg->start == 1286480) {
            std::cout << "";
        }

        if (prev->start == 1429536 || next->start == 1429536 || seg->start == 1429536) {
            std::cout << "";
        }
        if (beg < prev->start || (beg > next->start && next != tail)) {
            std::cout << "";
        }
        if (prev->start == 1389740 && next->start == 1420786) {
            std::cout << "";
        }

        if (prev->start == 1286480 && next->start == 1420786) {
            std::cout << "";
        }
		if (next->primary_path.prev.load(std::memory_order_acquire) != prev) {
            std::cout << "";
;		}


        seg->primary_path.prev.store(prev);
        if (seg->start == 1286480 && next->start == 1429536) {
            std::cout << "";
        }
        seg->primary_path.next.store(next);

        if (prev->start == 1286480 && seg->start == 1429536) {
            std::cout << "";
        }
        prev->primary_path.next.store(seg);
        next->primary_path.prev.store(seg);

        // 检验链表是否正确，检验seg，prev，next双向是否正确
		if (next->primary_path.prev.load(std::memory_order_acquire) != seg) {
			std::cout << "";
		}
		if (prev->primary_path.next.load(std::memory_order_acquire) != seg) {
			std::cout << "";
		}
		if (seg->primary_path.prev.load(std::memory_order_acquire) != prev) {
			std::cout << "";
		}
		if (seg->primary_path.next.load(std::memory_order_acquire) != next) {
			std::cout << "";
		}


        setToSampling(seg);
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
    void RaMesh::GenomeEnd::removeOverlap(bool if_ref)
    {
        std::unique_lock<std::shared_mutex> lk(rw);      // 独占写


        if (!head || !tail) return;

        SegPtr prev = head->primary_path.next.load(std::memory_order_acquire);
        if (!prev || prev == tail) return;

        SegPtr cur = prev->primary_path.next.load(std::memory_order_acquire);

        while (cur && cur != tail)
        {
            SegPtr next = cur;
            uint_t prev_end = prev->start + prev->length;
            uint_t cur_start = cur->start;
            if (cur_start == 583001 || cur_start == 2956412 || cur_start == 4488638) {
                std::cout << "";
            }
            /* -------- 检测交叠 -------- */
            while (prev_end > cur_start && cur && cur != tail && next && next != tail)
            {
                cur_start = cur->start;

                if (if_ref && prev->cigar.size() == 0 && cur->cigar.size() == 0) {
                    cur = cur->primary_path.next.load(std::memory_order_acquire);
                }
                else {
                    bool cur_longer = prev->length <= cur->length;
                    if (cur_longer) {
                        prev->parent_block->removeAllSegments();
                        break;
                    }
                    else {
                        if (next == cur) {
                            next = cur->primary_path.next.load(std::memory_order_acquire);
                        }
                        cur = cur->primary_path.next.load(std::memory_order_acquire);                   
                        SegPtr cur_prev = cur->primary_path.prev.load(std::memory_order_acquire);                       
                        cur_prev->parent_block->removeAllSegments();
                    }
                }
            }
            // setToSampling(next);
            /* -------- 无交叠，正常前进 -------- */
            prev = next;
            cur = prev->primary_path.next.load(std::memory_order_acquire);
        }

    }
    /* -------------------------------------------------------------
 *  移除同一染色体链表中的重叠 Segment
 *  规则：两条 Segment 区间有任何交叠时，删除 length 较小者
 * ------------------------------------------------------------*/
    //void GenomeEnd::removeOverlap(bool if_ref)
    //{
    //    std::unique_lock lk(rw);                 // 串行调用，独占即可

    //    SegPtr seg = head->primary_path.next.load(std::memory_order_relaxed);
    //    while (seg && !seg->isTail())
    //    {
    //        SegPtr nxt = seg->primary_path.next.load(std::memory_order_relaxed);

    //        /* 检查 seg 与后继是否重叠 —— 只要起点落在 seg 区间内就算重叠 */
    //        while (nxt && !nxt->isTail() &&
    //            nxt->start < seg->start + seg->length)
    //        {
    //            /* 如果是 ref 链且两段都是“纯 M”（cigar 为空），
    //               你曾选择保留两段，这里沿用原逻辑                 */
    //            if (if_ref && seg->cigar.empty() && nxt->cigar.empty())
    //            {
    //                break;                       // 不删除任何一条
    //            }

    //            /* 选出较短者作为待删对象 */
    //            SegPtr victim = (seg->length <= nxt->length) ? seg : nxt;
    //            SegPtr survivor = (victim == seg) ? nxt : seg;

    //            /* ---- 真正删除 ---- */
    //            Segment::deleteSegment(victim);  // 已包含 unlink + 从 block 的 anchors 擦除

    //            /* 维护采样表 */
    //            //invalidateSampling(victim->start,
    //                //victim->start + victim->length);

    //            /* 继续比较 survivor 与下一个 */
    //            seg = survivor;
    //            nxt = seg->primary_path.next.load(std::memory_order_relaxed);
    //        }

    //        /* 无重叠，正常前进 */
    //        seg = nxt;
    //    }
    //}




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
    
    //void Block::removeAllSegments() {
    //    std::unique_lock lk(rw);
    //    
    //    // 从所有链表中解除链接并清理anchors
    //    for (auto& [species_chr, segment] : anchors) {
    //        if (segment) {
    //            Segment::unlinkSegment(segment);
    //            segment->parent_block.reset();
    //            delete segment;
    //        }
    //    }
    //    anchors.clear();
    //}
    void Block::removeAllSegments()
    {
        std::unique_lock lk(rw);                     // 独占 Block

        std::vector<SegPtr> seg_list;
        seg_list.reserve(anchors.size());

        for (auto& [_, seg] : anchors)       // 只读，绝不修改 anchors
            if (seg) seg_list.push_back(seg);

        anchors.clear();                              // 现在可以一次性清空

        /* 真正断链 + 释放在容器之外进行，
           即使过程中触发再修改 anchors，也不会再有活迭代器。 */
        for (SegPtr seg : seg_list) {
            Segment::unlinkSegment(seg);
            seg->parent_block.reset();
        }
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
            }
            anchors.erase(it);
            return true;
        }
        return false;
    }

} // namespace RaMesh