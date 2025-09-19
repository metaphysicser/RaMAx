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

        SegPtr prev = seg->primary_path.prev.load(std::memory_order_acquire);
        SegPtr next = seg->primary_path.next.load(std::memory_order_acquire);

        uint_t cur_start = prev->start;

		prev->primary_path.next.store(next);
        next->primary_path.prev.store(prev);

		seg->primary_path.next.store(nullptr);
		seg->primary_path.prev.store(nullptr);
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

        if (!sample_vec[idx] || !sample_vec[idx]->parent_block || cur->start > sample_vec[idx]->start)
        {
            sample_vec[idx] = cur;
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
        //lk.unlock();                             // 之后只读链表，不再访问 sample_vec
        SegPtr hint = (slot < sample_vec.size() && sample_vec[slot])
            ? sample_vec[slot] : head;
   
        // 2) 保证 hint 在目标区间左侧
        while (!hint->isHead() && hint->start > range_start) {
            SegPtr nxt = hint->primary_path.prev.load(std::memory_order_acquire);
            if (!nxt) {                    // 保险：有人在并发删除
                hint = head;               // 回到链表起点重新来
                break;
            }
            hint = nxt;
        }

        SegPtr prev = hint;
        SegPtr curr = hint->primary_path.next.load(std::memory_order_acquire);

        // 3) 向右遍历，直到越过 range_start
        while (curr && !curr->isTail() && curr->start <= range_start) {
            prev = curr;
            curr = curr->primary_path.next.load(std::memory_order_acquire);
        }
        return prev;
    }




    void GenomeEnd::insertSegment(const SegPtr seg)
    {
        if (!seg) return;

        uint_t beg = seg->start;

        // 1) 找到目标区间的前驱/后继（只读操作）
        SegPtr prev = findSurrounding(beg);
        SegPtr next = prev->primary_path.next.load();

        seg->primary_path.prev.store(prev);
        
        seg->primary_path.next.store(next);
      
        prev->primary_path.next.store(seg);
        next->primary_path.prev.store(seg);
       
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
    void RaMesh::GenomeEnd::alignInterval(const SpeciesName ref_name,
        const SpeciesName query_name,
        const ChrName query_chr_name,
        SegPtr cur_node,
        std::map<SpeciesName, SeqPro::SharedManagerVariant> managers,
        bool left_extend,
        bool right_extend) {
        if (cur_node == head || cur_node == tail || cur_node == NULL) return;
		if (cur_node->right_extend) return; // 已经扩展过了
        auto fetchSeq = [](const SeqPro::ManagerVariant& mv,
            const ChrName& chr, Coord_t b, Coord_t l) {
                return std::visit([&](auto& p) {
                    using T = std::decay_t<decltype(p)>;
                    if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>)
                        return p->getSubSequence(chr, b, l);
                    else
                        return p->getOriginalManager().getSubSequence(chr, b, l);
                    }, mv);
            };

        auto getChrLen = [](const SeqPro::ManagerVariant& mv, const ChrName& chr) {
            return std::visit([&](auto& p) {
                using T = std::decay_t<decltype(p)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>)
                    return p->getSequenceLength(chr);
                else
                    return p->getOriginalManager().getSequenceLength(chr);
                }, mv);
            };

        // ---------------- 左扩展 ----------------
        if (left_extend && cur_node->left_extend == false) {
            
            int_t query_len = 0;
            Strand strand = cur_node->strand;
            int_t query_start = 0;

            if (strand == FORWARD) {
                // 正向：以 prev 的 end 作为左界
                SegPtr query_left_node = cur_node->primary_path.prev.load(std::memory_order_acquire);
                query_start = (!query_left_node->isHead()) ? query_left_node->start + query_left_node->length : 0;
                query_len = cur_node->start - query_start;
            }
            else {
                // 反向：以 next 的 start 作为左界
                SegPtr query_right_node = cur_node->primary_path.next.load(std::memory_order_acquire);
                query_start = cur_node->start + cur_node->length;
                if (!query_right_node->isTail())
                    query_len = query_right_node->start - query_start;
                else
                    query_len = getChrLen(*managers[query_name], query_chr_name) - query_start;
            }


            BlockPtr cur_block = cur_node->parent_block;
            ChrName ref_chr_name = cur_block->ref_chr;
            SegPtr ref_cur_node = cur_block->anchors[{ ref_name, cur_block->ref_chr }];

            cur_node->left_extend = true;
            ref_cur_node->left_extend = true;
            SegPtr ref_left_node = ref_cur_node->primary_path.prev.load(std::memory_order_acquire);
            while (!ref_left_node->isHead() && ref_left_node->parent_block->ref_chr != ref_chr_name) {
                ref_left_node = ref_left_node->primary_path.prev.load(std::memory_order_acquire);
            }

            int_t ref_start = (!ref_left_node->isHead()) ? ref_left_node->start + ref_left_node->length : 0;
            int_t ref_len = ref_cur_node->start - ref_start;

            // === 实际比对逻辑 ===
            if (query_len > 0 && ref_len > 0) {
                std::string query_seq = fetchSeq(*managers[query_name], query_chr_name, query_start, query_len);
                std::string ref_seq = fetchSeq(*managers[ref_name], cur_block->ref_chr, ref_start, ref_len);

                // TODO: 在这里执行左扩展比对，比如 Smith-Waterman 或者自定义比对逻辑
                std::reverse(ref_seq.begin(), ref_seq.end());
                if (strand == FORWARD) {
                    reverseSeq(query_seq);
                }
                else {
					baseComplement(query_seq);
                }

				Cigar_t result = extendAlignWFA2(ref_seq, query_seq);
				std::reverse(result.begin(), result.end());
				AlignCount cnt = countAlignedBases(result);
                if (strand == FORWARD) {
                    cur_node->start -= cnt.query_bases;
                }
				
				cur_node->length += cnt.query_bases;
                ref_cur_node->start -= cnt.ref_bases;
                ref_cur_node->length += cnt.ref_bases;
				prependCigar(cur_node->cigar, result);
					
            }
        }

        // ---------------- 右扩展 ----------------
        if (right_extend && cur_node->right_extend == false) {
            int_t query_len = 0;
            Strand strand = cur_node->strand;
            int_t query_start = 0;

            if (strand == FORWARD) {
                SegPtr query_right_node = cur_node->primary_path.next.load(std::memory_order_acquire);
                query_start = cur_node->start + cur_node->length;
                if (!query_right_node->isTail()) {
                    query_len = query_right_node->start - query_start;
                }
                else {
                    query_len = getChrLen(*managers[query_name], query_chr_name) - query_start;
                }
            }
            else {
                SegPtr query_left_node = cur_node->primary_path.prev.load(std::memory_order_acquire);
                if (!query_left_node->isHead()) {
                    query_start = query_left_node->start + query_left_node->length;
                    query_len = cur_node->start - query_start;
                }
                else {
                    query_start = 0;
                    query_len = cur_node->start;
                }
            }

            BlockPtr cur_block = cur_node->parent_block;
            ChrName ref_chr_name = cur_block->ref_chr;
            SegPtr ref_cur_node = cur_block->anchors[{ ref_name, cur_block->ref_chr }];

			cur_node->right_extend = true;
			ref_cur_node->right_extend = true;

            SegPtr ref_right_node = ref_cur_node->primary_path.next.load(std::memory_order_acquire);
            while (!ref_right_node->isTail() && ref_right_node->parent_block->ref_chr != ref_chr_name) {
                ref_right_node = ref_right_node->primary_path.next.load(std::memory_order_acquire);
            }          

            int_t ref_start = ref_cur_node->start + ref_cur_node->length;
            int_t ref_len = (!ref_right_node->isTail()) ? ref_right_node->start - ref_start : getChrLen(*managers[ref_name], cur_block->ref_chr) - ref_start;

            // === 实际比对逻辑 ===
            if (query_len > 0 && ref_len > 0) {
                std::string query_seq = fetchSeq(*managers[query_name], query_chr_name, query_start, query_len);
                std::string ref_seq = fetchSeq(*managers[ref_name], cur_block->ref_chr, ref_start, ref_len);

                // TODO: 在这里执行右扩展比对
                if (strand == FORWARD) {
                }
                else {
					reverseSeq(ref_seq);
                }

                Cigar_t result = extendAlignWFA2(ref_seq, query_seq);

                AlignCount cnt = countAlignedBases(result);
                if (strand == REVERSE) {
                    cur_node->start -= cnt.query_bases;
                }

                cur_node->length += cnt.query_bases;

                ref_cur_node->length += cnt.ref_bases;
                appendCigar(cur_node->cigar, result);
            }
        }
    }

    // 串行重新排序：按照 start 坐标
    void GenomeEnd::resortSegments() {
        std::vector<SegPtr> segs;
        SegPtr cur = head->primary_path.next.load(std::memory_order_acquire);
        while (cur && !cur->isTail()) {
            segs.push_back(cur);
            cur = cur->primary_path.next.load(std::memory_order_acquire);
        }

        if (segs.empty()) return;

        // 按 start 排序
        std::sort(segs.begin(), segs.end(), [](const SegPtr& a, const SegPtr& b) {
            return a->start < b->start;
            });

        // 重新链接：head -> segs[0] -> ... -> segs[n-1] -> tail
        Segment::linkChain(segs);
        head->primary_path.next.store(segs.front(), std::memory_order_release);
        segs.front()->primary_path.prev.store(head, std::memory_order_release);

        tail->primary_path.prev.store(segs.back(), std::memory_order_release);
        segs.back()->primary_path.next.store(tail, std::memory_order_release);

#ifdef _DEBUG_
        cur = head->primary_path.next.load(std::memory_order_acquire);
        uint_t prev_start = 0;
        bool first = true;
        while (cur && !cur->isTail()) {
            if (!first && cur->start < prev_start) {
                std::cerr << "[GenomeEnd::resortSegments] ERROR: "
                    << "segment order invalid: "
                    << cur->start << " < " << prev_start << std::endl;
                break; // 发现问题就退出
            }
            prev_start = cur->start;
            first = false;
            cur = cur->primary_path.next.load(std::memory_order_acquire);
        }
#endif
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
            if (cur_start == 1286480 || cur_start == 1429536) {
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
        uint_t match_len = match.match_len();
        SegPtr ref_seg = Segment::create(match.ref_start, match_len,
            Strand::FORWARD, Cigar_t{ cigarToInt('M', match_len) }, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Create qry segment
        SegPtr qry_seg = Segment::create(match.qry_start, match_len,
            match.strand(), Cigar_t{ cigarToInt('M', match_len) }, AlignRole::PRIMARY,
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
        SegPtr ref_seg = Segment::create(anchor.ref_start, anchor.ref_len,
            Strand::FORWARD, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Create qry segment
        SegPtr qry_seg = Segment::create(anchor.qry_start, anchor.qry_len,
            anchor.strand, anchor.cigar, AlignRole::PRIMARY,
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
        //std::unique_lock lk(rw);                     // 独占 Block

        std::vector<SegPtr> seg_list;
        seg_list.reserve(anchors.size());

        for (auto& [_, seg] : anchors)       // 只读，绝不修改 anchors
            if (seg) seg_list.push_back(seg);

        anchors.clear();                              // 现在可以一次性清空

        ///* 真正断链 + 释放在容器之外进行，
        //   即使过程中触发再修改 anchors，也不会再有活迭代器。 */
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