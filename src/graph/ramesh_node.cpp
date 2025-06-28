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
            segs[i + 1]->primary_path.prev.store(segs[i], std::memory_order_relaxed);
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
        tail->primary_path.prev.store(head, std::memory_order_relaxed);

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
        first->primary_path.prev.store(prev, std::memory_order_relaxed);
        last->primary_path.next.store(next, std::memory_order_relaxed);
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

    /* =============================================================
     * 2. Block factories
     * ===========================================================*/
    BlockPtr Block::make(std::size_t hint)
    {
        auto bp = std::make_shared<Block>();
        bp->anchors.reserve(hint);
        return bp;
    }

    BlockPtr Block::createEmpty(const ChrName& chr, std::size_t hint)
    {
        auto bp = Block::make(hint);
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

} // namespace RaMesh