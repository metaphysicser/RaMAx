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
    GenomeEnd::GenomeEnd()
    {
        head_holder.reset(Segment::createHead());
        tail_holder.reset(Segment::createTail());
        head = head_holder.get();
        tail = tail_holder.get();

        head->primary_path.next.store(tail, std::memory_order_relaxed);
        tail->primary_path.prev.store(head, std::memory_order_relaxed);
    }

    std::pair<SegPtr, SegPtr> GenomeEnd::findSurrounding(uint_t range_start,
        uint_t range_end)
    {
        SegPtr curr = head->primary_path.next.load(std::memory_order_acquire);
        SegPtr prev = head;

        if (curr == tail) return { head, tail }; // empty list

        while (curr && curr != tail)
        {
            uint_t seg_beg = curr->start;
            uint_t seg_end = curr->start + curr->length;

            if (seg_end <= range_start) {
                prev = curr;
                curr = curr->primary_path.next.load(std::memory_order_acquire);
            }
            else if (seg_beg >= range_end) {
                break;
            }
            else {
                break; // overlap – still return [prev, curr]
            }
        }
        return { prev, curr };
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

    void GenomeEnd::spliceSegmentChain(const std::vector<SegPtr>& segments,
        uint_t beg, uint_t end)
    {
        if (segments.empty()) return;

        for (;;) {
            auto [prev, next] = findSurrounding(beg, end);
            if (spliceRange(prev, next, segments.front(), segments.back()))
                break; // success
        }
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

    BlockPtr Block::createFromRegion(const Region& region, Strand sd, AlignRole rl)
    {
        auto bp = Block::make(1);
        SegPtr s = Segment::createFromRegion(const_cast<Region&>(region), sd,
            Cigar_t{}, rl, SegmentRole::SEGMENT, bp);
        bp->anchors[{ "", region.chr_name }] = s;
        bp->ref_chr = region.chr_name;
        return bp;
    }

    BlockPtr Block::createFromMatch(const Match& match)
    {
        auto bp = Block::make(2);

        SegPtr ref = Segment::createFromRegion(const_cast<Region&>(match.ref_region),
            match.strand, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, bp);

        SegPtr qry = Segment::createFromRegion(const_cast<Region&>(match.query_region),
            match.strand, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, bp);

        bp->anchors[{ "", match.ref_region.chr_name  }] = ref;
        bp->anchors[{ "", match.query_region.chr_name}] = qry;
        bp->ref_chr = match.ref_region.chr_name;
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
            match.strand, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Create qry segment
        SegPtr qry_seg = Segment::createFromRegion(const_cast<Region&>(match.query_region),
            match.strand, Cigar_t{}, AlignRole::PRIMARY,
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