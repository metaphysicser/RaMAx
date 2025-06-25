// ******************************************************
//  ramesh.cpp – 适配 ramesh.h (v0.5-beta)
//  ✧  无锁并行插入 + Sentinel + Block 工厂实现  ✧
//  2025-06-24
// ******************************************************
#include "ramesh.h"
#include <iomanip>
#include <shared_mutex>

namespace RaMesh {

    /* ======================================================
     * 0.  Segment 工厂实现
     * ====================================================*/
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

        /* 指针初始化为空 */
        s->primary_path.next.store(nullptr, std::memory_order_relaxed);
        s->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return s;
    }

    SegPtr Segment::create_from_region(Region& region, Strand sd,
        Cigar_t cg, AlignRole rl, SegmentRole sl,
        const BlockPtr& bp)
    {
        return create(region.start, region.length, sd, std::move(cg), rl, sl, bp);
    }

    SegPtr Segment::create_head()
    {
        auto* h = new Segment();
        h->seg_role = SegmentRole::HEAD;
        h->align_role = AlignRole::PRIMARY;
        h->primary_path.next.store(nullptr, std::memory_order_relaxed);
        h->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return h;
    }

    SegPtr Segment::create_tail()
    {
        auto* t = new Segment();
        t->seg_role = SegmentRole::TAIL;
        t->align_role = AlignRole::PRIMARY;
        t->primary_path.next.store(nullptr, std::memory_order_relaxed);
        t->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return t;
    }

    /* ======================================================
     * 1. GenomeEnd helper – find_surrounding
     * ====================================================*/
    std::pair<SegPtr, SegPtr> GenomeEnd::find_surrounding(uint_t range_start,
        uint_t range_end)
    {
        SegPtr headPtr = head;
        SegPtr tailPtr = tail;

        SegPtr curr = headPtr->primary_path.next.load(std::memory_order_acquire);
        SegPtr prev = headPtr;

        /* 空链表：head 直接指向 tail */
        if (curr == tailPtr) return { headPtr, tailPtr };

        while (curr && curr != tailPtr)
        {
            uint_t seg_beg = curr->start;
            uint_t seg_end = curr->start + curr->length;

            if (seg_end <= range_start) {
                // 仍在左边，前进
                prev = curr;
                curr = curr->primary_path.next.load(std::memory_order_acquire);

            }
            else if (seg_beg >= range_end) {
                // 已越过区间右端
                break;
            }
            else {
                // overlap
                break;
            }
        }


        return { prev, curr };
    }

    /* ======================================================
     * 2. Block 工厂
     * ====================================================*/
    BlockPtr Block::make(std::size_t hint)
    {
        auto bp = std::make_shared<Block>();
        bp->anchors.reserve(hint);
        return bp;
    }

    BlockPtr Block::create_empty(const ChrName& chr, std::size_t hint)
    {
        auto bp = Block::make(hint);
        bp->ref_chr = chr;
        return bp;
    }

    BlockPtr Block::create_from_region(const Region& region,
        Strand sd, AlignRole rl)
    {
        auto bp = Block::make(1);
        SegPtr s = Segment::create_from_region(const_cast<Region&>(region), sd,
            Cigar_t{}, rl, SegmentRole::SEGMENT, bp);
        bp->anchors[{"", region.chr_name}] = s;
        bp->ref_chr = region.chr_name;
        return bp;
    }

    BlockPtr Block::create_from_match(const Match& match)
    {
        auto bp = Block::make(2);

        SegPtr ref = Segment::create_from_region(const_cast<Region&>(match.ref_region),
            match.strand, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, bp);

        SegPtr qry = Segment::create_from_region(const_cast<Region&>(match.query_region),
            match.strand, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, bp);

        bp->anchors[{"", match.ref_region.chr_name }] = ref;
        bp->anchors[{"", match.query_region.chr_name}] = qry;
        bp->ref_chr = match.ref_region.chr_name;
        return bp;
    }

    

} // namespace RaMesh
