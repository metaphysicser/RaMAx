/* ramesh.cpp —— patch 2025-06-24 */
#include "ramesh.h"

namespace RaMesh {

    /* ========================================================= *
 *  Segment 带 SegmentRole 的工厂实现                          *
 *  （放在 namespace RaMesh 内，与其他实现并列）              *
 * ========================================================= */
    SegPtr Segment::create(uint_t        start,
        uint_t        len,
        Strand        sd,
        Cigar_t       cg,
        AlignRole     rl,
        SegmentRole   sl,
        const BlockPtr& bp)
    {
        auto* seg = new Segment;

        /* —— 基本字段 —— */
        seg->start = start;
        seg->length = len;
        seg->cigar = cg;
        seg->strand = sd;

        /* —— 角色标记 —— */
        seg->align_role = rl;         // PRIMARY / SECONDARY
        seg->seg_role = sl;         // SEGMENT / HEAD / TAIL

        /* —— 所属 Block —— */
        if (bp) seg->parent_block = bp;

        /* —— 链表指针初始为空 —— */
        seg->primary_path.next = nullptr;
        seg->primary_path.prev = nullptr;
        seg->secondary_path.next = nullptr;
        seg->secondary_path.prev = nullptr;

        return seg;                   // SegPtr = Segment*
    }

    SegPtr Segment::create_from_region(Region& region,
        Strand      sd,
        Cigar_t     cg,
        AlignRole   rl,
        SegmentRole sl,
        const BlockPtr& bp)
    {
        /* 统一转调到上一个工厂，防止代码重复 */
        return Segment::create(region.start,
            region.length,
            sd,
            cg,
            rl,
            sl,
            bp);
    }


    /******************  Sentinel factories  ******************/
    SegPtr Segment::create_head()
    {
        auto* seg = new Segment;
        seg->seg_role = SegmentRole::HEAD;
        seg->align_role = AlignRole::PRIMARY;   // 无意义，占位
        // 链表指针置空
        seg->primary_path.next = nullptr;
        seg->primary_path.prev = nullptr;
        return seg;
    }

    SegPtr Segment::create_tail()
    {
        auto* seg = new Segment;
        seg->seg_role = SegmentRole::TAIL;
        seg->align_role = AlignRole::PRIMARY;
        seg->primary_path.next = nullptr;
        seg->primary_path.prev = nullptr;
        return seg;
    }

    ///******************  GenomeEnd ctor  **********************/
    //GenomeEnd::GenomeEnd()
    //{
    //    SegPtr headSeg = Segment::create_head();
    //    SegPtr tailSeg = Segment::create_tail();

    //    head = make_atom(headSeg);
    //    tail = make_atom(tailSeg);

    //    headSeg->primary_path.next = tail;
    //    tailSeg->primary_path.prev = head;
    //}


    /* ========================================================= *
     * 4. 区间定位  prev / next                                   *
     * ========================================================= */
    std::pair<SegAtomPtr, SegAtomPtr>
        GenomeEnd::findSurroundingSegmentRange(uint_t range_start,
            uint_t range_end)
    {
        SegAtomPtr headAtom = head;
        SegAtomPtr tailAtom = tail;
        SegPtr     tailPtr = tail->load();


        bool head_next_tail = headAtom->load()->primary_path.next->load()->is_tail();
        /* 空链：head 直接指向 tail */
        if (head_next_tail) {
            return { headAtom, tailAtom };
        }


        SegAtomPtr prevAtom = headAtom;
        SegAtomPtr currAtom =(headAtom->load())->primary_path.next;

        while (currAtom && currAtom->load() && currAtom->load() != tailPtr)
        {
			SegPtr seg = currAtom->load();
            uint_t s_beg = seg->start;
            uint_t s_end = seg->start + seg->length;

            if (s_end <= range_start) {                   // 仍在区间左侧
                prevAtom = currAtom;
                currAtom = seg->primary_path.next ? seg->primary_path.next
                    : tailAtom;
            }
            else if (s_beg >= range_end) {                // 已到右侧
                break;
            }
            else {                                        // overlap
                break;
            }
        }
        return { prevAtom, currAtom };
    }

    /* ========================================================= *
     * 5. Block 工厂                                              *
     * ========================================================= */
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
        Strand   sd,
        AlignRole rl)
    {
        auto bp = Block::make(1);
        auto* s = Segment::create_from_region(
            const_cast<Region&>(region), sd, Cigar_t{}, rl, SegmentRole::SEGMENT, bp);

        bp->anchors[{ "", region.chr_name }] = make_atom(s);
        bp->ref_chr = region.chr_name;
        return bp;
    }

    BlockPtr Block::create_from_match(const Match& match)
    {
        auto bp = Block::make(2);

        Segment* ref = Segment::create_from_region(
            const_cast<Region&>(match.ref_region),
            match.strand,
            Cigar_t{}, AlignRole::PRIMARY, SegmentRole::SEGMENT, bp);

        Segment* qry = Segment::create_from_region(
            const_cast<Region&>(match.query_region),
            match.strand,
            Cigar_t{}, AlignRole::PRIMARY, SegmentRole::SEGMENT, bp);

        bp->anchors[{ "", match.ref_region.chr_name  }] = make_atom(ref);
        bp->anchors[{ "", match.query_region.chr_name}] = make_atom(qry);

        bp->ref_chr = match.ref_region.chr_name;
        return bp;
    }

} // namespace RaMesh
