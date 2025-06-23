/* ramesh.cpp  (与 ramesh.h v0.4 对应) */
#include "ramesh.h"

namespace RaMesh {

    /* ------------------------------------------------------------------ */
        /* SegPtr 工厂                                                         */
        /* ------------------------------------------------------------------ */
    SegPtr Segment::create(uint_t        start,
        uint_t        length,
        const Cigar_t& cigar,
        Strand        strand,
        AlignRole     role,
        const BlockPtr& parent)
    {
        auto* seg = new Segment;
        seg->start = start;
        seg->length = length;
        seg->cigar = cigar;
        seg->strand = strand;
        seg->role = role;
        if (parent) seg->parent_block = parent;

        seg->primary_path.next = nullptr;
        seg->primary_path.prev = nullptr;
        seg->secondary_path.next = nullptr;
        seg->secondary_path.prev = nullptr;
        return seg;
    }

    SegPtr Segment::create(uint_t        start,
        uint_t        length,
        Strand        strand,
        Cigar_t       cigar,
        AlignRole     role,
        const BlockPtr& parent)
    {
        return Segment::create(start, length, cigar, strand, role, parent);
    }

    SegPtr Segment::create_from_region(Region& region,
        Strand      strand,
        Cigar_t     cigar,
        AlignRole   role,
        const BlockPtr& parent)
    {
        return Segment::create(region.start,
            region.length,
            cigar,
            strand,
            role,
            parent);
    }

    /* ------------------------------------------------------------------ */
    /* Sentinel 节点                                                       */
    /* ------------------------------------------------------------------ */
    std::shared_ptr<HeadSegment> HeadSegment::create()
    {
        auto h = std::make_shared<HeadSegment>();
        h->primary_path.next = nullptr;
        h->primary_path.prev = nullptr;
        return h;
    }

    std::shared_ptr<TailSegment> TailSegment::create()
    {
        auto t = std::make_shared<TailSegment>();
        t->primary_path.next = nullptr;
        t->primary_path.prev = nullptr;
        return t;
    }



    /* ------------------------------------------------------------------ */
    /* GenomeEnd 构造 & 查找                                               */
    /* ------------------------------------------------------------------ */
    GenomeEnd::GenomeEnd()
        : head(HeadSegment::create()),
        tail(TailSegment::create())
    {
        head->primary_path.next = std::make_shared<SegAtom>(
            reinterpret_cast<SegPtr>(tail.get()));
        tail->primary_path.prev = std::make_shared<SegAtom>(
            reinterpret_cast<SegPtr>(head.get()));
    }

    GenomeEnd::GenomeEnd(const std::shared_ptr<HeadSegment>& h,
        const std::shared_ptr<TailSegment>& t)
        : head(h), tail(t)
    {
        head->primary_path.next = std::make_shared<SegAtom>(
            reinterpret_cast<SegPtr>(t.get()));
        t->primary_path.prev = std::make_shared<SegAtom>(
            reinterpret_cast<SegPtr>(h.get()));
    }

    std::pair<SegPtr, SegPtr>
        GenomeEnd::findSurroundingSegmentRange(uint_t range_start,
            uint_t range_end)
    {
        SegPtr prev = reinterpret_cast<SegPtr>(head.get());
        SegPtr curr = head->primary_path.next
            ? head->primary_path.next->load()
            : reinterpret_cast<SegPtr>(tail.get());

        while (curr && curr != reinterpret_cast<SegPtr>(tail.get()))
        {
            auto* seg = dynamic_cast<Segment*>(curr);
            if (!seg) {
                curr = reinterpret_cast<RaMeshNode*>(curr)
                    ->primary_path.next
                    ? reinterpret_cast<RaMeshNode*>(curr)
                    ->primary_path.next->load()
                    : reinterpret_cast<SegPtr>(tail.get());
                continue;
            }

            uint_t seg_start = seg->start;
            uint_t seg_end = seg->start + seg->length;

            if (seg_end <= range_start) {
                prev = curr;
                curr = seg->primary_path.next
                    ? seg->primary_path.next->load()
                    : reinterpret_cast<SegPtr>(tail.get());
            }
            else if (seg_start >= range_end) {
                break;
            }
            else {                            // overlap
                break;
            }
        }
        return { prev, curr };
    }

    /* ------------------------------------------------------------------ */
    /* Block 工厂                                                         */
    /* ------------------------------------------------------------------ */
    BlockPtr Block::make(std::size_t genome_hint)
    {
        auto bp = std::make_shared<Block>();
        bp->anchors.reserve(genome_hint);
        return bp;
    }

    BlockPtr Block::create_empty(const ChrName& ref_chr,
        std::size_t    genome_hint)
    {
        auto bp = Block::make(genome_hint);
        bp->ref_chr = ref_chr;
        return bp;
    }

    BlockPtr Block::create_from_region(const Region& region,
        Strand        strand,
        AlignRole     role)
    {
        auto bp = Block::make(1);

        Segment* seg = Segment::create_from_region(
            const_cast<Region&>(region), strand, Cigar_t{}, role, bp);

        SpeciesChrPair key{ "", region.chr_name };
        bp->anchors[key] = std::make_shared<SegAtom>(seg);
        bp->ref_chr = region.chr_name;
        return bp;
    }

    BlockPtr Block::create_from_match(const Match& match)
    {
        auto bp = Block::make(2);

        Segment* ref = Segment::create_from_region(
            const_cast<Region&>(match.ref_region),
            match.strand,
            Cigar_t{},
            AlignRole::PRIMARY,
            bp);

        Segment* qry = Segment::create_from_region(
            const_cast<Region&>(match.query_region),
            match.strand,
            Cigar_t{},
            AlignRole::SECONDARY,
            bp);

        bp->anchors[{ "", match.ref_region.chr_name  }] =
            std::make_shared<SegAtom>(ref);
        bp->anchors[{ "", match.query_region.chr_name}] =
            std::make_shared<SegAtom>(qry);

        bp->ref_chr = match.ref_region.chr_name;
        return bp;
    }


} // namespace RaMesh
