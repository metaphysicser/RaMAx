#include "ramesh.h"

namespace RaMesh {
    SegPtr Segment::create(uint_t          start,
        uint_t          length,
        const Cigar_t& cigar,
        Strand          strand,
        AlignRole       role,
        const BlockPtr& parent)
    {
        auto* seg = new Segment();          // 1️⃣ 分配对象

        seg->start = start;
        seg->length = length;
        seg->cigar = cigar;                // 2️⃣ 填充字段
        seg->strand = strand;
        seg->role = role;

        seg->parent_block = parent;         // 3️⃣ 传入则转弱引用

        // 4️⃣ 清空链表指针（用 relaxed 免栅栏）
        seg->primary_path.next.store(nullptr, std::memory_order_relaxed);
        seg->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        seg->secondary_path.next.store(nullptr, std::memory_order_relaxed);
        seg->secondary_path.prev.store(nullptr, std::memory_order_relaxed);

        return seg;
    }

    SegPtr Segment::create(const Region& region,
        Strand strand,
        AlignRole role,
        const BlockPtr& parent)
    {
        // 在堆上 new 一个 Segment
        Segment* s = new Segment;
        // 初始化字段
        s->start = region.start;
        s->length = region.length;
        // 如果想存储 chr_name，需要在 Segment 定义中添加成员：
        //    ChrName chr_name;
        // 并在此设置：
        //    s->chr_name = region.chr_name;
        s->strand = strand;
        s->role = role;
        // cigar 初始留空，调用者可自行修改：
        s->cigar.clear();

        // parent_block 只设置弱引用，不自动插入
        if (parent) {
            s->parent_block = parent;
        }
        // RaMeshPath next/prev 默认 nullptr，无需额外操作
        return s;
    }

    // 头哨兵
    std::shared_ptr<HeadSegment> HeadSegment::create()
    {
        // 用 make_shared 一步完成分配与构造
        return std::shared_ptr<HeadSegment>(new HeadSegment());
    }

    // 尾哨兵
    std::shared_ptr<TailSegment> TailSegment::create()
    {
        return std::shared_ptr<TailSegment>(new TailSegment());
    }

    /* 已有：Block::make 实现 */
    BlockPtr Block::make(std::size_t genome_hint) {
        auto bp = BlockPtr(new Block);
        bp->anchors.reserve(genome_hint);
        return bp;
    }

    /* -------------- create_empty -------------- */
    BlockPtr Block::create_empty(const ChrName& ref_chr,
        std::size_t    genome_hint)
    {
        auto bp = Block::make(genome_hint);
        {
            std::unique_lock lk(bp->rw);
            bp->ref_chr = ref_chr;
            // anchors 不插任何 Segment，仅预留
        }
        return bp;
    }

    /* -------------- create_from_region -------------- */
    BlockPtr Block::create_from_region(const Region& region,
        Strand         strand,
        AlignRole      role)
    {
        auto bp = Block::make(1);
        {
            std::unique_lock lk(bp->rw);
            bp->ref_chr = region.chr_name;
            // 1) 在堆上 new Segment
            Segment* s = Segment::create(region, strand, role, bp);
            // 2) 插入 anchors
            auto& vec = bp->anchors[region.chr_name];
            // 注意: anchors 存储 SegAtom (atomic<Segment*>)
            vec.push_back(SegAtom(s));
        }
        return bp;
    }

    /* -------------- create_from_match -------------- */
    BlockPtr Block::create_from_match(const Match& match)
    {
        // 基于 ref_region 初始化。若需处理 query_region，请在此函数外或自行扩展
        auto bp = Block::make(2);
        {
            std::unique_lock lk(bp->rw);
            bp->ref_chr = match.ref_region.chr_name;
            // 新 Segment
            Segment* ref = Segment::create(match.ref_region,
                match.strand == Strand::FORWARD ? Strand::FORWARD : Strand::REVERSE,
                AlignRole::PRIMARY,
                bp);
            auto& vec = bp->anchors[match.ref_region.chr_name];
            vec.push_back(SegAtom(ref));

			Segment* query = Segment::create(match.query_region,
				match.strand == Strand::FORWARD ? Strand::FORWARD : Strand::REVERSE,
				AlignRole::PRIMARY,
				bp);
			// 如果需要处理 query_region，可以在这里插入
			auto& qvec = bp->anchors[match.query_region.chr_name];
			qvec.push_back(SegAtom(query));
        }
        return bp;
    }




}