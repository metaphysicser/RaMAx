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

} // namespace RaMesh