
#ifndef RAMESH_H
#define RAMESH_H

#include <atomic>
#include <memory>
#include <unordered_map>
#include <vector>
#include <map>
#include <shared_mutex>
#include <iostream>
#include "config.hpp"
#include "anchor.h"


namespace RaMesh {

    // ────────────────────────────────────────────────
    // 前向声明
    // ────────────────────────────────────────────────
    class  RaMeshNode;
    class  Block;
    struct Segment;

    using BlockPtr = std::shared_ptr<Block>;
    using WeakBlock = std::weak_ptr<Block>;
    using SegPtr = Segment*;

    // ────────────────────────────────────────────────
    // 辅助类型 / 哈希
    // ────────────────────────────────────────────────
    using SpeciesChrPair = std::pair<SpeciesName, ChrName>;


    // TODO 可以换成更快的hash
    struct SpeciesChrPairHash {
        std::size_t operator()(const SpeciesChrPair& k) const noexcept {
            std::size_t h1 = std::hash<std::string>{}(k.first);
            std::size_t h2 = std::hash<std::string>{}(k.second);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };


    using ChrHeadMap = std::unordered_map<SpeciesChrPair, SegPtr, SpeciesChrPairHash>;

    // ────────────────────────────────────────────────
    // 并发链表结构
    // ────────────────────────────────────────────────
    struct RaMeshPath {
        std::atomic<SegPtr> next{ nullptr };
        std::atomic<SegPtr> prev{ nullptr };
    };

    // ────────────────────────────────────────────────
    // 业务枚举
    // ────────────────────────────────────────────────
    enum class AlignRole : uint8_t { PRIMARY = 0, SECONDARY = 1 };
    enum class SegmentRole : uint8_t { SEGMENT = 0, HEAD = 1, TAIL = 2 };

    // ────────────────────────────────────────────────
    // Segment / Sentinel
    // ────────────────────────────────────────────────
    class Segment {
    public:
        uint_t      start{ 0 };
        uint_t      length{ 0 };
        Cigar_t     cigar;
        Strand      strand{ Strand::FORWARD };
        AlignRole   align_role{ AlignRole::PRIMARY };
        SegmentRole seg_role{ SegmentRole::SEGMENT };

        RaMeshPath  primary_path;
        // RaMeshPath  secondary_path; // 之后加入这个数据结构

        WeakBlock   parent_block;
        mutable std::shared_mutex rw; // 保护非链表字段

        bool is_head()    const noexcept { return seg_role == SegmentRole::HEAD; }
        bool is_tail()    const noexcept { return seg_role == SegmentRole::TAIL; }
        bool is_segment() const noexcept { return seg_role == SegmentRole::SEGMENT; }
        bool is_primary() const noexcept { return align_role == AlignRole::PRIMARY; }
        bool is_secondary() const noexcept { return align_role == AlignRole::SECONDARY; }

        // ——— 工厂 ———
        static SegPtr create(uint_t start, uint_t len, Strand sd,
            Cigar_t cg, AlignRole rl, SegmentRole sl,
            const BlockPtr& bp);

        static SegPtr create_from_region(Region& region, Strand sd,
            Cigar_t cg, AlignRole rl, SegmentRole sl,
            const BlockPtr& bp);

        // Sentinel 生成
        static SegPtr create_head();
        static SegPtr create_tail();
    };

    //// ────────────────────────────────────────────────
    //// 无锁链表辅助函数 (已去除 detail 命名空间)
    //// ────────────────────────────────────────────────
    //inline void link_node(SegPtr new_node, SegPtr pred, SegPtr succ) noexcept {
    //    new_node->primary_path.prev.store(pred, std::memory_order_relaxed);
    //    new_node->primary_path.next.store(succ, std::memory_order_relaxed);
    //}

    //inline bool cas_insert_after(SegPtr pred, SegPtr new_node) noexcept {
    //    SegPtr succ = pred->primary_path.next.load(std::memory_order_acquire);
    //    link_node(new_node, pred, succ);

    //    if (pred->primary_path.next.compare_exchange_strong(
    //        succ, new_node,
    //        std::memory_order_release,
    //        std::memory_order_relaxed)) {
    //        succ->primary_path.prev.store(new_node, std::memory_order_release);
    //        return true;
    //    }
    //    return false;
    //}

    //// 根据 start 位置插入到以 head 为首的有序链表
    //inline void concurrent_insert_segment(SegPtr head, SegPtr new_node) {
    //    while (true) {
    //        SegPtr pred = head;
    //        SegPtr succ = pred->primary_path.next.load(std::memory_order_acquire);
    //        while (!succ->is_tail() && succ->start < new_node->start) {
    //            pred = succ;
    //            succ = succ->primary_path.next.load(std::memory_order_acquire);
    //        }
    //        if (cas_insert_after(pred, new_node))
    //            break;
    //    }
    //}

    // ────────────────────────────────────────────────
    // Block
    // ────────────────────────────────────────────────
    class Block : public std::enable_shared_from_this<Block> {
    public:
        ChrName   ref_chr;          // guard: rw
        ChrHeadMap anchors;         // guard: rw (每条 chr 的 head)
        mutable std::shared_mutex rw;

        static BlockPtr make(std::size_t hint = 1);
        static BlockPtr create_empty(const ChrName& chr, std::size_t hint = 1);
        static BlockPtr create_from_region(const Region& region,
            Strand sd = Strand::FORWARD,
            AlignRole rl = AlignRole::PRIMARY);
        static BlockPtr create_from_match(const Match& match);

        Block() = default;
        ~Block() = default;
    };

    // ────────────────────────────────────────────────
    // GenomeEnd
    // ────────────────────────────────────────────────
    struct GenomeEnd {
        std::unique_ptr<Segment> head_holder;
        std::unique_ptr<Segment> tail_holder;

        SegPtr head;
        SegPtr tail;

        mutable std::shared_mutex rw;

        GenomeEnd() {
            head_holder.reset(Segment::create_head());
            tail_holder.reset(Segment::create_tail());
            head = head_holder.get();
            tail = tail_holder.get();
            head->primary_path.next.store(tail, std::memory_order_relaxed);
            tail->primary_path.prev.store(head, std::memory_order_relaxed);
        }

        GenomeEnd(GenomeEnd&&)            noexcept = default;
        GenomeEnd& operator=(GenomeEnd&&) noexcept = default;
        GenomeEnd(const GenomeEnd&) = delete;
        GenomeEnd& operator=(const GenomeEnd&) = delete;

        std::pair<SegPtr, SegPtr> find_surrounding(uint_t beg, uint_t end);
    };

    // ────────────────────────────────────────────────
    // RaMeshGenomeGraph / RaMeshMultiGenomeGraph
    // ────────────────────────────────────────────────
    class RaMeshGenomeGraph {
    public:
        RaMeshGenomeGraph() = default;
        explicit RaMeshGenomeGraph(const SpeciesName& sp);
        RaMeshGenomeGraph(const SpeciesName& sp, const std::vector<ChrName>& chrs);

        SpeciesName                             species_name;
        std::unordered_map<ChrName, GenomeEnd>  chr2end;   // guard: rw
        mutable std::shared_mutex               rw;        // 多读单写

        size_t debug_print(std::ostream& os = std::cout, bool show_secondary = false) const;
    };

    class RaMeshMultiGenomeGraph {
    public:
        std::unordered_map<SpeciesName, RaMeshGenomeGraph> species_graphs; // guard: rw
        std::vector<WeakBlock> blocks;                                     // guard: rw
        mutable std::shared_mutex rw;            // 多读单写

        RaMeshMultiGenomeGraph() = default;
        explicit RaMeshMultiGenomeGraph(std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map);

        void insertClusterIntoGraph(SpeciesName ref, SpeciesName qry, const MatchCluster& cluster);
        void debug_print(std::ostream& os = std::cout, bool show_secondary = false) const;
    };

} // namespace RaMesh
#endif /* RAMESH_H */
