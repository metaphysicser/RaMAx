/*******************************************************
*  RaMesh Graph – Segment-centric model  (v0.4-rev)    *
*  — 2025-06-24 —                                      *
*******************************************************/
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

namespace RaMesh
{

    /* ===== 前向声明 ===== */
    class  RaMeshNode;
    class  Block;
    struct Segment;

    using BlockPtr = std::shared_ptr<Block>;
    using WeakBlock = std::weak_ptr<Block>;
    using SegPtr = Segment*;

    /* ===================================================== *
     * 1. SegAtom / SegAtomPtr                               *
     * ===================================================== */
    struct SegAtom {
        std::atomic<SegPtr> p{ nullptr };

        SegAtom() = default;
        explicit constexpr SegAtom(SegPtr s) noexcept : p(s) {}

        SegAtom(const SegAtom&) = delete;
        SegAtom& operator=(const SegAtom&) = delete;

        SegAtom(SegAtom&& other) noexcept {
            p.store(other.p.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
        }
        SegAtom& operator=(SegAtom&& other) noexcept {
            p.store(other.p.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
            return *this;
        }

        /* 与 std::atomic<SegPtr> 接口一致 */
        inline void   store(SegPtr v,
            std::memory_order mo = std::memory_order_seq_cst) noexcept {
            p.store(v, mo);
        }
        [[nodiscard]] inline SegPtr load(
            std::memory_order mo = std::memory_order_seq_cst) const noexcept {
            return p.load(mo);
        }
        inline SegPtr exchange(SegPtr v,
            std::memory_order mo = std::memory_order_seq_cst) noexcept {
            return p.exchange(v, mo);
        }

        /* 方便读取：允许隐式转成裸指针 */
        constexpr operator SegPtr() const noexcept { return p.load(std::memory_order_relaxed); }
    };

    using SegAtomPtr = std::shared_ptr<SegAtom>;

    [[nodiscard]] inline SegAtomPtr make_atom(SegPtr p) { return std::make_shared<SegAtom>(p); }

    /* ===================================================== *
     * 2. 辅助别名 / 哈希                                     *
     * ===================================================== */
    using SpeciesChrPair = std::pair<SpeciesName, ChrName>;

    struct SpeciesChrPairHash {
        constexpr std::size_t operator()(const SpeciesChrPair& k) const noexcept {
            std::size_t h1 = std::hash<std::string>{}(k.first);
            std::size_t h2 = std::hash<std::string>{}(k.second);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };

    using ChrSegMap = std::unordered_map<SpeciesChrPair,
        SegAtomPtr,
        SpeciesChrPairHash>;

    /* ===================================================== *
     * 3. 链表节点基类                                       *
     * ===================================================== */
    struct RaMeshPath { SegAtomPtr next{ nullptr }; SegAtomPtr prev{ nullptr }; };

    enum class AlignRole : uint8_t { PRIMARY = 0, SECONDARY = 1 };
    enum class SegmentRole: uint8_t {SEGMENT = 0, HEAD = 1, TAIL = 2};


    /* ===================================================== *
     * 4. Segment / Sentinel                                 *
     * ===================================================== */
    class Segment{
    public:
        uint_t    start{ 0 };
        uint_t    length{ 0 };
        Cigar_t   cigar;
        Strand    strand{ Strand::FORWARD };
        AlignRole align_role{ AlignRole::PRIMARY };
		SegmentRole seg_role{ SegmentRole::SEGMENT };

        RaMeshPath  primary_path;
        // RaMeshPath  secondary_path;

        WeakBlock   parent_block;              // guard: anchors / graph
        mutable std::shared_mutex rw;          // guard: fields inside Segment

        bool is_head()    const noexcept { return seg_role == SegmentRole::HEAD; }
        bool is_tail()    const noexcept { return seg_role == SegmentRole::TAIL; }
        bool is_segment() const noexcept { return seg_role == SegmentRole::SEGMENT; }
		bool is_primary() const noexcept { return align_role == AlignRole::PRIMARY; }
		bool is_secondary() const noexcept { return align_role == AlignRole::SECONDARY; }

        /* factory helpers */
/* factory helpers（原来三条保留） */
        static SegPtr create(uint_t start, uint_t len, Strand sd,
            Cigar_t cg, AlignRole rl, SegmentRole sl, const BlockPtr& bp);
        static SegPtr create_from_region(Region& region, Strand sd,
            Cigar_t cg, AlignRole rl, SegmentRole sl, const BlockPtr& bp);

        /* ★ 新增哨兵工厂 ★ */
        static SegPtr create_head();
        static SegPtr create_tail();

    };

    
    /* ===================================================== *
     * 5. Block                                              *
     * ===================================================== */
    class Block : public std::enable_shared_from_this<Block>
    {
    public:
        ChrName   ref_chr;                       // guard: rw
        ChrSegMap anchors;                       // guard: rw
        mutable std::shared_mutex rw;            // guard: ref_chr / anchors

        static BlockPtr make(std::size_t hint = 1);
        static BlockPtr create_empty(const ChrName& chr, std::size_t hint = 1);
        static BlockPtr create_from_region(const Region& region,
            Strand sd = Strand::FORWARD,
            AlignRole rl = AlignRole::PRIMARY);
        static BlockPtr create_from_match(const Match& match);

        Block() = default;
        ~Block() = default;
    };

    /* ===================================================== *
     * 6. GenomeEnd                                          *
     * ===================================================== */
    struct GenomeEnd {
        SegAtomPtr head;
        SegAtomPtr tail;

        mutable std::shared_mutex rw;

        GenomeEnd() {
            head = make_atom(Segment::create_head());
            tail = make_atom(Segment::create_tail());
            head->load()->primary_path.next = tail;
            tail->load()->primary_path.prev = head;
        }

        /* ★ 显式移动，只移动指针，不动锁 ★ */
        GenomeEnd(GenomeEnd&& other) noexcept
            : head(std::move(other.head)),
            tail(std::move(other.tail)) {
        }

        GenomeEnd& operator=(GenomeEnd&& other) noexcept {
            if (this != &other) {
                head = std::move(other.head);
                tail = std::move(other.tail);
            }
            return *this;
        }

        /* 禁止拷贝 */
        GenomeEnd(const GenomeEnd&) = delete;
        GenomeEnd& operator=(const GenomeEnd&) = delete;

        std::pair<SegAtomPtr, SegAtomPtr>
            findSurroundingSegmentRange(uint_t beg, uint_t end);
    };

    /* ===================================================== *
     * 7. GenomeGraph / MultiGenomeGraph                     *
     * ===================================================== */
    class RaMeshGenomeGraph {
    public:
        RaMeshGenomeGraph() = default;
        explicit RaMeshGenomeGraph(const SpeciesName& sp);
        RaMeshGenomeGraph(const SpeciesName& sp,
            const std::vector<ChrName>& chrs);

        SpeciesName                             species_name;
        std::unordered_map<ChrName, GenomeEnd>  chr2end;   // guard: rw
        mutable std::shared_mutex               rw;        // 多读单写

        void debug_print(std::ostream& os = std::cout,
            bool show_secondary = false) const;
    };

    class RaMeshMultiGenomeGraph {
    public:
        std::unordered_map<SpeciesName, RaMeshGenomeGraph> species_graphs; // guard: rw
        std::vector<WeakBlock> blocks;                                     // guard: rw
        mutable std::shared_mutex rw;            // 多读单写

        RaMeshMultiGenomeGraph() = default;
        RaMeshMultiGenomeGraph(std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map);

        void insertClusterIntoGraph(SpeciesName ref, SpeciesName qry,
            const MatchCluster& cluster);
        void debug_print(std::ostream& os = std::cout,
            bool show_secondary = false) const;
    };

} // namespace RaMesh
#endif   // RAMESH_H
