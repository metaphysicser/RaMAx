// =============================================================
//  File: ramesh.h   –  Public interface (v0.6‑alpha, 2025‑06‑25)
//  ✧  Header for the lock‑free multi‑genome graph  ✧
// =============================================================
#ifndef RAMESH_H
#define RAMESH_H

#include <atomic>
#include <memory>
#include <unordered_map>
#include <vector>
#include <map>
#include <shared_mutex>
#include <iostream>
#include <chrono>
#include <algorithm>

#include "config.hpp"
#include "anchor.h"

namespace RaMesh {

    // ────────────────────────────────────────────────
    // Forward declarations
    // ────────────────────────────────────────────────
    class  RaMeshNode;
    class  Block;
    struct Segment;

    using BlockPtr = std::shared_ptr<Block>;
    using WeakBlock = std::weak_ptr<Block>;
    using SegPtr = Segment*;

    // ────────────────────────────────────────────────
    // Helper types / hashing
    // ────────────────────────────────────────────────
    using SpeciesChrPair = std::pair<SpeciesName, ChrName>;

    struct SpeciesChrPairHash {
        std::size_t operator()(const SpeciesChrPair& k) const noexcept {
            std::size_t h1 = std::hash<std::string>{}(k.first);
            std::size_t h2 = std::hash<std::string>{}(k.second);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };

    using ChrHeadMap = std::unordered_map<SpeciesChrPair, SegPtr, SpeciesChrPairHash>;

    // ────────────────────────────────────────────────
    // Intrusive concurrent list node
    // ────────────────────────────────────────────────
    struct RaMeshPath {
        std::atomic<SegPtr> next{ nullptr };
        std::atomic<SegPtr> prev{ nullptr };
    };

    // ────────────────────────────────────────────────
    // Domain enums
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
        BlockPtr   parent_block;

        mutable std::shared_mutex rw;        // guards non‑list fields

        // ――― predicates ―――
        [[nodiscard]] bool isHead()      const noexcept { return seg_role == SegmentRole::HEAD; }
        [[nodiscard]] bool isTail()      const noexcept { return seg_role == SegmentRole::TAIL; }
        [[nodiscard]] bool isSegment()   const noexcept { return seg_role == SegmentRole::SEGMENT; }
        [[nodiscard]] bool isPrimary()   const noexcept { return align_role == AlignRole::PRIMARY; }
        [[nodiscard]] bool isSecondary() const noexcept { return align_role == AlignRole::SECONDARY; }

        // ――― factories ―――
        static SegPtr create(uint_t start, uint_t len, Strand sd,
            Cigar_t cg, AlignRole rl, SegmentRole sl,
            const BlockPtr& bp);

        static SegPtr createFromRegion(Region& region, Strand sd,
            Cigar_t cg, AlignRole rl, SegmentRole sl,
            const BlockPtr& bp);

        static SegPtr createHead();
        static SegPtr createTail();

        // ――― utilities ―――
        static void  linkChain(const std::vector<SegPtr>& segments);
        
        // ――― deletion utilities ―――
        static void unlinkSegment(SegPtr segment);
        static void deleteSegment(SegPtr segment);
        static void deleteBatch(const std::vector<SegPtr>& segments);
    };

    // ────────────────────────────────────────────────
    // Block – a bi‑species alignment block
    // ────────────────────────────────────────────────
    class Block : public std::enable_shared_from_this<Block> {
    public:
        ChrName   ref_chr;      // guard: rw
        ChrHeadMap anchors;     // guard: rw (head sentinel of every chr)
        mutable std::shared_mutex rw;

        static BlockPtr create(std::size_t hint = 1);
        static BlockPtr createEmpty(const ChrName& chr, std::size_t hint = 1);

        // Convenience helper – create both ref&qry segments, register anchors
        static std::pair<SegPtr, SegPtr> createSegmentPair(const Match& match,
            const SpeciesName& ref_name,
            const SpeciesName& qry_name,
            const ChrName& ref_chr,
            const ChrName& qry_chr,
            const BlockPtr& blk);

        static std::pair<SegPtr, SegPtr> createSegmentPair(const Anchor& anchor,
            const SpeciesName& ref_name,
            const SpeciesName& qry_name,
            const ChrName& ref_chr,
            const ChrName& qry_chr,
            const BlockPtr& blk);
            
        // ――― deletion methods ―――
        void removeAllSegments();
        void removeSegmentsBySpecies(const SpeciesName& species);
        bool removeSegment(const SpeciesName& species, const ChrName& chr);

        Block() = default;
        ~Block() = default;
    };

    // ────────────────────────────────────────────────
    // GenomeEnd  – head/tail sentinels per chromosome
    // ────────────────────────────────────────────────
    class GenomeEnd {
    public:
        static constexpr uint_t kSampleStep = 10000;
        GenomeEnd();
        GenomeEnd(GenomeEnd&&)            noexcept = default;
        GenomeEnd& operator=(GenomeEnd&&) noexcept = default;
        GenomeEnd(const GenomeEnd&) = delete;
        GenomeEnd& operator=(const GenomeEnd&) = delete;

        std::pair<SegPtr, SegPtr> findSurrounding(uint_t beg, uint_t end);

        // Atomically splice an already linked chain [segments.front(), segments.back()] into list
        void spliceSegmentChain(const std::vector<SegPtr>& segments,
            uint_t beg, uint_t end);

        void clearAllSegments();
        
        // ――― deletion methods ―――
        bool removeSegment(SegPtr segment);
        bool removeSegmentRange(uint_t start, uint_t end);
        void removeBatch(const std::vector<SegPtr>& segments);
        void invalidateSampling(uint_t start, uint_t end);

        void removeOverlap();

        SegPtr head{ nullptr };
        SegPtr tail{ nullptr };

        mutable std::shared_mutex rw;


    private:
        std::unique_ptr<Segment> head_holder;
        std::unique_ptr<Segment> tail_holder;

        std::vector<SegPtr> sample_vec;


        static bool spliceRange(SegPtr prev, SegPtr next,
            SegPtr first, SegPtr last);

        void ensureSampleSize(uint_t pos);          // 自动扩容
    
    public:
        void updateSampling(const std::vector<SegPtr>& segs); // 插入后修补
        
    private:
    };

    // ────────────────────────────────────────────────
    // RaMeshGenomeGraph / RaMeshMultiGenomeGraph
    // ────────────────────────────────────────────────
    class RaMeshGenomeGraph {
    public:
        RaMeshGenomeGraph() = default;
        explicit RaMeshGenomeGraph(const SpeciesName& sp);
        RaMeshGenomeGraph(const SpeciesName& sp, const std::vector<ChrName>& chrs);

        size_t debugPrint(bool show_detail) const;

        SpeciesName                             species_name;
        std::unordered_map<ChrName, GenomeEnd>  chr2end;   // guard: rw
        mutable std::shared_mutex               rw;        // multi‑reader / single‑writer
    };

    class RaMeshMultiGenomeGraph {
    public:
        // RaMeshMultiGenomeGraph() = default;
        explicit RaMeshMultiGenomeGraph(std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_map);

        explicit RaMeshMultiGenomeGraph(std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_map);

        void insertClusterIntoGraph(SpeciesName ref_name, SpeciesName qry_name,
            const MatchCluster& cluster);

        void insertAnchorIntoGraph(SpeciesName ref_name, SpeciesName qry_name,
			const AnchorVec& anchor_vec);

        void debugPrint(bool show_detail) const;
        
        // 图正确性验证函数
        bool verifyGraphCorrectness(bool verbose = false) const;

        void mergeMultipleGraphs(const SpeciesName& ref_name, uint_t thread_num);

        std::unordered_map<SpeciesName, RaMeshGenomeGraph> species_graphs; // guard: rw
        std::vector<WeakBlock>                             blocks;         // guard: rw
        mutable std::shared_mutex                          rw;             // multi‑reader / single‑writer

        void exportToMaf(const FilePath& maf_path, const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers, bool only_primary, bool is_pairwise) const;

        
        // ――― high-performance deletion methods ―――
        bool removeBlock(const BlockPtr& block);
        bool removeBlockById(size_t block_index);
        void removeBlocksBatch(const std::vector<BlockPtr>& blocks);
        size_t removeExpiredBlocks();
        
        void removeSpecies(const SpeciesName& species);
        void removeChromosome(const SpeciesName& species, const ChrName& chr);
        void clearAllGraphs();
        
        // ――― garbage collection and maintenance ―――
        size_t compactBlockPool();
        void optimizeGraphStructure();
        
        // ――― deletion statistics and diagnostics ―――
        struct DeletionStats {
            size_t segments_removed = 0;
            size_t blocks_removed = 0;
            size_t expired_blocks_cleaned = 0;
            size_t sampling_updates = 0;
            std::chrono::microseconds total_time{0};
        };
        
        DeletionStats performMaintenance(bool full_gc = false);
    };

} // namespace RaMesh

#endif /* RAMESH_H */