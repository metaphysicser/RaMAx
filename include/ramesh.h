/*******************************************************
*  RaMesh Graph – Segment-centric model  (v0.4)        *
*  — 2025-06-20 —                                      *
*  1) Segment/Head/Tail 全部是 RaMeshNode 派生类。     *
*  2) Block 退化为“Segment 池”容器，仅管理内存 & 索引。 *
*  3) 主(Primary)/次(Secondary) 比对仍保留双链。        *
*******************************************************/

#ifndef RAMESH_H
#define RAMESH_H
#include <memory>
#include <atomic>
#include <unordered_map>
#include <vector>
#include <memory_resource>
#include <shared_mutex>
#include "config.hpp"
#include "anchor.h"

namespace RaMesh {

    /* ---------- 别名 ---------- */
    class  RaMeshNode;
    class  Block;                    // 前置声明
    struct Segment;

    using BlockPtr = std::shared_ptr<Block>;
    using WeakBlock = std::weak_ptr<Block>;

    using SegPtr = Segment*;                 // 生命周期由所属 Block 控制
    using SegAtom = std::atomic<SegPtr>;      // 原子可见裸指针

    using SegVec = std::vector<SegAtom>;
    using ChrSegVecMap = std::unordered_map<ChrName, SegVec>;


    /* ---------- 链表节点 ---------- */
    struct RaMeshPath {
        SegAtom next{ nullptr };
        SegAtom prev{ nullptr };
    };

    /* ---------- 主/次 比对标记 ---------- */
    enum class AlignRole : uint8_t { PRIMARY = 0, SECONDARY = 1 };

    /* ---------- 抽象节点基类 ---------- */
    class RaMeshNode : public std::enable_shared_from_this<RaMeshNode> {
    public:
        [[nodiscard]] virtual bool is_head()  const noexcept = 0;
        [[nodiscard]] virtual bool is_tail()  const noexcept = 0;
        [[nodiscard]] virtual bool is_segment() const noexcept = 0;
        virtual ~RaMeshNode() = default;

    protected:                      // 禁意外复制
        RaMeshNode() = default;
        RaMeshNode(const RaMeshNode&) = delete;
        RaMeshNode& operator=(const RaMeshNode&) = delete;
    };

    /* ---------- Segment (普通节点) ---------- */
    class Segment final : public RaMeshNode {
    public:
        uint_t start{ 0 };
        uint_t length{ 0 };
        Cigar_t  cigar;
        Strand   strand{ Strand::FORWARD };
        AlignRole role{ AlignRole::PRIMARY };

        RaMeshPath primary_path;
        RaMeshPath secondary_path;

        WeakBlock parent_block;                 // ← 由裸指针改为弱引用
        mutable std::shared_mutex rw;

        [[nodiscard]] bool is_head()    const noexcept override { return false; }
        [[nodiscard]] bool is_tail()    const noexcept override { return false; }
        [[nodiscard]] bool is_segment() const noexcept override { return true; }

        static SegPtr create(uint_t           start,
            uint_t           length,
            const Cigar_t& cigar,
            Strand           strand = Strand::FORWARD,
            AlignRole        role = AlignRole::PRIMARY,
            const BlockPtr& parent = nullptr);

        static SegPtr create(const Region& region,
            Strand strand = Strand::FORWARD,
            AlignRole role = AlignRole::PRIMARY,
            const BlockPtr& parent = nullptr);
    };


    /* ---------- Sentinel 节点 ---------- */
    class HeadSegment final : public RaMeshNode {
    public:
        HeadSegment() {}               // 任意默认
        static std::shared_ptr<HeadSegment> create();

        [[nodiscard]] bool is_head()    const noexcept override { return true; }
        [[nodiscard]] bool is_tail()    const noexcept override { return false; }
        [[nodiscard]] bool is_segment() const noexcept override { return false; }
    };

    class TailSegment final : public RaMeshNode {
    public:
        TailSegment() {}
        static std::shared_ptr<TailSegment> create();
		[[nodiscard]] bool is_tail()    const noexcept override { return true; }
		[[nodiscard]] bool is_head()    const noexcept override { return false; }
        [[nodiscard]] bool is_segment() const noexcept override { return false; }
    };

    /* ---------- Block : Segment 容器 ---------- */
    class Block : public std::enable_shared_from_this<Block> {
    public:
        ChrName ref_chr;

        ChrSegVecMap anchors;

        mutable std::shared_mutex rw;

        /* 工厂：返回 BlockPtr，便于统一持有 */
        static BlockPtr make(std::size_t genome_hint = 1);

        /**
 * @brief 创建一个空 Block，仅设置 ref_chr
 * @param ref_chr     参考染色体名称 (可用于后续插 Segment)
 * @param genome_hint anchors 容器预留容量提示
 * @return BlockPtr   新建 shared_ptr<Block>
 */
        static BlockPtr create_empty(const ChrName& ref_chr,
            std::size_t    genome_hint = 1);

        /**
         * @brief 从一个 Region 初始化 Block：设置 ref_chr 并向 anchors[ref_chr] 插入对应 Segment
         * @param region      参考基因组区间，用于初始化一个 Segment
         * @param strand      Segment 链方向
         * @param role        AlignRole（Primary/Secondary）
         * @param parent      如果非空，则设置新 Segment 的 parent_block，为弱引用
         * @param genome_hint anchors map 预留容量
         * @return BlockPtr   新建 shared_ptr<Block>，其中已插入一个基于 region 的 Segment
         */
        static BlockPtr create_from_region(const Region& region,
            Strand         strand = Strand::FORWARD,
            AlignRole      role = AlignRole::PRIMARY);

        /**
         * @brief 从一个 Match 初始化 Block：
         *  - 设置 ref_chr = match.ref_region.chr_name
         *  - 向 anchors[ref_chr] 插入一个 Segment（基于 ref_region）
         *  - （可扩展：如需同时处理 query_region，可自行再插入）
         * @param match       Match 对象
         * @param parent      如果非空，则设置新 Segment 的 parent_block
         * @param genome_hint anchors map 预留容量
         * @return BlockPtr   新建 shared_ptr<Block>，其中已插入 ref_region 对应 Segment
         */
        static BlockPtr create_from_match(const Match& match);

    private:
        Block() = default;
        ~Block() = default;

    };

    /* ---------- 染色体端点索引 ---------- */
    struct GenomeEnd {
        std::shared_ptr<HeadSegment> head;
        std::shared_ptr<TailSegment> tail;

        /* ctor ①：空端点 */
        GenomeEnd();

        /* ctor ②：直接传入 head / tail */
        GenomeEnd(const std::shared_ptr<HeadSegment>& h,
            const std::shared_ptr<TailSegment>& t);
    };
    /* ---------- 物种级图 ---------- */
    class RaMeshGenomeGraph {
    public:
        /* Constructors */
        RaMeshGenomeGraph() = default;                                                  // ① 空
        explicit RaMeshGenomeGraph(const SpeciesName& sp);                      // ② 物种名
        RaMeshGenomeGraph(const SpeciesName& sp,                                // ③ 物种 + 染色体表
            const std::vector<ChrName>& chromosomes);

        /* data */
        SpeciesName                             species_name;
        std::unordered_map<ChrName, GenomeEnd>  chr2end;
        mutable std::shared_mutex               rw;
    };

    /* ---------- 跨物种多基因组图 ---------- */
    class RaMeshMultiGenomeGraph {
    public:
        /* 物种 → 基因组图 */
        std::unordered_map<SpeciesName, RaMeshGenomeGraph> species_graphs;
        mutable std::shared_mutex                          species_mtx;

        /* 所有 Block 的弱引用（用于跨物种共享） */
        std::vector<WeakBlock> blocks;
        mutable std::mutex     rw;

        RaMeshMultiGenomeGraph() = default;
    };


    /* ---------- 辅助函数 ---------- */
/* 判断主 / 次比对 */
    [[nodiscard]] inline bool is_primary(const Segment& s) noexcept
    {
        return s.role == AlignRole::PRIMARY;
    }

    [[nodiscard]] inline bool is_secondary(const Segment& s) noexcept
    {
        return s.role == AlignRole::SECONDARY;
    }

}

#endif