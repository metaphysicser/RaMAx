#pragma once
#include <memory>
#include <unordered_map>
#include <atomic>
#include "config.hpp"
#include "anchor.h"

/* ---------- 基础别名 ---------- */
using NodePtr = std::shared_ptr<class RaMeshNode>;
using WeakNodePtr = std::weak_ptr<RaMeshNode>;
using PathPtr = std::shared_ptr<class RaMeshPath>;

/* ---------- 基类 ---------- */
struct RaMeshNode : std::enable_shared_from_this<RaMeshNode> {
    [[nodiscard]] virtual bool is_head()  const noexcept = 0;
    [[nodiscard]] virtual bool is_tail()  const noexcept = 0;
    [[nodiscard]] virtual bool is_block() const noexcept = 0;
    virtual ~RaMeshNode() = default;

protected:
    // 禁用拷贝，防止无意复制 node
    RaMeshNode() = default;
    RaMeshNode(const RaMeshNode&) = delete;
    RaMeshNode& operator=(const RaMeshNode&) = delete;
};

/* ---------- 路径节点 ---------- */
struct RaMeshPath {
    // 不拥有邻接节点，只做链式导航；用 weak_ptr 防环并省原子计数
    std::atomic<WeakNodePtr> next_block;
    std::atomic<WeakNodePtr> last_block;

    RaMeshPath() = default;
};

/* ---------- Segment ---------- */
struct Segment {
    Region     region;      // 该基因组在此块上的坐标
    Cigar_t    cigar;
    Strand     strand{};
    RaMeshPath path;        // 在每个基因组链中的位置

    Segment() = default;
};

/* ---------- Block & 哨兵 ---------- */
struct Block final : RaMeshNode {
    ChrName ref_chr;

    // 基因组数通常可估计 → 提前 reserve 减少 rehash
    std::pmr::unordered_map<ChrName, Segment> anchors;

    [[nodiscard]] bool is_head()  const noexcept override { return false; }
    [[nodiscard]] bool is_tail()  const noexcept override { return false; }
    [[nodiscard]] bool is_block() const noexcept override { return true; }

    Block() = default;
    ~Block() = default;
};

struct HeadBlock final : RaMeshNode {
    std::atomic<WeakNodePtr> next_block;

    [[nodiscard]] bool is_head()  const noexcept override { return true; }
    [[nodiscard]] bool is_tail()  const noexcept override { return false; }
    [[nodiscard]] bool is_block() const noexcept override { return false; }

    HeadBlock() = default;
    ~HeadBlock() = default;
};

struct TailBlock final : RaMeshNode {
    std::atomic<WeakNodePtr> last_block;

    [[nodiscard]] bool is_head()  const noexcept override { return false; }
    [[nodiscard]] bool is_tail()  const noexcept override { return true; }
    [[nodiscard]] bool is_block() const noexcept override { return false; }

    TailBlock() = default;
    ~TailBlock() = default;
};

/* ---------- 染色体端点索引 ---------- */
struct GenomeEndBlock {
    WeakNodePtr head;   // HeadBlock
    WeakNodePtr tail;   // TailBlock
};

/* ---------- 物种级图 ---------- */
struct RaMeshGenomeGraph {
    SpeciesName species_name;
    std::unordered_map<ChrName, GenomeEndBlock> end_points;

    // reserve 可根据已知染色体数优化
    explicit RaMeshGenomeGraph(SpeciesName sp, std::size_t chr_est = 64)
        : species_name(std::move(sp)) {
        end_points.reserve(chr_est);
    }
};
