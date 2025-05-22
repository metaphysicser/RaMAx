#include "anchor.h"

// 构造函数（注释掉了）：初始化参考物种、查询物种及其锚点集合，并重建 R 树。
//PairGenomeAnchor::PairGenomeAnchor(SpeciesName ref,
//    SpeciesName query,
//    AnchorPtrListVec init_anchors)
//    : ref_species(ref)
//    , query_species(query)
//    , anchors(std::move(init_anchors))
//{
//    rebuild(); // 构建 R 树索引
//}

// 重建 R 树索引：将 anchors 中所有 Anchor 加入 R 树中。
//
//void PairGenomeAnchor::rebuild() {
//    // 1) 清空已有的 R 树索引
//    anchor_rtree.clear();
//
//    // 2) 遍历 anchors 的每个 list，插入 AnchorPtr 到 R-tree 中
//    for (auto const& list_group : anchors) {
//        for (auto const& ap : list_group) {
//            anchor_rtree.insert({ make_box(*ap), ap });
//        }
//    }
//}

// 将一个 Anchor 插入 R 树中
//void PairGenomeAnchor::insertAnchorRtree(AnchorPtr ap) {
//    anchor_rtree.insert({ make_box(*ap), ap });
//}

// 从 R 树中移除一个 Anchor
//void PairGenomeAnchor::removeAnchorRtree(AnchorPtr ap) {
//    // 从 R 树中删除 Anchor
//    anchor_rtree.remove({ make_box(*ap), ap });
//}

// 查询与给定区域相交的所有锚点
//std::vector<AnchorPtr> PairGenomeAnchor::query(const AnchorBox& region) const {
//    // 在 R 树中查找与 region 相交的 AnchorEntry
//    std::vector<AnchorEntry> hits;
//    anchor_rtree.query(bgi::intersects(region),
//        std::back_inserter(hits));
//
//    // 提取命中的 AnchorPtr 到输出向量中
//    std::vector<AnchorPtr> out;
//    out.reserve(hits.size());
//    for (auto const& e : hits) {
//        out.push_back(e.second);
//    }
//    return out;
//}

// ------------------------------------------------------------------
// 保存锚点数据到文件：将二维列表结构 anchors 写入二进制文件
// ------------------------------------------------------------------
bool saveAnchors(const std::string& filename,
    const AnchorPtrListVec& anchors)
{
    std::ofstream os(filename, std::ios::binary); // 以二进制方式打开文件
    if (!os) return false; // 打开失败

    cereal::BinaryOutputArchive oar(os); // 创建 cereal 的输出归档对象

    // 1) 写入外层 vector 的大小
    uint64_t outer = anchors.size();
    oar(outer);

    // 2) 对每个内部 list：
    for (const auto& lst : anchors) {
        uint64_t inner = lst.size(); // 写入 list 大小
        oar(inner);

        // 写入每个 Anchor 对象的内容
        for (const auto& ap : lst) {
            oar(*ap); // 保存对象内容，而不是指针本身
        }
    }
    return static_cast<bool>(os); // 返回写入是否成功
}

// ------------------------------------------------------------------
// 从文件读取锚点数据：与 saveAnchors 格式对应
// ------------------------------------------------------------------
bool loadAnchors(const std::string& filename,
    AnchorPtrListVec& anchors)
{
    std::ifstream is(filename, std::ios::binary); // 打开二进制文件
    if (!is) return false;

    cereal::BinaryInputArchive iar(is); // 创建 cereal 的输入归档对象

    uint64_t outer;
    iar(outer); // 读取外层 vector 大小

    anchors.clear();
    anchors.reserve(static_cast<size_t>(outer)); // 预分配空间

    // 读取每个内部 list
    for (uint64_t i = 0; i < outer; ++i) {
        uint64_t inner;
        iar(inner); // 读取内部 list 大小

        AnchorPtrList lst;
        for (uint64_t j = 0; j < inner; ++j) {
            auto ap = std::make_shared<Anchor>(); // 创建新的 Anchor 智能指针
            iar(*ap); // 从文件读取 Anchor 内容到 ap 中
            lst.push_back(std::move(ap));
        }
        anchors.emplace_back(std::move(lst));
    }
    return static_cast<bool>(is); // 返回读取是否成功
}
