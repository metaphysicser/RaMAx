#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <filesystem>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <unordered_map>
#include <set>

#include <spdlog/spdlog.h>
#include "../submodule/hal/api/inc/hal.h"

#include "ramesh.h"
#include "SeqPro.h"
#include "data_process.h"  // 使用现有的NewickParser

namespace RaMesh {

namespace hal_converter {

    // ========================================
    // 祖先节点信息结构
    // ========================================

    struct AncestorNode {
        std::string node_name;                              // 祖先节点名称
        std::vector<std::string> descendant_leaves;         // 该祖先下的所有叶节点
        std::map<std::string, std::string> reconstructed_sequences;  // chr_name -> sequence
        std::map<std::string, hal_size_t> sequence_lengths; // chr_name -> length

        // 树结构信息
        int tree_depth;                                     // 在树中的深度
        std::string parent_name;                            // 父节点名称
        std::vector<std::string> children_names;            // 子节点名称
        bool is_generated_root;                             // 是否为自动生成的根节点
        double branch_length;                               // 到父节点的分支长度

        AncestorNode() : tree_depth(0), is_generated_root(false), branch_length(0.0) {}
    };

    // ========================================
    // 系统发育树解析和处理
    // ========================================

    /**
     * 解析Newick树并提取祖先节点信息
     * @param newick_tree Newick格式的系统发育树字符串
     * @param seqpro_managers 序列管理器映射，用于验证叶节点名称
     * @return 祖先节点列表
     */
    std::vector<AncestorNode> parsePhylogeneticTree(
        const std::string& newick_tree,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 检查树是否有显式根节点，如果没有则自动生成
     * @param parser 已解析的Newick树
     * @param seqpro_managers 序列管理器映射
     * @return 是否添加了新的根节点
     */
    bool ensureRootNode(NewickParser& parser,
                       const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 从NewickParser提取祖先节点信息
     * @param parser 已解析的Newick树
     * @param seqpro_managers 序列管理器映射
     * @return 祖先节点列表
     */
    std::vector<AncestorNode> extractAncestorNodes(
        const NewickParser& parser,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 验证叶节点名称是否与提供的物种名称匹配
     * @param parser 已解析的Newick树
     * @param seqpro_managers 序列管理器映射
     * @return 验证结果和错误信息
     */
    std::pair<bool, std::string> validateLeafNames(
        const NewickParser& parser,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    // ========================================
    // HAL基础结构创建
    // ========================================

    /**
     * 创建基本的基因组结构（仅叶节点）
     * @param alignment HAL alignment对象
     * @param seqpro_managers 序列管理器映射
     */
    void createBasicGenomeStructure(
        hal::AlignmentPtr alignment,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 设置基因组的序列维度和DNA数据
     * @param alignment HAL alignment对象
     * @param seqpro_managers 序列管理器映射
     */
    void setupGenomeSequences(
        hal::AlignmentPtr alignment,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);

    /**
     * 根据祖先节点信息创建祖先基因组
     * @param alignment HAL alignment对象
     * @param ancestor_nodes 祖先节点列表
     * @param parser NewickParser对象，用于获取叶节点的分支长度
     */
    void createAncestorGenomes(
        hal::AlignmentPtr alignment,
        const std::vector<AncestorNode>& ancestor_nodes,
        const NewickParser& parser);

    /**
     * 应用系统发育树结构到HAL alignment
     * @param alignment HAL alignment对象
     * @param parser 已解析的Newick树
     */
    void applyPhylogeneticTree(
        hal::AlignmentPtr alignment,
        const NewickParser& parser);

    // ========================================
    // 验证和工具函数
    // ========================================

    /**
     * 验证HAL文件的完整性
     * @param alignment HAL alignment对象
     * @param hal_path HAL文件路径
     */
    void validateHalFile(hal::AlignmentPtr alignment, const FilePath& hal_path);

    /**
     * 生成默认的根节点名称
     * @param leaf_names 叶节点名称列表
     * @return 生成的根节点名称
     */
    std::string generateRootName(const std::vector<std::string>& leaf_names);

} // namespace hal_converter

} // namespace RaMesh
