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

#include <spdlog/spdlog.h>
#include "../submodule/hal/api/inc/hal.h"

#include "ramesh.h"
#include "SeqPro.h"

namespace RaMesh {

// ========================================
// HAL导出辅助函数
// ========================================

namespace hal_converter {

    // 检查并修正Newick树的根节点命名
    std::string correctNewickTree(const std::string& newick_tree);
    
    // 创建基本的基因组结构
    void createBasicGenomeStructure(hal::AlignmentPtr alignment, 
                                  const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);
    
    // 设置基因组的序列维度和DNA数据
    void setupGenomeSequences(hal::AlignmentPtr alignment, 
                             const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);
    
    // 为内部节点设置基本数据
    void setupInternalNodes(hal::AlignmentPtr alignment);
    
    // 应用Newick树到已创建的基因组结构
    void applyNewickTree(hal::AlignmentPtr alignment, const std::string& newick_tree);
    
    // ========================================
    // 祖先重建相关数据结构和函数
    // ========================================
    
    // 祖先节点信息
    struct AncestorNode {
        std::string node_name;
        std::vector<std::string> descendant_leaves;  // 该祖先下的所有叶节点
        std::map<std::string, std::string> reconstructed_sequences;  // chr_name -> sequence
        std::map<std::string, hal_size_t> sequence_lengths;  // chr_name -> length
    };
    
    // 祖先比对块
    struct AncestorBlock {
        std::string ancestor_name;
        std::map<std::string, std::vector<RaMesh::SegPtr>> descendant_segments;  // species -> segments
        std::map<std::string, std::string> ancestor_sequences;  // chr_name -> consensus sequence
    };
    
    // 解析系统发育树，获取所有内部节点
    std::vector<AncestorNode> parsePhylogeneticTree(const std::string& newick_tree, 
                                                   const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);
    
    // 收集与特定祖先相关的比对块
    std::vector<AncestorBlock> collectAncestorBlocks(const std::vector<RaMesh::WeakBlock>& blocks,
                                                    const std::vector<AncestorNode>& ancestor_nodes);
    
    // Majority rule consensus重建
    std::string reconstructConsensusSequence(const std::vector<std::string>& sequences);
    
    // 重建祖先序列
    void reconstructAncestorSequences(std::vector<AncestorBlock>& ancestor_blocks,
                                     const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers);
    
    // 将祖先数据写入HAL
    void writeAncestorDataToHAL(hal::AlignmentPtr alignment,
                               const std::vector<AncestorBlock>& ancestor_blocks);
    
    // 验证和记录HAL文件状态
    void validateHalFile(hal::AlignmentPtr alignment, const FilePath& hal_path);

} // namespace hal_converter

} // namespace RaMesh
