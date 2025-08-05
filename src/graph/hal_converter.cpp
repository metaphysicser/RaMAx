#include "hal_converter.h"
#include "align.h"
#include <spdlog/spdlog.h>
#include "../../submodule/hal/api/inc/hal.h"

namespace RaMesh {
namespace hal_converter {

    // 检查并修正Newick树的根节点命名
    std::string correctNewickTree(const std::string& newick_tree) {
        if (newick_tree.empty()) {
            spdlog::warn("Empty Newick tree provided");
            return "";
        }
        
        std::string corrected_newick = newick_tree;
        
        // 去除首尾空白字符
        corrected_newick.erase(0, corrected_newick.find_first_not_of(" \t\r\n"));
        corrected_newick.erase(corrected_newick.find_last_not_of(" \t\r\n") + 1);
        
        spdlog::debug("Original Newick tree: '{}'", corrected_newick);
        
        // 检查括号匹配
        int paren_count = 0;
        for (char c : corrected_newick) {
            if (c == '(') paren_count++;
            else if (c == ')') paren_count--;
            
            if (paren_count < 0) {
                spdlog::error("Unmatched closing parenthesis in Newick tree");
                throw std::runtime_error("Invalid Newick tree: unmatched closing parenthesis");
            }
        }
        
        if (paren_count != 0) {
            spdlog::error("Unmatched opening parentheses in Newick tree (count: {})", paren_count);
            throw std::runtime_error("Invalid Newick tree: unmatched opening parentheses");
        }
        
        // 检查是否以分号结尾
        if (!corrected_newick.empty() && corrected_newick.back() != ';') {
            spdlog::warn("Newick tree does not end with semicolon, adding one");
            corrected_newick += ';';
        }
        
        // 检查根节点命名
        if (!corrected_newick.empty() && corrected_newick.back() == ';') {
            // 查找最后一个右括号的位置
            size_t last_paren = corrected_newick.find_last_of(')', corrected_newick.length() - 2);
            if (last_paren != std::string::npos) {
                // 检查右括号后到分号前是否有根节点名称
                std::string after_paren = corrected_newick.substr(last_paren + 1, corrected_newick.length() - last_paren - 2);
                
                // 去除空白字符
                after_paren.erase(0, after_paren.find_first_not_of(" \t\r\n"));
                after_paren.erase(after_paren.find_last_not_of(" \t\r\n") + 1);
                
                // 如果只包含分支长度（冒号开头）或为空，则需要添加根节点名称
                if (after_paren.empty() || after_paren[0] == ':') {
                    // 在分号前插入根节点名称
                    corrected_newick.insert(corrected_newick.length() - 1, "root");
                    spdlog::info("Added root node name to Newick tree");
                }
            }
        }
        
        spdlog::debug("Corrected Newick tree: '{}'", corrected_newick);
        return corrected_newick;
    }
    
    // 创建基本的基因组结构
    void createBasicGenomeStructure(hal::AlignmentPtr alignment, 
                                  const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {
        spdlog::info("Creating basic genome structure for {} species...", seqpro_managers.size());
        
        // 创建根节点基因组
        std::string temp_root_name = "Ancestral";
        hal::Genome* root_genome = alignment->addRootGenome(temp_root_name);
        if (!root_genome) {
            throw std::runtime_error("Failed to create root genome: " + temp_root_name);
        }
        spdlog::info("Created root genome: {}", temp_root_name);
        
        // 设置根基因组的维度
        hal_size_t root_length = 1000;  // 虚拟根序列长度
        std::vector<hal::Sequence::Info> root_dimensions;
        root_dimensions.emplace_back(temp_root_name, root_length, 0, 0);
        root_genome->setDimensions(root_dimensions);
        
        // 为根基因组添加DNA数据
        std::string root_dna(root_length, 'N');
        root_genome->setString(root_dna);
        alignment->closeGenome(root_genome);
        spdlog::info("Root genome setup completed");
        
        // 为每个物种创建叶节点基因组
        for (const auto& [species_name, seq_mgr] : seqpro_managers) {
            hal::Genome* leaf_genome = alignment->addLeafGenome(species_name, temp_root_name, 1.0);
            if (!leaf_genome) {
                throw std::runtime_error("Failed to create leaf genome: " + species_name);
            }
            spdlog::info("Created leaf genome: {}", species_name);
        }
    }
    
    // 设置基因组的序列维度和DNA数据
    void setupGenomeSequences(hal::AlignmentPtr alignment, 
                             const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {
        spdlog::info("Setting up sequence data for all genomes...");
        
        for (const auto& [species_name, seq_mgr_variant] : seqpro_managers) {
            hal::Genome* genome = alignment->openGenome(species_name);
            if (!genome) {
                throw std::runtime_error("Failed to open genome: " + species_name);
            }
            
            // 收集序列信息
            std::vector<hal::Sequence::Info> sequence_dimensions;
            std::visit([&](const auto& seq_mgr) {
                auto seq_names = seq_mgr->getSequenceNames();
                spdlog::info("  Found {} sequences in {}", seq_names.size(), species_name);
                
                for (const auto& seq_name : seq_names) {
                    hal_size_t seq_length = seq_mgr->getSequenceLength(seq_name);
                    sequence_dimensions.emplace_back(seq_name, seq_length, 0, 0);
                    spdlog::debug("    Sequence: {} (length: {})", seq_name, seq_length);
                }
            }, *seq_mgr_variant);
            
            // 设置基因组维度和DNA数据
            if (!sequence_dimensions.empty()) {
                genome->setDimensions(sequence_dimensions);
                
                // 添加DNA数据（暂时用N填充）
                std::string dna_data(genome->getSequenceLength(), 'N');
                genome->setString(dna_data);
                
                spdlog::info("  Setup completed for genome: {} ({} sequences, {} bp)", 
                           species_name, sequence_dimensions.size(), genome->getSequenceLength());
            } else {
                spdlog::warn("  No sequences found for genome: {}", species_name);
            }
            
            alignment->closeGenome(genome);
        }
    }
    
    // 为内部节点设置基本数据
    void setupInternalNodes(hal::AlignmentPtr alignment) {
        spdlog::info("Setting up internal node data...");
        
        std::string root_name = alignment->getRootName();
        if (root_name.empty()) {
            spdlog::warn("No root genome found");
            return;
        }
        
        // 获取所有基因组名称
        std::vector<std::string> all_genomes;
        hal::Genome* root = alignment->openGenome(root_name);
        if (root) {
            // 递归遍历所有基因组
            std::function<void(hal::Genome*)> traverse = [&](hal::Genome* genome) {
                if (!genome) return;
                
                std::string genome_name = genome->getName();
                all_genomes.push_back(genome_name);
                
                // 检查是否需要设置数据
                if (genome->getSequenceLength() == 0) {
                    spdlog::info("  Setting up internal node: {}", genome_name);
                    
                    // 为内部节点创建一个虚拟序列
                    hal_size_t internal_length = 1000;
                    std::vector<hal::Sequence::Info> dimensions;
                    dimensions.emplace_back(genome_name, internal_length, 0, 0);
                    
                    genome->setDimensions(dimensions);
                    
                    // 添加DNA数据
                    std::string dna_data(internal_length, 'N');
                    genome->setString(dna_data);
                    
                    spdlog::debug("    Internal node {} setup: {} bp", genome_name, internal_length);
                }
                
                // 遍历子节点
                for (hal_size_t i = 0; i < genome->getNumChildren(); ++i) {
                    hal::Genome* child = genome->getChild(i);
                    if (child) {
                        traverse(child);
                        alignment->closeGenome(child);
                    }
                }
            };
            
            traverse(root);
            alignment->closeGenome(root);
        }
        
        spdlog::info("Internal node setup completed for {} total genomes", all_genomes.size());
    }

    // 应用Newick树到已创建的基因组结构
    void applyNewickTree(hal::AlignmentPtr alignment, const std::string& newick_tree) {
        spdlog::info("Applying phylogenetic tree...");
        
        if (newick_tree.empty()) {
            spdlog::warn("Empty Newick tree provided, skipping tree application");
            return;
        }
        
        try {
            std::string corrected_newick = correctNewickTree(newick_tree);
            spdlog::info("Using corrected Newick tree: {}", corrected_newick);
            
            if (corrected_newick.empty()) {
                spdlog::warn("Empty corrected Newick tree, skipping tree application");
                return;
            }
            
            alignment->replaceNewickTree(corrected_newick);
            spdlog::info("Phylogenetic tree applied successfully");
            
            // 应用树后，需要为新创建的内部节点设置数据
            setupInternalNodes(alignment);
            
        } catch (const std::exception& e) {
            spdlog::error("Failed to apply Newick tree: {}", e.what());
            spdlog::error("Original tree: '{}'", newick_tree);
            throw std::runtime_error("Invalid Newick tree format: " + std::string(e.what()));
        }
    }
    
    // 解析系统发育树，获取所有内部节点
    std::vector<AncestorNode> parsePhylogeneticTree(const std::string& newick_tree, 
                                                   const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {
        spdlog::info("Parsing phylogenetic tree to identify ancestor nodes...");
        
        std::vector<AncestorNode> ancestor_nodes;
        
        if (newick_tree.empty()) {
            spdlog::info("Empty Newick tree provided, creating star topology");
        } else {
            spdlog::debug("Parsing Newick tree: '{}'", newick_tree);
        }
        
        // TODO: 实现完整的Newick树解析
        // 这里需要解析树结构，识别所有内部节点及其后代叶节点
        // 暂时创建一个简单的根祖先节点作为示例
        
        AncestorNode root_ancestor;
        root_ancestor.node_name = "root";
        
        // 收集所有叶节点作为根的后代
        for (const auto& [species_name, _] : seqpro_managers) {
            root_ancestor.descendant_leaves.push_back(species_name);
        }
        
        ancestor_nodes.push_back(root_ancestor);
        
        spdlog::info("Found {} ancestor nodes to reconstruct", ancestor_nodes.size());
        for (const auto& ancestor : ancestor_nodes) {
            spdlog::debug("  Ancestor '{}' with {} descendants", 
                         ancestor.node_name, ancestor.descendant_leaves.size());
        }
        
        return ancestor_nodes;
    }
    
    // 收集与特定祖先相关的比对块
    std::vector<AncestorBlock> collectAncestorBlocks(const std::vector<RaMesh::WeakBlock>& blocks,
                                                    const std::vector<AncestorNode>& ancestor_nodes) {
        spdlog::info("Collecting alignment blocks for ancestor reconstruction...");
        
        std::vector<AncestorBlock> ancestor_blocks;
        
        for (const auto& ancestor : ancestor_nodes) {
            spdlog::debug("Processing ancestor: {}", ancestor.node_name);
            
            std::map<std::string, std::vector<AncestorBlock>> chr_blocks;  // chr_name -> blocks
            
            // 遍历所有比对块
            for (const auto& weak_block : blocks) {
                auto block = weak_block.lock();
                if (!block) continue;
                
                // 检查这个block是否包含该祖先的后代物种
                std::map<std::string, std::vector<RaMesh::SegPtr>> relevant_segments;
                
                for (const auto& [species_chr, segment] : block->anchors) {
                    const std::string& species = species_chr.first;
                    const std::string& chr = species_chr.second;
                    
                    // 检查是否是该祖先的后代
                    if (std::find(ancestor.descendant_leaves.begin(), 
                                ancestor.descendant_leaves.end(), species) != ancestor.descendant_leaves.end()) {
                        relevant_segments[species].push_back(segment);
                    }
                }
                
                // 如果有足够的物种参与，创建祖先块
                if (relevant_segments.size() >= 2) {  // 至少需要2个物种
                    AncestorBlock ancestor_block;
                    ancestor_block.ancestor_name = ancestor.node_name;
                    ancestor_block.descendant_segments = relevant_segments;
                    
                    ancestor_blocks.push_back(ancestor_block);
                }
            }
            
            spdlog::debug("  Found {} relevant blocks for ancestor '{}'", 
                         ancestor_blocks.size(), ancestor.node_name);
        }
        
        spdlog::info("Collected {} total ancestor blocks", ancestor_blocks.size());
        return ancestor_blocks;
    }
    
    // Majority rule consensus重建
    std::string reconstructConsensusSequence(const std::vector<std::string>& sequences) {
        if (sequences.empty()) return "";
        if (sequences.size() == 1) return sequences[0];
        
        // 找到最长序列的长度
        size_t max_length = 0;
        for (const auto& seq : sequences) {
            max_length = std::max(max_length, seq.length());
        }
        
        std::string consensus;
        consensus.reserve(max_length);
        
        // 对每个位置进行majority voting
        for (size_t pos = 0; pos < max_length; ++pos) {
            std::map<char, int> base_counts;
            
            // 统计每个碱基的出现次数
            for (const auto& seq : sequences) {
                if (pos < seq.length()) {
                    char base = std::toupper(seq[pos]);
                    if (base == 'A' || base == 'T' || base == 'G' || base == 'C') {
                        base_counts[base]++;
                    } else {
                        base_counts['N']++;  // 对于ambiguous碱基
                    }
                } else {
                    base_counts['-']++;  // gap
                }
            }
            
            // 选择出现次数最多的碱基
            char consensus_base = 'N';
            int max_count = 0;
            for (const auto& [base, count] : base_counts) {
                if (count > max_count && base != '-') {  // 忽略gap
                    max_count = count;
                    consensus_base = base;
                }
            }
            
            consensus.push_back(consensus_base);
        }
        
        return consensus;
    }
    
    // 重建祖先序列
    void reconstructAncestorSequences(std::vector<AncestorBlock>& ancestor_blocks,
                                     const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {
        spdlog::info("Reconstructing ancestor sequences...");
        
        for (auto& ancestor_block : ancestor_blocks) {
            // spdlog::debug("Reconstructing sequences for ancestor: {}", ancestor_block.ancestor_name);
            
            // 为每个比对块重建祖先序列
            std::map<std::string, std::vector<std::string>> chr_sequences;  // chr_name -> sequences from different species
            
            // 收集所有后代在该块中的序列
            for (const auto& [species, segments] : ancestor_block.descendant_segments) {
                auto seq_mgr_it = seqpro_managers.find(species);
                if (seq_mgr_it == seqpro_managers.end()) continue;
                
                for (const auto& segment : segments) {
                    std::string chr_name = "unknown";  // TODO: 从segment获取染色体名称
                    
                    // 获取序列数据
                    std::string sequence = std::visit([&](const auto& mgr) -> std::string {
                        // TODO: 使用segment的坐标获取实际序列
                        // return mgr->getSubSequence(chr_name, segment->start, segment->length);
                        return std::string(segment->length, 'N');  // 暂时用N填充
                    }, *seq_mgr_it->second);
                    
                    // 处理反向链
                    if (segment->strand == Strand::REVERSE) {
                        reverseComplement(sequence);
                    }
                    
                    chr_sequences[chr_name].push_back(sequence);
                }
            }
            
            // 对每个染色体重建祖先序列
            for (const auto& [chr_name, sequences] : chr_sequences) {
                if (sequences.empty()) continue;
                
                // 使用majority rule consensus重建祖先序列
                std::string consensus = reconstructConsensusSequence(sequences);
                ancestor_block.ancestor_sequences[chr_name] = consensus;
                
                spdlog::debug("  Reconstructed {} bp for chromosome {} in ancestor {}", 
                             consensus.length(), chr_name, ancestor_block.ancestor_name);
            }
        }
        
        spdlog::info("Ancestor sequence reconstruction completed");
    }
    
    // 将祖先数据写入HAL
    void writeAncestorDataToHAL(hal::AlignmentPtr alignment,
                               const std::vector<AncestorBlock>& ancestor_blocks) {
        spdlog::info("Writing ancestor data to HAL file...");
        
        // 收集所有祖先的序列数据
        std::map<std::string, std::map<std::string, std::string>> ancestor_data;  // ancestor_name -> chr_name -> sequence
        
        for (const auto& block : ancestor_blocks) {
            for (const auto& [chr_name, sequence] : block.ancestor_sequences) {
                ancestor_data[block.ancestor_name][chr_name] = sequence;
            }
        }
        
        // 为每个祖先节点写入数据
        for (const auto& [ancestor_name, chr_sequences] : ancestor_data) {
            hal::Genome* ancestor_genome = alignment->openGenome(ancestor_name);
            if (!ancestor_genome) {
                spdlog::warn("Cannot open ancestor genome: {}", ancestor_name);
                continue;
            }
            
            // 设置祖先基因组的维度
            std::vector<hal::Sequence::Info> dimensions;
            for (const auto& [chr_name, sequence] : chr_sequences) {
                dimensions.emplace_back(chr_name, sequence.length(), 0, 0);
            }
            
            if (!dimensions.empty()) {
                ancestor_genome->setDimensions(dimensions);
                
                // 连接所有染色体序列
                std::string full_sequence;
                for (const auto& [chr_name, sequence] : chr_sequences) {
                    full_sequence += sequence;
                }
                
                ancestor_genome->setString(full_sequence);
                
                spdlog::info("  Set ancestor '{}' with {} chromosomes, {} bp total", 
                            ancestor_name, dimensions.size(), full_sequence.length());
            }
            
            alignment->closeGenome(ancestor_genome);
        }
        
        spdlog::info("Ancestor data written to HAL successfully");
    }

    // 验证和记录HAL文件状态
    void validateHalFile(hal::AlignmentPtr alignment, const FilePath& hal_path) {
        spdlog::info("Validating HAL file structure...");
        spdlog::info("Total genomes created: {}", alignment->getNumGenomes());
        
        try {
            std::string root_name = alignment->getRootName();
            if (!root_name.empty()) {
                spdlog::info("Root genome: {}", root_name);
            }
            
            std::string tree = alignment->getNewickTree();
            if (!tree.empty()) {
                spdlog::info("Phylogenetic tree: {}", tree);
            }
        } catch (const std::exception& e) {
            spdlog::warn("Error accessing HAL structure: {}", e.what());
        }
        
        // 验证文件
        if (std::filesystem::exists(hal_path) && std::filesystem::file_size(hal_path) > 0) {
            spdlog::info("HAL file created successfully: {} bytes", std::filesystem::file_size(hal_path));
        } else {
            throw std::runtime_error("HAL file was not created or is empty!");
        }
    }

} // namespace hal_converter
} // namespace RaMesh 