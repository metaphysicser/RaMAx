#include "hal_converter.h"
#include "align.h"
#include <spdlog/spdlog.h>

namespace RaMesh {
namespace hal_converter {

    // ========================================
    // 系统发育树解析和处理
    // ========================================

    std::vector<AncestorNode> parsePhylogeneticTree(
        const std::string& newick_tree,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        spdlog::info("Parsing phylogenetic tree to identify ancestor nodes...");

        std::vector<AncestorNode> ancestor_nodes;

        if (newick_tree.empty()) {
            spdlog::info("Empty Newick tree provided, creating star topology");

            // 创建星形拓扑：一个根节点连接所有叶节点
            AncestorNode root_ancestor;
            root_ancestor.node_name = "ancestor";  // 默认根节点名称
            root_ancestor.is_generated_root = true;
            root_ancestor.tree_depth = 0;
            root_ancestor.branch_length = 0.0;

            // 收集所有叶节点作为根的后代
            for (const auto& [species_name, _] : seqpro_managers) {
                root_ancestor.descendant_leaves.push_back(species_name);
                root_ancestor.children_names.push_back(species_name);
            }

            ancestor_nodes.push_back(root_ancestor);
            spdlog::info("Created star topology with root 'ancestor' and {} leaves",
                        root_ancestor.descendant_leaves.size());

            return ancestor_nodes;
        }

        // 解析Newick树
        try {
            NewickParser parser(newick_tree);

            // 验证叶节点名称
            auto [is_valid, error_msg] = validateLeafNames(parser, seqpro_managers);
            if (!is_valid) {
                spdlog::warn("Leaf name validation failed: {}", error_msg);
                spdlog::warn("Proceeding with available species...");
            }

            // 检查并确保有根节点
            NewickParser mutable_parser = parser;  // 创建可修改的副本
            bool added_root = ensureRootNode(mutable_parser, seqpro_managers);
            if (added_root) {
                spdlog::info("Added artificial root node 'ancestor'");
            }

            // 提取祖先节点信息
            ancestor_nodes = extractAncestorNodes(mutable_parser, seqpro_managers);

            spdlog::info("Found {} ancestor nodes from tree", ancestor_nodes.size());
            for (const auto& ancestor : ancestor_nodes) {
                spdlog::debug("  Ancestor '{}' with {} descendants (depth: {})",
                             ancestor.node_name, ancestor.descendant_leaves.size(), ancestor.tree_depth);
            }

        } catch (const std::exception& e) {
            spdlog::error("Failed to parse Newick tree: {}", e.what());
            spdlog::info("Falling back to star topology");

            // 回退到星形拓扑
            AncestorNode root_ancestor;
            root_ancestor.node_name = "ancestor";
            root_ancestor.is_generated_root = true;
            root_ancestor.tree_depth = 0;

            for (const auto& [species_name, _] : seqpro_managers) {
                root_ancestor.descendant_leaves.push_back(species_name);
                root_ancestor.children_names.push_back(species_name);
            }

            ancestor_nodes.push_back(root_ancestor);
        }

        return ancestor_nodes;
    }

    bool ensureRootNode(NewickParser& parser,
                       const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        const auto& nodes = parser.getNodes();
        if (nodes.empty()) {
            return false;
        }

        // 找到根节点（father == -1的节点）
        int root_id = -1;
        for (const auto& node : nodes) {
            if (node.father == -1) {
                root_id = node.id;
                break;
            }
        }

        if (root_id == -1) {
            spdlog::error("No root node found in parsed tree");
            return false;
        }

        const auto& root_node = nodes[root_id];

        // 检查根节点是否是叶节点
        if (root_node.isLeaf) {
            spdlog::info("Root is a leaf node, adding artificial root");

            // 创建新的根节点
            NewickTreeNode new_root;
            new_root.id = parser.currentIndex_++;
            new_root.name = "ancestor";
            new_root.father = -1;
            new_root.isLeaf = false;
            new_root.branchLength = 0.0;
            new_root.leftChild = root_id;
            new_root.rightChild = -1;  // 只有一个子节点

            // 更新原根节点
            parser.nodes_[root_id].father = new_root.id;
            parser.nodes_[root_id].branchLength = 1.0;  // 设置默认分支长度

            // 添加新根节点
            parser.nodes_.push_back(new_root);

            return true;
        }

        // 检查根节点是否有名称，如果没有则设置默认名称
        if (root_node.name.empty()) {
            parser.nodes_[root_id].name = "ancestor";
            spdlog::info("Set root node name to 'ancestor'");
            return true;
        }

        return false;
    }

    std::vector<AncestorNode> extractAncestorNodes(
        const NewickParser& parser,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        std::vector<AncestorNode> ancestor_nodes;
        const auto& nodes = parser.getNodes();

        if (nodes.empty()) {
            return ancestor_nodes;
        }

        // 为每个内部节点创建AncestorNode
        for (const auto& node : nodes) {
            if (!node.isLeaf) {  // 只处理内部节点
                AncestorNode ancestor;
                ancestor.node_name = node.name.empty() ? ("internal_" + std::to_string(node.id)) : node.name;
                ancestor.branch_length = node.branchLength;
                ancestor.is_generated_root = (node.name == "ancestor" && node.father == -1);

                // 计算树深度
                ancestor.tree_depth = 0;
                int current_id = node.id;
                while (current_id != -1) {
                    bool found = false;
                    for (const auto& n : nodes) {
                        if (n.id == current_id && n.father != -1) {
                            ancestor.tree_depth++;
                            current_id = n.father;
                            found = true;
                            break;
                        }
                    }
                    if (!found) break;
                }

                // 设置父节点名称
                if (node.father != -1) {
                    for (const auto& n : nodes) {
                        if (n.id == node.father) {
                            ancestor.parent_name = n.name.empty() ? ("internal_" + std::to_string(n.id)) : n.name;
                            break;
                        }
                    }
                }

                // 收集所有后代叶节点
                std::function<void(int)> collectLeaves = [&](int node_id) {
                    for (const auto& n : nodes) {
                        if (n.father == node_id) {
                            if (n.isLeaf) {
                                // 只添加在seqpro_managers中存在的叶节点
                                if (seqpro_managers.find(n.name) != seqpro_managers.end()) {
                                    ancestor.descendant_leaves.push_back(n.name);
                                }
                            } else {
                                collectLeaves(n.id);
                            }

                            // 添加到直接子节点列表
                            std::string child_name = n.name.empty() ? ("internal_" + std::to_string(n.id)) : n.name;
                            ancestor.children_names.push_back(child_name);
                        }
                    }
                };

                collectLeaves(node.id);

                // 只有当祖先节点有后代叶节点时才添加
                if (!ancestor.descendant_leaves.empty()) {
                    ancestor_nodes.push_back(ancestor);
                }
            }
        }

        return ancestor_nodes;
    }

    std::pair<bool, std::string> validateLeafNames(
        const NewickParser& parser,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        std::vector<std::string> tree_leaves = parser.getLeafNames();
        std::set<std::string> manager_species;

        for (const auto& [species_name, _] : seqpro_managers) {
            manager_species.insert(species_name);
        }

        std::vector<std::string> missing_in_tree;
        std::vector<std::string> missing_in_managers;

        // 检查树中的叶节点是否都在managers中
        for (const auto& leaf : tree_leaves) {
            if (manager_species.find(leaf) == manager_species.end()) {
                missing_in_managers.push_back(leaf);
            }
        }

        // 检查managers中的物种是否都在树中
        std::set<std::string> tree_leaf_set(tree_leaves.begin(), tree_leaves.end());
        for (const auto& species : manager_species) {
            if (tree_leaf_set.find(species) == tree_leaf_set.end()) {
                missing_in_tree.push_back(species);
            }
        }

        if (missing_in_tree.empty() && missing_in_managers.empty()) {
            return {true, "All leaf names match"};
        }

        std::string error_msg = "Leaf name mismatches: ";
        if (!missing_in_tree.empty()) {
            error_msg += "Missing in tree: ";
            for (const auto& name : missing_in_tree) {
                error_msg += name + " ";
            }
        }
        if (!missing_in_managers.empty()) {
            error_msg += "Missing in managers: ";
            for (const auto& name : missing_in_managers) {
                error_msg += name + " ";
            }
        }

        return {false, error_msg};
    }

    // ========================================
    // HAL基础结构创建
    // ========================================
    void setupGenomeSequences(
        hal::AlignmentPtr alignment,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers) {

        spdlog::info("Setting up sequence data for all genomes...");

        for (const auto& [species_name, seq_mgr] : seqpro_managers) {
            hal::Genome* genome = alignment->openGenome(species_name);
            if (!genome) {
                spdlog::warn("Cannot open genome: {}", species_name);
                continue;
            }

            // 获取序列信息并创建维度
            std::vector<hal::Sequence::Info> sequence_dimensions;

            std::visit([&](const auto& mgr) {
                auto chr_names = mgr->getSequenceNames();
                for (const auto& chr_name : chr_names) {
                    hal_size_t length = mgr->getSequenceLength(chr_name);
                    sequence_dimensions.emplace_back(chr_name, length, 0, 0);
                }
            }, *seq_mgr);

            // 设置基因组维度和DNA数据
            if (!sequence_dimensions.empty()) {
                genome->setDimensions(sequence_dimensions);

                // 添加DNA数据（暂时用N填充）
                // TODO 更换真实DNA数据
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

    void createAncestorGenomes(
        hal::AlignmentPtr alignment,
        const std::vector<AncestorNode>& ancestor_nodes,
        const NewickParser& parser) {

        spdlog::info("Creating all genomes (ancestors and leaves)...");

        // 首先创建根节点
        for (const auto& ancestor : ancestor_nodes) {
            if (ancestor.parent_name.empty()) {
                // 这是根节点
                hal::Genome* existing = alignment->openGenome(ancestor.node_name);
                if (existing) {
                    alignment->closeGenome(existing);
                    spdlog::debug("Root genome '{}' already exists", ancestor.node_name);
                    continue;
                }

                hal::Genome* root_genome = alignment->addRootGenome(ancestor.node_name);
                if (!root_genome) {
                    throw std::runtime_error("Failed to create root genome: " + ancestor.node_name);
                }
                spdlog::info("Created root genome: {}", ancestor.node_name);
                break; // 只应该有一个根节点
            }
        }

        // 然后创建内部节点（祖先）
        for (const auto& ancestor : ancestor_nodes) {
            if (!ancestor.parent_name.empty()) {
                // 跳过已经存在的基因组
                hal::Genome* existing = alignment->openGenome(ancestor.node_name);
                if (existing) {
                    alignment->closeGenome(existing);
                    spdlog::debug("Ancestor genome '{}' already exists", ancestor.node_name);
                    continue;
                }

                // 这是内部节点
                hal::Genome* ancestor_genome = alignment->addLeafGenome(ancestor.node_name, ancestor.parent_name, ancestor.branch_length);
                if (!ancestor_genome) {
                    spdlog::error("Failed to create ancestor genome: {}", ancestor.node_name);
                    continue;
                }

                spdlog::info("Created ancestor genome: {}", ancestor.node_name);
            }
        }

        // 最后创建叶节点基因组
        for (const auto& ancestor : ancestor_nodes) {
            for (const auto& child_name : ancestor.children_names) {
                // 检查这个child是否是叶节点（不在ancestor_nodes中）
                bool is_leaf = true;
                for (const auto& other_ancestor : ancestor_nodes) {
                    if (other_ancestor.node_name == child_name) {
                        is_leaf = false;
                        break;
                    }
                }

                if (is_leaf) {
                    // 跳过已经存在的基因组
                    hal::Genome* existing = alignment->openGenome(child_name);
                    if (existing) {
                        alignment->closeGenome(existing);
                        spdlog::debug("Leaf genome '{}' already exists", child_name);
                        continue;
                    }

                    // 创建叶节点基因组
                    // 从NewickParser中获取叶节点的正确分支长度
                    double leaf_branch_length = 1.0; // 默认分支长度
                    const auto& nodes = parser.getNodes();
                    for (const auto& node : nodes) {
                        if (node.isLeaf && node.name == child_name) {
                            leaf_branch_length = node.branchLength;
                            break;
                        }
                    }
                    hal::Genome* leaf_genome = alignment->addLeafGenome(child_name, ancestor.node_name, leaf_branch_length);
                    if (!leaf_genome) {
                        spdlog::error("Failed to create leaf genome: {}", child_name);
                        continue;
                    }

                    spdlog::info("Created leaf genome: {}", child_name);
                }
            }
        }
    }

    void applyPhylogeneticTree(
        hal::AlignmentPtr alignment,
        const NewickParser& parser) {

        spdlog::info("Applying phylogenetic tree structure to HAL alignment...");

        // 注意：NewickParser没有直接获取原始字符串的方法
        // 这个函数需要在调用时传递原始的newick字符串
        // 这里我们先跳过树的应用，因为我们已经在parsePhylogeneticTree中处理了树结构

        spdlog::info("Tree structure already processed during parsing phase");
    }

    // ========================================
    // 验证和工具函数
    // ========================================

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

        spdlog::info("HAL file validation completed: {}", hal_path.string());
    }

    std::string generateRootName(const std::vector<std::string>& leaf_names) {
        if (leaf_names.empty()) {
            return "ancestor";
        }

        // 生成一个不与任何叶节点冲突的根节点名称
        std::set<std::string> leaf_set(leaf_names.begin(), leaf_names.end());

        std::string base_name = "ancestor";
        std::string root_name = base_name;
        int counter = 1;

        while (leaf_set.find(root_name) != leaf_set.end()) {
            root_name = base_name + "_" + std::to_string(counter);
            counter++;
        }

        return root_name;
    }

} // namespace hal_converter
} // namespace RaMesh