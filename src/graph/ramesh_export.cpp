// RaMeshExport.cpp
// ------------------------------------------------------------
// 统一导出实现：emitMafBlock + 三个外层接口
// ------------------------------------------------------------
#include "ramesh.h"
#include "hal_converter.h"
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <unordered_map>
#include <exception>
#include <map>
#include <vector>
#include <spdlog/spdlog.h>
#include <optional>
#include <variant>
#include <functional>
#include "../../submodule/hal/api/inc/hal.h"

// ============================================================
// emitMafBlock —— 所有导出函数共享的“写一个 MAF 块”实现
// ============================================================
static bool emitMafBlock(std::ostream& os,
    const RaMesh::BlockPtr& blk,
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
    bool only_primary,
    bool pairwise_mode,
    bool allow_reverse,
    const SpeciesName* first_sp = nullptr,
    bool ensure_forward = true)
{
    if (!blk) return false;

    // ---------- 1. 收集 segment ----------
    struct Rec { SpeciesName sp; ChrName chr; RaMesh::SegPtr seg; };
    std::vector<Rec> recs;
    recs.reserve(blk->anchors.size());
    bool have_reverse = false;
    for (auto& [sp_chr, seg] : blk->anchors) {
        if (only_primary && !seg->isPrimary()) continue;
        if (!allow_reverse && seg->strand == Strand::REVERSE) have_reverse = true;
        recs.push_back({ sp_chr.first, sp_chr.second, seg });
    }
    if (recs.size() < 2 || (!allow_reverse && have_reverse)) return false;

    // ---------- 2. 决定首行 ----------
    if (first_sp) {
        auto it_first = std::find_if(recs.begin(), recs.end(),
            [&](auto& r) { return r.sp == *first_sp; });
        if (it_first == recs.end()) {
            auto it_ref = std::find_if(recs.begin(), recs.end(),
                [&](auto& r) { return r.chr == blk->ref_chr; });
            if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);
        }
        else {
            if (it_first != recs.begin()) std::swap(*recs.begin(), *it_first);
        }
        
    }
    else {
        auto it_ref = std::find_if(recs.begin(), recs.end(),
            [&](auto& r) { return r.chr == blk->ref_chr; });
        if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);
    }

    // ---------- 3. lambda: fetchSeq / fetchLen ----------
    auto fetchSeq = [](const SeqPro::ManagerVariant& mv,
        const ChrName& chr, Coord_t b, Coord_t l) {
            return std::visit([&](auto& p) {
                using T = std::decay_t<decltype(p)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>)
                    return p->getSubSequence(chr, b, l);
                else
                    return p->getOriginalManager().getSubSequence(chr, b, l);
                }, mv);
        };
    auto fetchLen = [](const SeqPro::ManagerVariant& mv,
        const ChrName& chr) {
            return std::visit([&](auto& p) {
                using T = std::decay_t<decltype(p)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>)
                    return p->getSequenceLength(chr);
                else
                    return p->getOriginalManager().getSequenceLength(chr);
                }, mv);
        };

    // ---------- 4. 准备原始序列 ----------
    std::unordered_map<ChrName, std::string> seqs;
    std::unordered_map<ChrName, Cigar_t>     cigars;
    for (auto& r : recs) {
        auto mit = seq_mgrs.find(r.sp);
        if (mit == seq_mgrs.end()) { seqs.clear(); break; }
        std::string raw = fetchSeq(*mit->second, r.chr, r.seg->start, r.seg->length);
        if (r.seg->strand == Strand::REVERSE) reverseComplement(raw);
        ChrName key = pairwise_mode ? r.chr : r.sp + "." + r.chr;
        seqs.emplace(key, std::move(raw));
        cigars.emplace(key, r.seg->cigar);
    }
    if (seqs.empty()) return false;

    // ---------- 5. 归并对齐 ----------
    ChrName ref_key = blk->ref_chr;
    if (pairwise_mode) {
        ref_key = blk->ref_chr;
    }
    else {
        auto& rr = *std::find_if(recs.begin(), recs.end(),
            [&](auto& r) { return r.chr == blk->ref_chr; });
        ref_key = rr.sp + "." + rr.chr;
    }
    try {
        mergeAlignmentByRef(ref_key, seqs, cigars);
    }
    catch (const std::exception& e) {
        spdlog::warn("mergeAlignmentByRef failed: {}", e.what());
        return false;
    }

    // ---------- 6. 判断是否整体翻转 ----------
    bool need_flip = false;
    if (first_sp && ensure_forward)
        need_flip = (recs.front().seg->strand == Strand::REVERSE);
    std::unordered_map<ChrName, std::string> flipped;
    
    if (need_flip) {
        //auto rc = [](auto& s) { std::string t(s.rbegin(), s.rend()); for (char& c : t) { switch (c) { case 'A':c = 'T';break;case 'a':c = 't';break;case 'T':c = 'A';break;case 't':c = 'a';break;case 'C':c = 'G';break;case 'c':c = 'g';break;case 'G':c = 'C';break;case 'g':c = 'c';break;default:; } } return t; };
        for (auto& kv : seqs) reverseComplement(kv.second);
    }
    const auto* view = &seqs;

    // ---------- 7. 写块 ----------
    os << "a score=0\n";
    for (size_t i = 0;i < recs.size();++i) {
        auto& r = recs[i]; auto& mgr = *seq_mgrs.at(r.sp);
        uint64_t chr_len = fetchLen(mgr, r.chr);
        bool orig_rev = (r.seg->strand == Strand::REVERSE);
        bool final_rev = need_flip ? !orig_rev : orig_rev;
        uint64_t start = need_flip ? (chr_len - r.seg->start - r.seg->length) : r.seg->start;
        ChrName key = pairwise_mode ? r.chr : r.sp + "." + r.chr;
        os << "s " << std::left << std::setw(20) << key << std::right << std::setw(12) << start << std::setw(12) << r.seg->length
            << ' ' << (final_rev ? '-' : '+') << std::setw(12) << chr_len << ' ' << view->at(key) << "\n";
    }
    os << "\n";
    return true;
}

// ============================================================
// RaMeshMultiGenomeGraph 导出接口
// ============================================================
namespace RaMesh {

    void RaMeshMultiGenomeGraph::exportToMaf(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        if (!maf_path.parent_path().empty()) fs::create_directories(maf_path.parent_path());
        std::ofstream ofs(maf_path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open: " + maf_path.string());
        ofs << "##maf version=1 scoring=none\n";
        for (auto& wblk : blocks) emitMafBlock(ofs, wblk.lock(), seq_mgrs, only_primary, pairwise_mode, true, nullptr, true);
    }

    void RaMeshMultiGenomeGraph::exportToMafWithoutReverse(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        if (!maf_path.parent_path().empty()) fs::create_directories(maf_path.parent_path());
        std::ofstream ofs(maf_path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open: " + maf_path.string());
        ofs << "##maf version=1 scoring=none\n";
        for (auto& wblk : blocks) emitMafBlock(ofs, wblk.lock(), seq_mgrs, only_primary, pairwise_mode, false, nullptr, true);
    }

    void RaMeshMultiGenomeGraph::exportToMultipleMaf(
        const std::vector<std::pair<SpeciesName, FilePath>>& outs,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        struct Out { SpeciesName sp; std::ofstream ofs; }; std::vector<Out> dests;
        dests.reserve(outs.size());
        for (auto& pr : outs) {
            if (!pr.second.parent_path().empty()) fs::create_directories(pr.second.parent_path());
            std::ofstream f(pr.second, std::ios::binary | std::ios::trunc);
            if (!f) throw std::runtime_error("Cannot open: " + pr.second.string());
            f << "##maf version=1 scoring=none\n";
            dests.push_back({ pr.first,std::move(f) });
        }
        for (auto& wblk : blocks) {
            auto blk = wblk.lock(); if (!blk) continue;
            for (auto& out : dests) {
                emitMafBlock(out.ofs, blk, seq_mgrs, only_primary, pairwise_mode, true, &out.sp, true);
            }
        }
    }

    void RaMeshMultiGenomeGraph::exportToHal(const FilePath& hal_path,
                                            const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
                                            const std::string& newick_tree,
                                            bool only_primary) const
    {
        spdlog::info("Starting HAL export to: {}", hal_path.string());
        spdlog::info("Export configuration: {} species, tree={}, only_primary={}",
                     seqpro_managers.size(), newick_tree.empty() ? "none" : "provided", only_primary);

        // ========================================
        // 基础验证和准备
        // ========================================
        if (seqpro_managers.empty()) {
            throw std::runtime_error("No sequence managers provided for HAL export");
        }

        if (blocks.empty()) {
            spdlog::warn("No alignment blocks found - creating HAL with genomes only");
        }

        try {
            // ========================================
            // 第一阶段：基础结构建立
            // ========================================
            spdlog::info("Phase 1: Establishing basic HAL structure...");

            // 1.1 创建HAL文件和alignment对象
            std::filesystem::path abs_hal_path = std::filesystem::absolute(hal_path);
            if (abs_hal_path.has_parent_path() && !abs_hal_path.parent_path().empty()) {
                std::filesystem::create_directories(abs_hal_path.parent_path());
            }

            hal::AlignmentPtr alignment = hal::openHalAlignment(abs_hal_path.string(), nullptr, hal::CREATE_ACCESS);
            if (!alignment) {
                throw std::runtime_error("Failed to create HAL alignment instance");
            }
            spdlog::info("  HAL file created: {}", abs_hal_path.string());

            // 1.2 解析系统发育树并识别祖先节点
            std::vector<hal_converter::AncestorNode> ancestor_nodes;
            NewickParser parser; // 创建NewickParser对象用于获取分支长度

            spdlog::info("  Parsing phylogenetic tree...");
            ancestor_nodes = hal_converter::parsePhylogeneticTree(newick_tree, seqpro_managers);
            spdlog::info("  Found {} ancestor nodes from tree", ancestor_nodes.size());

            // 如果有Newick树，解析它以获取分支长度信息
            if (!newick_tree.empty()) {
                try {
                    parser = NewickParser(newick_tree);
                } catch (const std::exception& e) {
                    spdlog::warn("Failed to parse Newick tree for branch lengths: {}", e.what());
                }
            }

            // 1.3 创建叶节点基因组（祖先基因组稍后创建）
            if (!ancestor_nodes.empty()) {
                hal_converter::createLeafGenomes(alignment, ancestor_nodes, parser);
                spdlog::info("  Created leaf genomes");
            }
            spdlog::info("Phase 1 completed successfully");

            // ========================================
            // 第二阶段：祖先序列重建
            // ========================================
            spdlog::info("Phase 2: Reconstructing ancestor sequences...");

            // 声明祖先序列变量，供后续阶段使用
            std::map<std::string, std::string> ancestor_sequences;

            if (!ancestor_nodes.empty() && !blocks.empty()) {
                // ========================================
                // 子阶段 2.1：确定重建顺序与参考叶
                // ========================================
                spdlog::info("  Phase 2.1: Determining reconstruction order and reference leaves...");

                auto reconstruction_plan = planAncestorReconstruction(ancestor_nodes, parser);
                spdlog::info("  Planned reconstruction for {} ancestors", reconstruction_plan.size());

                for (const auto& [ancestor_name, ref_leaf] : reconstruction_plan) {
                    spdlog::debug("    Ancestor '{}' -> Reference leaf '{}'", ancestor_name, ref_leaf);
                }

                spdlog::info("  Phase 2.1 completed successfully");

                // ========================================
                // 子阶段 2.2：执行祖先节点重建
                // ========================================
                spdlog::info("  Phase 2.2: Executing ancestor sequence reconstruction...");

                auto ancestor_reconstruction_data = hal_converter::reconstructAncestorSequences(
                    reconstruction_plan, ancestor_nodes, seqpro_managers, const_cast<RaMeshMultiGenomeGraph&>(*this));

                spdlog::info("  Phase 2.2 completed successfully");

                // ========================================
                // 子阶段 2.3：构建祖先序列（投票法）
                // ========================================
                spdlog::info("  Phase 2.3: Building ancestor sequences using voting method...");

                ancestor_sequences = hal_converter::buildAllAncestorSequencesByVoting(
                    ancestor_reconstruction_data, ancestor_nodes, seqpro_managers);
                spdlog::info("  Phase 2.3 completed successfully");

                // ========================================
                // 子阶段 2.4：用正确维度创建祖先基因组
                // ========================================
                spdlog::info("  Phase 2.4: Creating ancestor genomes with correct dimensions...");
                hal_converter::createAncestorGenomesWithCorrectDimensions(
                    alignment, ancestor_sequences, ancestor_nodes, seqpro_managers);
                spdlog::info("  Phase 2.4 completed successfully");

                // ========================================
                // 子阶段 2.5：应用系统发育树结构
                // ========================================
                spdlog::info("  Phase 2.5: Applying phylogenetic tree structure to HAL alignment...");
                hal_converter::applyPhylogeneticTree(alignment, parser);
                spdlog::info("  Phase 2.5 completed successfully");

            } else {
                spdlog::info("  Skipping ancestor reconstruction (no ancestors or blocks)");
            }
            spdlog::info("Phase 2 completed successfully");

            // ========================================
            // 第三阶段：遍历blocks分析segment映射关系
            // ========================================
            spdlog::info("Phase 3: Analyzing blocks and building segment mappings...");

            size_t processed_blocks = 0;
            std::map<std::string, size_t> species_block_count;
            std::map<std::string, std::map<std::string, size_t>> top_segment_counts;
            std::map<std::string, std::map<std::string, size_t>> bottom_segment_counts;
            // std::map<std::string, std::vector<SegmentMappingInfo>> segment_mappings; // TODO: 定义SegmentMappingInfo结构

            if (!blocks.empty()) {
                spdlog::info("  Analyzing {} alignment blocks for segment relationships...", blocks.size());

                // 3.1 遍历所有blocks，分析segment关系和映射
                for (const auto& weak_block : blocks) {
                    auto block = weak_block.lock();
                    if (!block) continue;

                    // TODO: 分析这个block中的segment关系
                    // - 确定每个基因组在这个block中的segment
                    // - 计算parent-child segment映射关系
                    // - 统计每个基因组需要的TopSegment和BottomSegment数量
                    // analyzeBlockSegmentRelationships(block, segment_mappings, top_segment_counts, bottom_segment_counts);

                    for (const auto& [species_chr, segment] : block->anchors) {
                        species_block_count[species_chr.first]++;
                    }
                    processed_blocks++;
                }

                spdlog::info("  Processed {} alignment blocks", processed_blocks);
                for (const auto& [species, count] : species_block_count) {
                    spdlog::info("    Species {}: {} segments", species, count);
                }

                // 3.2 建立完整的parent-child segment映射
                spdlog::info("  Building parent-child segment mappings...");
                // TODO: buildParentChildSegmentMappings(ancestor_nodes, segment_mappings);

                spdlog::info("  TODO: Implement block segment relationship analysis");
                spdlog::info("  TODO: Calculate TopSegment and BottomSegment counts for each genome");
                spdlog::info("  TODO: Build parent-child segment mapping tables");
            } else {
                spdlog::info("  No alignment blocks to analyze");
            }
            spdlog::info("Phase 3 completed successfully");

            // ========================================
            // 第四阶段：更新HAL维度为正确维度
            // ========================================
            spdlog::info("Phase 4: Updating HAL dimensions to correct values...");

            if (!blocks.empty() && !ancestor_nodes.empty()) {
                spdlog::info("  Updating genome segment dimensions based on block analysis...");

                // 4.1 更新叶节点基因组的segment维度
                spdlog::info("  4.1: Updating leaf genome segment dimensions...");
                // TODO: 遍历所有叶节点基因组，使用updateTopDimensions()和updateBottomDimensions()
                // for (const auto& [genome_name, segment_counts] : top_segment_counts) {
                //     if (isLeafGenome(genome_name, ancestor_nodes)) {
                //         hal_converter::updateGenomeTopSegmentDimensions(alignment, genome_name, segment_counts);
                //         hal_converter::updateGenomeBottomSegmentDimensions(alignment, genome_name, bottom_segment_counts[genome_name]);
                //     }
                // }
                spdlog::info("    TODO: Implement leaf genome segment dimension updates");

                // 4.2 更新祖先基因组的segment维度
                spdlog::info("  4.2: Updating ancestor genome segment dimensions...");
                // TODO: 遍历所有祖先基因组，使用updateTopDimensions()和updateBottomDimensions()
                // for (const auto& ancestor : ancestor_nodes) {
                //     hal_converter::updateGenomeTopSegmentDimensions(alignment, ancestor.node_name, top_segment_counts[ancestor.node_name]);
                //     hal_converter::updateGenomeBottomSegmentDimensions(alignment, ancestor.node_name, bottom_segment_counts[ancestor.node_name]);
                // }
                spdlog::info("    TODO: Implement ancestor genome segment dimension updates");

                // 4.3 验证维度更新结果
                spdlog::info("  4.3: Validating dimension updates...");
                // TODO: 验证所有基因组的segment维度是否正确设置
                // hal_converter::validateSegmentDimensions(alignment, top_segment_counts, bottom_segment_counts);
                spdlog::info("    TODO: Validate segment dimension updates");

                spdlog::info("  TODO: Implement updateGenomeTopSegmentDimensions() using genome->updateTopDimensions()");
                spdlog::info("  TODO: Implement updateGenomeBottomSegmentDimensions() using genome->updateBottomDimensions()");
                spdlog::info("  TODO: Calculate correct segment counts for each chromosome based on block analysis");
            } else {
                spdlog::info("  Skipping dimension updates (no blocks or ancestors)");
            }
            spdlog::info("Phase 4 completed successfully");

            // ========================================
            // 第五阶段：填充segment数据
            // ========================================
            spdlog::info("Phase 5: Populating segment data...");

            if (!blocks.empty() && !ancestor_nodes.empty()) {
                spdlog::info("  Creating and populating HAL segments...");

                // 5.1 创建并填充TopSegment和BottomSegment
                spdlog::info("  5.1: Creating and populating TopSegment and BottomSegment objects...");
                // TODO: 遍历所有基因组，创建实际的segment对象
                // for (const auto& [genome_name, segments] : segment_mappings) {
                //     hal_converter::populateGenomeSegments(alignment, genome_name, segments);
                // }
                spdlog::info("    TODO: Create TopSegment objects for each genome");
                spdlog::info("    TODO: Create BottomSegment objects for each genome");
                spdlog::info("    TODO: Set segment coordinates and lengths");

                // 5.2 建立segment之间的parent-child关系
                spdlog::info("  5.2: Establishing segment parent-child relationships...");
                // TODO: 根据之前分析的映射关系，建立segment之间的连接
                // hal_converter::establishSegmentRelationships(alignment, parent_child_mappings);
                spdlog::info("    TODO: Link TopSegments to their parent BottomSegments");
                spdlog::info("    TODO: Link BottomSegments to their child TopSegments");
                spdlog::info("    TODO: Set segment indices and array positions");

                // 5.3 填充比对信息和坐标
                spdlog::info("  5.3: Populating alignment information and coordinates...");
                // TODO: 将RaMesh的比对块信息转换为HAL segment的比对信息
                // hal_converter::populateSegmentAlignments(alignment, blocks);
                spdlog::info("    TODO: Convert RaMesh blocks to HAL segment alignments");
                spdlog::info("    TODO: Set segment start/end coordinates");
                spdlog::info("    TODO: Handle gaps and insertions");

                spdlog::info("  TODO: Implement populateGenomeSegments()");
                spdlog::info("  TODO: Implement establishSegmentRelationships()");
                spdlog::info("  TODO: Implement populateSegmentAlignments()");
            } else {
                spdlog::info("  Skipping segment data population (no blocks or ancestors)");
            }
            spdlog::info("Phase 5 completed successfully");

            // ========================================
            // 第六阶段：验证和完成
            // ========================================
            spdlog::info("Phase 6: Final validation and completion...");

            // 6.1 验证HAL结构完整性
            spdlog::info("  6.1: Validating HAL structure integrity...");
            try {
                // TODO: 实现完整的HAL结构验证
                // hal_converter::validateHalStructure(alignment);
                spdlog::info("    TODO: Validate segment continuity");
                spdlog::info("    TODO: Verify parent-child relationships");
                spdlog::info("    TODO: Ensure coordinate consistency");
                spdlog::info("    TODO: Check DNA sequence integrity");

                hal_converter::validateHalFile(alignment, hal_path);
                spdlog::info("    HAL file validation passed");
            } catch (const std::exception& e) {
                spdlog::error("    HAL file validation failed: {}", e.what());
                throw;
            }

            // 6.2 生成导出统计报告
            spdlog::info("  6.2: Generating export statistics...");
            // TODO: 生成详细的导出统计信息
            // hal_converter::generateExportStatistics(alignment, processed_blocks);
            spdlog::info("    TODO: Generate genome statistics");
            spdlog::info("    TODO: Generate segment statistics");
            spdlog::info("    TODO: Generate alignment coverage statistics");

            // 6.3 关闭文件和清理
            spdlog::info("  6.3: Closing files and cleanup...");
            alignment->close();
            spdlog::info("    HAL file closed successfully");

            spdlog::info("Phase 6 completed successfully");

            spdlog::info("HAL export completed successfully: {}", abs_hal_path.string());

        } catch (const hal_exception& e) {
            spdlog::error("HAL API error during export: {}", e.what());
            spdlog::error("This usually indicates a problem with HAL file structure or API usage");
            throw;
        } catch (const std::filesystem::filesystem_error& e) {
            spdlog::error("Filesystem error during export: {}", e.what());
            spdlog::error("Path: {}", e.path1().string());
            throw;
        } catch (const std::runtime_error& e) {
            spdlog::error("Runtime error during export: {}", e.what());
            throw;
        } catch (const std::exception& e) {
            spdlog::error("Unexpected error during HAL export: {}", e.what());
            spdlog::error("Export failed - HAL file may be incomplete or corrupted");
            throw;
        }
    }
} // namespace RaMesh
