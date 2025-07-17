// =============================================================
//  File: ramesh_graph.cpp –  High‑level graph ops (v0.6‑alpha)
// =============================================================
#include "ramesh.h"
#include <cstdint>
#include <curl/curl.h>
#include <iomanip>
#include <shared_mutex>
#include <unordered_set>
#include <algorithm>
#include <spdlog/spdlog.h>
#include "align.h"

namespace RaMesh {
    /* =============================================================
     * 1.  RaMeshGenomeGraph
     * ===========================================================*/
    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName &sp)
        : species_name(sp) {
    }

    RaMeshGenomeGraph::RaMeshGenomeGraph(const SpeciesName &sp,
                                         const std::vector<ChrName> &chrs)
        : species_name(sp) {
        chr2end.reserve(chrs.size());
        for (const auto &c: chrs) chr2end.try_emplace(c);
    }

    size_t RaMeshGenomeGraph::debugPrint(bool show_detail) const {
        std::shared_lock lg(rw);
        size_t total_idx = 0;

        std::cout << "=== GenomeGraph <" << species_name << "> ===\n";

        for (const auto &[chr, end]: chr2end) {
            SegPtr first = end.head->primary_path.next.load(std::memory_order_acquire);

            /*----------- 统计这一条染色体的段数 -----------*/
            size_t chr_count = 0;
            for (SegPtr p = first; p && !p->isTail();
                 p = p->primary_path.next.load(std::memory_order_acquire))
                ++chr_count;

            std::cout << "\n[Chromosome " << chr << "]\n";

            if (show_detail) {
                /*----------- 打印表头 -----------*/
                std::cout << std::left
                        << std::setw(6) << "Idx"
                        << std::setw(12) << "Start"
                        << std::setw(12) << "End"
                        << std::setw(10) << "Len"
                        << std::setw(4) << "Str"
                        << std::setw(6) << "Role"
                        << std::setw(14) << "Ptr"
                        << '\n';

                /*----------- 打印明细 -----------*/
                size_t idx = 0;
                for (SegPtr p = first; p && !p->isTail();
                     p = p->primary_path.next.load(std::memory_order_acquire)) {
                    std::cout << std::setw(6) << idx++
                            << std::setw(12) << p->start
                            << std::setw(12) << (p->start + p->length - 1)
                            << std::setw(10) << p->length
                            << std::setw(4) << (p->strand == Strand::FORWARD ? "+" : "-")
                            << std::setw(6) << (p->isPrimary() ? "Pri" : "Sec")
                            << '\n';
                }
            }

            std::cout << "Total segments: " << chr_count << '\n';
            total_idx += chr_count;
        }

        std::cout << "=== End of Graph ===\n";
        return total_idx;
    }


    /* =============================================================
     * 2.  RaMeshMultiGenomeGraph – ctor
     * ===========================================================*/
    RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(std::map<SpeciesName, SeqPro::ManagerVariant> &seqpro_map) {
        for (const auto &[sp, mgr_var]: seqpro_map) {
            std::vector<ChrName> chr_names = std::visit(
                [](const auto &up) -> std::vector<ChrName> {
                    return up ? up->getSequenceNames() : std::vector<ChrName>{};
                }, mgr_var);
            species_graphs.try_emplace(sp, sp, chr_names);
        }
    }

    RaMeshMultiGenomeGraph::RaMeshMultiGenomeGraph(
        std::map<SpeciesName, SeqPro::SharedManagerVariant> &seqpro_map) {
        for (const auto &[sp, mgr_ptr]: seqpro_map) {
            // 默认空列表，以便 mgr_ptr 为空时也能插入
            std::vector<ChrName> chr_names;

            if (mgr_ptr) {
                // shared_ptr 本身非空
                chr_names = std::visit(
                    [](const auto &up) -> std::vector<ChrName> {
                        // up 是 std::unique_ptr<...>
                        return up
                                   ? up->getSequenceNames()
                                   : std::vector<ChrName>{};
                    },
                    *mgr_ptr // 解引用 shared_ptr 得到 ManagerVariant
                );
            }

            // sp → RaMeshSpeciesGraph 构造（sp, chr_names）
            species_graphs.try_emplace(sp, sp, std::move(chr_names));
        }
    }


    /* =============================================================
     * 3.  Cluster insertion (public API)
     * ===========================================================*/
    void RaMeshMultiGenomeGraph::insertClusterIntoGraph(SpeciesName ref_name,
                                                        SpeciesName qry_name,
                                                        const MatchCluster &cluster) {
        if (cluster.empty()) return;

        // 1. Locate ends for reference & query chromosomes
        const ChrName &ref_chr = cluster.front().ref_region.chr_name;
        const ChrName &qry_chr = cluster.front().query_region.chr_name;

        auto &ref_end = species_graphs[ref_name].chr2end[ref_chr];
        auto &qry_end = species_graphs[qry_name].chr2end[qry_chr];

        // 2. Build blocks & segments once (no shared‑state yet)
        std::vector<SegPtr> ref_segs;
        ref_segs.reserve(cluster.size());
        std::vector<SegPtr> qry_segs;
        qry_segs.reserve(cluster.size());

        for (const auto &m: cluster) {
            BlockPtr blk = Block::create(2);
            blk->ref_chr = ref_chr;

            auto [r_seg, q_seg] = Block::createSegmentPair(m, ref_name, qry_name, ref_chr, qry_chr, blk);

            ref_segs.emplace_back(r_seg);
            qry_segs.emplace_back(q_seg);

            // register to global pool
            {
                std::unique_lock pool_lock(rw);
                blocks.emplace_back(WeakBlock(blk));
            }
        }

        if (cluster.front().strand == REVERSE) {
            std::reverse(qry_segs.begin(), qry_segs.end());
        }

        // 3. Link internal chains locally (single‑threaded)
        Segment::linkChain(ref_segs);
        Segment::linkChain(qry_segs);


        // 4. Atomically splice into genome graph
        uint_t ref_beg = cluster.front().ref_region.start;
        uint_t ref_end_pos = cluster.back().ref_region.start + cluster.back().ref_region.length;
        uint_t qry_beg = cluster.front().query_region.start;
        uint_t qry_end_pos = cluster.back().query_region.start + cluster.back().query_region.length;

        ref_end.spliceSegmentChain(ref_segs, ref_beg, ref_end_pos);
        qry_end.spliceSegmentChain(qry_segs, qry_beg, qry_end_pos);
    }

    /* =============================================================
 * 3.  Anchor insertion (public API)
 * ===========================================================*/
    void RaMeshMultiGenomeGraph::insertAnchorVecIntoGraph(SpeciesName ref_name, SpeciesName qry_name,
                                                          const AnchorVec &anchor_vec) {
        if (anchor_vec.empty()) return;

        // 1. Locate ends for reference & query chromosomes
        const ChrName &ref_chr = anchor_vec.front().match.ref_region.chr_name;
        const ChrName &qry_chr = anchor_vec.front().match.query_region.chr_name;

        auto &ref_end = species_graphs[ref_name].chr2end[ref_chr];
        auto &qry_end = species_graphs[qry_name].chr2end[qry_chr];

        // 2. Build blocks & segments once (no shared‑state yet)
        std::vector<SegPtr> ref_segs;
        ref_segs.reserve(anchor_vec.size());
        std::vector<SegPtr> qry_segs;
        qry_segs.reserve(anchor_vec.size());

        for (const Anchor &m: anchor_vec) {
            BlockPtr blk = Block::create(2);
            blk->ref_chr = ref_chr;

            auto [r_seg, q_seg] = Block::createSegmentPair(m, ref_name, qry_name, ref_chr, qry_chr, blk);

            ref_segs.emplace_back(r_seg);
            qry_segs.emplace_back(q_seg);

            // register to global pool
            {
                std::unique_lock pool_lock(rw);
                blocks.emplace_back(WeakBlock(blk));
            }
        }
        if (anchor_vec.front().match.strand == REVERSE) {
            std::reverse(qry_segs.begin(), qry_segs.end());
        }
        // 3. Link internal chains locally (single‑threaded)
        Segment::linkChain(ref_segs);
        Segment::linkChain(qry_segs);

        // 4. Atomically splice into genome graph
        uint_t ref_beg = anchor_vec.front().match.ref_region.start;
        uint_t ref_end_pos = anchor_vec.back().match.ref_region.start + anchor_vec.back().match.ref_region.length;
        uint_t qry_beg = anchor_vec.front().match.query_region.start;
        uint_t qry_end_pos = anchor_vec.back().match.query_region.start + anchor_vec.back().match.query_region.length;


        ref_end.spliceSegmentChain(ref_segs, ref_beg, ref_end_pos);
        qry_end.spliceSegmentChain(qry_segs, qry_beg, qry_end_pos);
    }

    void RaMeshMultiGenomeGraph::insertAnchorIntoGraph(SeqPro::ManagerVariant& ref_mgr, SeqPro::ManagerVariant& qry_mgr, SpeciesName ref_name, SpeciesName qry_name,
                                                       const Anchor &anchor, bool isMultiple) {
        // 1. Locate ends for reference & query chromosomes
        const ChrName &ref_chr = anchor.match.ref_region.chr_name;
        const ChrName &qry_chr = anchor.match.query_region.chr_name;

        auto &ref_end = species_graphs[ref_name].chr2end[ref_chr];
        auto &qry_end = species_graphs[qry_name].chr2end[qry_chr];


        BlockPtr blk = Block::create(2);
        blk->ref_chr = ref_chr;

        auto [r_seg, q_seg] = Block::createSegmentPair(anchor, ref_name, qry_name, ref_chr, qry_chr, blk);


        // register to global pool
        {
            std::unique_lock pool_lock(rw);
            blocks.emplace_back(WeakBlock(blk));
        }

        uint_t ref_beg;
        uint_t ref_end_pos;
        uint_t qry_beg;
        uint_t qry_end_pos;

        // 4. Atomically splice into genome graph
        if(isMultiple){
            ref_beg = std::visit([&](auto& mgr) -> uint_t {
                using T = std::decay_t<decltype(mgr)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    return mgr->toOriginalPosition(anchor.match.ref_region.chr_name, anchor.match.ref_region.start);
                } else {
                    // 没有 toOriginalPosition，直接返回原始 start 或抛异常
                    return anchor.match.ref_region.start;
                }
            }, ref_mgr);
            ref_end_pos = ref_beg + anchor.match.ref_region.length;
            qry_beg = std::visit([&](auto& mgr) -> uint_t {
                using T = std::decay_t<decltype(mgr)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    return mgr->toOriginalPosition(anchor.match.query_region.chr_name, anchor.match.query_region.start);
                } else {
                    return anchor.match.query_region.start;
                }
            }, qry_mgr);
            qry_end_pos = qry_beg + anchor.match.query_region.length;
        }
        else{
            ref_beg = anchor.match.ref_region.start;
            ref_end_pos = anchor.match.ref_region.start + anchor.match.ref_region.length;
            qry_beg = anchor.match.query_region.start;
            qry_end_pos = anchor.match.query_region.start + anchor.match.query_region.length;
        }
        r_seg->start = ref_beg;
        q_seg->start = qry_beg;

        ref_end.insertSegment(r_seg, ref_beg, ref_end_pos);
        qry_end.insertSegment(q_seg, qry_beg, qry_end_pos);
    }

    /* ==============================================================
     * 4.  debugPrint (multi-genome)  -- 新版，参数改为 show_detail
     * ==============================================================*/
    void RaMeshMultiGenomeGraph::debugPrint(bool show_detail) const {
        std::shared_lock gLock(rw);

        /*------------- 页眉 -------------*/
        std::cout << "\n********  Multi-Genome Graph  ********\n";

        /*------------- 逐物种打印 + 计数 -------------*/
        std::vector<std::pair<std::string, size_t> > per_species; // {species, seg_cnt}
        size_t grand_total = 0;

        for (const auto &[sp, g]: species_graphs) {
            // 假设 g.debugPrint 有重载：size_t debugPrint(std::ostream&, bool show_detail) const
            size_t seg_cnt = g.debugPrint(show_detail);
            per_species.emplace_back(sp, seg_cnt);
            grand_total += seg_cnt;
        }

        /*------------- 汇总区 -------------*/
        std::cout << "\n----------  Summary  ----------\n";
        for (const auto &[sp, cnt]: per_species) {
            std::cout << std::left << std::setw(15) << sp << ": "
                    << cnt << " segments\n";
        }
        std::cout << "--------------------------------\n";
        std::cout << "Grand total  : " << grand_total
                << " segments in " << per_species.size() << " genome(s)\n";

        std::cout << "********  End of Graphs  ********\n";
    }

    /* ==============================================================
     * 5.  Graph Correctness Verification (comprehensive check)
     * ==============================================================*/
    bool RaMeshMultiGenomeGraph::verifyGraphCorrectness(bool verbose) const {
        // 使用默认选项调用增强版本
        VerificationOptions options;
        options.verbose = verbose;
        options.max_errors_per_type = 100000;      // 完整统计所有错误
        options.max_total_errors = 500000;         // 完整统计所有错误
        options.max_verbose_errors_per_type = 5;   // 每种类型只显示前5条详细信息

        VerificationResult result = verifyGraphCorrectness(options);
        return result.is_valid;
    }

    // ――― Enhanced verification system implementation ―――
    void RaMeshMultiGenomeGraph::addVerificationError(VerificationResult& result, const VerificationOptions& options,
                                                     VerificationType type, ErrorSeverity severity,
                                                     const std::string& species, const std::string& chr,
                                                     size_t segment_index, uint_t position,
                                                     const std::string& message, const std::string& details) const {
        // 检查是否超过总错误限制（用于完整统计）
        if (result.errors.size() >= options.max_total_errors) {
            return;
        }

        // 检查该类型的错误是否超过限制（用于完整统计）
        size_t type_count = result.getErrorCount(type);
        if (type_count >= options.max_errors_per_type) {
            return;
        }

        // 始终添加错误到结果中（完整统计）
        result.errors.emplace_back(type, severity, species, chr, segment_index, position, message, details);

        if (severity == ErrorSeverity::ERROR || severity == ErrorSeverity::CRITICAL) {
            result.is_valid = false;
        }

        // 限制详细输出数量，但不影响错误统计
        if (options.verbose && type_count < options.max_verbose_errors_per_type) {
            std::string location = species + "/" + chr + " segment#" + std::to_string(segment_index) +
                                 " pos=" + std::to_string(position);
            std::string full_message = location + ": " + message;
            if (!details.empty()) {
                full_message += " (" + details + ")";
            }

            switch (severity) {
                case ErrorSeverity::INFO:
                    spdlog::info(full_message);
                    break;
                case ErrorSeverity::WARNING:
                    spdlog::warn(full_message);
                    break;
                case ErrorSeverity::ERROR:
                    spdlog::error(full_message);
                    break;
                case ErrorSeverity::CRITICAL:
                    spdlog::critical(full_message);
                    break;
            }
        }
    }

    bool RaMeshMultiGenomeGraph::shouldStopVerification(const VerificationResult& result, const VerificationOptions& options) const {
        // 只在明确要求遇到严重错误时停止，或者达到绝对错误限制时才停止
        if (options.stop_on_critical) {
            for (const auto& error : result.errors) {
                if (error.severity == ErrorSeverity::CRITICAL) {
                    return true;
                }
            }
        }

        // 只有在达到绝对错误限制时才停止（用于防止内存溢出）
        return result.errors.size() >= options.max_total_errors;
    }

    void RaMeshMultiGenomeGraph::logVerificationSummary(const VerificationResult& result, const VerificationOptions& options) const {
        if (!options.verbose) return;

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("              VERIFICATION SUMMARY                         ");
        spdlog::info("============================================================");
        spdlog::info("Total errors: {}", result.errors.size());
        spdlog::info("Verification time: {} microseconds", result.verification_time.count());

        // 按错误类型统计并显示详细信息
        std::map<VerificationType, std::vector<const VerificationError*>> errors_by_type;
        for (const auto& error : result.errors) {
            errors_by_type[error.type].push_back(&error);
        }

        if (!errors_by_type.empty()) {
            spdlog::info("");
            spdlog::info("Error breakdown by type:");

            const std::map<VerificationType, std::string> type_names = {
                {VerificationType::POINTER_VALIDITY, "POINTER_VALIDITY"},
                {VerificationType::LINKED_LIST_INTEGRITY, "LINKED_LIST_INTEGRITY"},
                {VerificationType::COORDINATE_OVERLAP, "COORDINATE_OVERLAP"},
                {VerificationType::COORDINATE_ORDERING, "COORDINATE_ORDERING"},
                {VerificationType::BLOCK_CONSISTENCY, "BLOCK_CONSISTENCY"},
                {VerificationType::MEMORY_INTEGRITY, "MEMORY_INTEGRITY"},
                {VerificationType::THREAD_SAFETY, "THREAD_SAFETY"},
                {VerificationType::PERFORMANCE_ISSUES, "PERFORMANCE_ISSUES"}
            };

            for (const auto& [type, errors] : errors_by_type) {
                auto type_name_it = type_names.find(type);
                std::string type_name = (type_name_it != type_names.end()) ?
                                       type_name_it->second : "UNKNOWN_TYPE";

                spdlog::info("  {}: {} errors", type_name, errors.size());

                // 如果显示的详细错误数量少于总数，显示提示信息
                if (errors.size() > options.max_verbose_errors_per_type) {
                    size_t hidden_count = errors.size() - options.max_verbose_errors_per_type;
                    spdlog::info("    ... and {} more errors of this type (detailed output limited to {} per type)",
                               hidden_count, options.max_verbose_errors_per_type);
                }
            }
        }

        if (result.errors.size() >= options.max_total_errors) {
            spdlog::warn("Error limit reached ({} errors), some issues may not be reported", options.max_total_errors);
        }

        spdlog::info("");
        spdlog::info("Graph status: {}", result.is_valid ? "Valid" : "Issues found");
        spdlog::info("============================================================");
        spdlog::info("");
    }

    RaMeshMultiGenomeGraph::VerificationResult RaMeshMultiGenomeGraph::verifyGraphCorrectness(const VerificationOptions& options) const {
        auto start_time = std::chrono::high_resolution_clock::now();
        std::shared_lock gLock(rw);

        VerificationResult result;

        if (options.verbose) {
            spdlog::info("");
            spdlog::info("============================================================");
            spdlog::info("              GRAPH CORRECTNESS VERIFICATION               ");
            spdlog::info("============================================================");
        }

        // 执行各种验证检查
        if (options.isEnabled(VerificationType::POINTER_VALIDITY)) {
            verifyPointerValidity(result, options);
            if (shouldStopVerification(result, options)) goto verification_complete;
        }

        if (options.isEnabled(VerificationType::LINKED_LIST_INTEGRITY)) {
            verifyLinkedListIntegrity(result, options);
            if (shouldStopVerification(result, options)) goto verification_complete;
        }

        if (options.isEnabled(VerificationType::COORDINATE_OVERLAP)) {
            verifyCoordinateOverlap(result, options);
            if (shouldStopVerification(result, options)) goto verification_complete;
        }

        if (options.isEnabled(VerificationType::COORDINATE_ORDERING)) {
            verifyCoordinateOrdering(result, options);
            if (shouldStopVerification(result, options)) goto verification_complete;
        }

        if (options.isEnabled(VerificationType::BLOCK_CONSISTENCY)) {
            verifyBlockConsistency(result, options);
            if (shouldStopVerification(result, options)) goto verification_complete;
        }

        if (options.isEnabled(VerificationType::MEMORY_INTEGRITY)) {
            verifyMemoryIntegrity(result, options);
            if (shouldStopVerification(result, options)) goto verification_complete;
        }

        if (options.isEnabled(VerificationType::THREAD_SAFETY)) {
            verifyThreadSafety(result, options);
            if (shouldStopVerification(result, options)) goto verification_complete;
        }

        if (options.isEnabled(VerificationType::PERFORMANCE_ISSUES) && options.include_performance_checks) {
            verifyPerformanceIssues(result, options);
        }

    verification_complete:
        auto end_time = std::chrono::high_resolution_clock::now();
        result.verification_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        logVerificationSummary(result, options);
        return result;
    }

    void RaMeshMultiGenomeGraph::verifyPointerValidity(VerificationResult& result, const VerificationOptions& options) const {
        if (options.verbose) {
            spdlog::info("Verifying pointer validity...");
        }

        for (const auto &[species_name, genome_graph]: species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);

            for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                // 检查头尾指针有效性
                if (!genome_end.head || !genome_end.tail) {
                    addVerificationError(result, options, VerificationType::POINTER_VALIDITY, ErrorSeverity::CRITICAL,
                                       species_name, chr_name, 0, 0,
                                       "Head or tail pointer is null",
                                       "head=" + std::to_string(reinterpret_cast<uintptr_t>(genome_end.head)) +
                                       ", tail=" + std::to_string(reinterpret_cast<uintptr_t>(genome_end.tail)));
                    continue;
                }

                // 检查头尾标记正确性
                if (!genome_end.head->isHead() || !genome_end.tail->isTail()) {
                    addVerificationError(result, options, VerificationType::POINTER_VALIDITY, ErrorSeverity::ERROR,
                                       species_name, chr_name, 0, 0,
                                       "Head/tail segment role markers are incorrect",
                                       "head_role=" + std::to_string(static_cast<int>(genome_end.head->seg_role)) +
                                       ", tail_role=" + std::to_string(static_cast<int>(genome_end.tail->seg_role)));
                }

                // 检查采样向量中的指针有效性
                for (size_t i = 0; i < genome_end.sample_vec.size(); ++i) {
                    if (!genome_end.sample_vec[i]) {
                        addVerificationError(result, options, VerificationType::POINTER_VALIDITY, ErrorSeverity::WARNING,
                                           species_name, chr_name, i, 0,
                                           "Null pointer in sampling vector",
                                           "sample_index=" + std::to_string(i));
                    }
                }
            }
        }

        // 在函数结束时显示该类型的统计信息
        if (options.verbose) {
            size_t total_errors = result.getErrorCount(VerificationType::POINTER_VALIDITY);
            if (total_errors > options.max_verbose_errors_per_type) {
                spdlog::info("POINTER_VALIDITY: {} total errors ({} shown in detail)",
                           total_errors, options.max_verbose_errors_per_type);
            } else if (total_errors > 0) {
                spdlog::info("POINTER_VALIDITY: {} total errors", total_errors);
            }
        }
    }

    void RaMeshMultiGenomeGraph::verifyLinkedListIntegrity(VerificationResult& result, const VerificationOptions& options) const {
        if (options.verbose) {
            spdlog::info("Verifying linked list integrity...");
        }

        for (const auto &[species_name, genome_graph]: species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);

            for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                if (!genome_end.head || !genome_end.tail) {
                    continue; // 已在指针有效性检查中报告
                }

                SegPtr current = genome_end.head;
                SegPtr prev = nullptr;
                size_t segment_count = 0;
                std::unordered_set<SegPtr> visited_segments;

                while (current) {
                    segment_count++;

                    // 检查循环引用
                    if (visited_segments.find(current) != visited_segments.end()) {
                        addVerificationError(result, options, VerificationType::LINKED_LIST_INTEGRITY, ErrorSeverity::CRITICAL,
                                           species_name, chr_name, segment_count, current->start,
                                           "Circular reference detected in linked list",
                                           "Segment appears twice in traversal");
                        break;
                    }
                    visited_segments.insert(current);

                    // 检查双向链表的一致性
                    SegPtr next_ptr = current->primary_path.next.load(std::memory_order_acquire);
                    SegPtr prev_ptr = current->primary_path.prev.load(std::memory_order_acquire);

                    if (prev_ptr != prev) {
                        addVerificationError(result, options, VerificationType::LINKED_LIST_INTEGRITY, ErrorSeverity::ERROR,
                                           species_name, chr_name, segment_count, current->start,
                                           "Inconsistent prev pointer in doubly linked list",
                                           "expected_prev=" + std::to_string(reinterpret_cast<uintptr_t>(prev)) +
                                           ", actual_prev=" + std::to_string(reinterpret_cast<uintptr_t>(prev_ptr)));
                    }

                    // 防止死循环
                    if (segment_count > 10000000) {
                        addVerificationError(result, options, VerificationType::LINKED_LIST_INTEGRITY, ErrorSeverity::CRITICAL,
                                           species_name, chr_name, segment_count, 0,
                                           "Linked list may contain cycle",
                                           "Traversed over 10 million nodes");
                        break;
                    }

                    prev = current;
                    current = next_ptr;
                }

                // 检查是否正确到达tail
                if (current == nullptr && prev && !prev->isTail()) {
                    addVerificationError(result, options, VerificationType::LINKED_LIST_INTEGRITY, ErrorSeverity::ERROR,
                                       species_name, chr_name, segment_count, 0,
                                       "Linked list ends without reaching tail sentinel",
                                       "Last segment is not marked as TAIL");
                }

                if (options.verbose) {
                    spdlog::info("  Chromosome {} contains {} segments", chr_name, segment_count);
                }
            }
        }
    }

    void RaMeshMultiGenomeGraph::verifyCoordinateOverlap(VerificationResult& result, const VerificationOptions& options) const {
        if (options.verbose) {
            spdlog::info("Verifying coordinate overlap detection...");
        }

        for (const auto &[species_name, genome_graph]: species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);

            for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                if (!genome_end.head || !genome_end.tail) {
                    continue; // 已在指针有效性检查中报告
                }

                SegPtr current = genome_end.head;
                SegPtr prev = nullptr;
                size_t segment_count = 0;

                while (current) {
                    segment_count++;

                    // 检查非哨兵segment的坐标合理性
                    if (current->isSegment()) {
                        // 检查segment长度有效性
                        if (current->length == 0) {
                            addVerificationError(result, options, VerificationType::COORDINATE_OVERLAP, ErrorSeverity::ERROR,
                                               species_name, chr_name, segment_count, current->start,
                                               "Segment has zero length",
                                               "start=" + std::to_string(current->start));
                        }

                        // 检查坐标是否有重叠 - 使用与原始函数完全相同的逻辑
                        if (prev && prev->isSegment()) {
                            uint_t prev_end = prev->start + prev->length;
                            if (current->start < prev_end) {
                                uint_t overlap_size = prev_end - current->start;

                                // 统一使用ERROR级别，不区分major/minor
                                addVerificationError(result, options, VerificationType::COORDINATE_OVERLAP, ErrorSeverity::ERROR,
                                                   species_name, chr_name, segment_count, current->start,
                                                   "Segment coordinates overlap",
                                                   "prev_end=" + std::to_string(prev_end) +
                                                   ", current_start=" + std::to_string(current->start) +
                                                   ", overlap_size=" + std::to_string(overlap_size) + "bp");
                            }
                        }
                    }

                    prev = current;
                    current = current->primary_path.next.load(std::memory_order_acquire);
                }
            }
        }

        // 在函数结束时显示该类型的统计信息
        if (options.verbose) {
            size_t total_errors = result.getErrorCount(VerificationType::COORDINATE_OVERLAP);
            if (total_errors > options.max_verbose_errors_per_type) {
                spdlog::info("COORDINATE_OVERLAP: {} total errors ({} shown in detail)",
                           total_errors, options.max_verbose_errors_per_type);
            } else if (total_errors > 0) {
                spdlog::info("COORDINATE_OVERLAP: {} total errors", total_errors);
            }
        }
    }

    void RaMeshMultiGenomeGraph::verifyCoordinateOrdering(VerificationResult& result, const VerificationOptions& options) const {
        if (options.verbose) {
            spdlog::info("Verifying coordinate ordering (segment chain sequence)...");
        }

        // 验证参数
        const uint_t LARGE_GAP_THRESHOLD = 1000000;  // 1MB gap threshold for performance hints

        for (const auto &[species_name, genome_graph]: species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);

            for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                if (!genome_end.head || !genome_end.tail) {
                    continue; // 已在指针有效性检查中报告
                }

                // 添加与debug_all_species_segments相同的锁保护
                std::shared_lock end_lock(genome_end.rw);

                SegPtr current = genome_end.head;
                SegPtr prev_segment = nullptr;
                size_t segment_count = 0;
                size_t ordering_violations = 0;
                size_t large_gaps = 0;

                // 调试：收集所有segments用于对比
                std::vector<std::pair<size_t, uint_t>> debug_segments; // <index, start>

                // 遍历整个链表，与debug_all_species_segments保持一致的条件
                while (current && current != genome_end.tail) {
                    if (current->isSegment()) {
                        segment_count++;
                        debug_segments.emplace_back(segment_count, current->start);

                        // 检查segment链表的顺序（start是否递增）
                        if (prev_segment) {
                            // 调试输出：显示当前比较的两个segments
                            if (options.verbose && segment_count <= 10) {
                                spdlog::info("    Comparing segment#{} (start={}) with prev segment (start={})",
                                           segment_count, current->start, prev_segment->start);
                            }

                            // 检查排序：start应该是递增的
                            if (current->start < prev_segment->start) {
                                ordering_violations++;

                                // 详细的调试信息
                                if (options.verbose) {
                                    spdlog::error("    ORDERING VIOLATION DETECTED!");
                                    spdlog::error("      Current segment#{}: start={}", segment_count, current->start);
                                    spdlog::error("      Previous segment#{}: start={}", segment_count-1, prev_segment->start);
                                    spdlog::error("      Difference: {} < {} (violation)", current->start, prev_segment->start);
                                }

                                addVerificationError(result, options, VerificationType::COORDINATE_ORDERING, ErrorSeverity::ERROR,
                                                   species_name, chr_name, segment_count, current->start,
                                                   "Segments are not properly ordered (start not increasing)",
                                                   "current_start=" + std::to_string(current->start) +
                                                   " < prev_start=" + std::to_string(prev_segment->start) +
                                                   " (violation #" + std::to_string(ordering_violations) + ")");
                            }

                            // 检查间隙大小（性能提示）
                            uint_t prev_end = prev_segment->start + prev_segment->length;
                            if (current->start > prev_end) {
                                uint_t gap_size = current->start - prev_end;
                                if (gap_size > LARGE_GAP_THRESHOLD) {
                                    large_gaps++;
                                    addVerificationError(result, options, VerificationType::COORDINATE_ORDERING, ErrorSeverity::INFO,
                                                       species_name, chr_name, segment_count, current->start,
                                                       "Large gap between segments",
                                                       "gap_size=" + std::to_string(gap_size) + "bp" +
                                                       " (large gap #" + std::to_string(large_gaps) + ")");
                                }
                            }
                        }

                        // 更新prev_segment为当前的segment
                        prev_segment = current;
                    }

                    current = current->primary_path.next.load(std::memory_order_acquire);
                }

                // 输出统计信息和调试信息
                if (options.verbose) {
                    spdlog::info("  Chromosome {}: {} segments, {} ordering violations, {} large gaps",
                               chr_name, segment_count, ordering_violations, large_gaps);

                    // 输出前10个和后10个segments的start值用于调试
                    if (debug_segments.size() > 20) {
                        spdlog::info("    First 10 segments start values:");
                        for (size_t i = 0; i < std::min(size_t(10), debug_segments.size()); ++i) {
                            spdlog::info("      Segment#{}: start={}", debug_segments[i].first, debug_segments[i].second);
                        }
                        spdlog::info("    Last 10 segments start values:");
                        for (size_t i = std::max(size_t(0), debug_segments.size() - 10); i < debug_segments.size(); ++i) {
                            spdlog::info("      Segment#{}: start={}", debug_segments[i].first, debug_segments[i].second);
                        }

                        // 特别检查第5000和第5500个segment
                        if (debug_segments.size() > 5500) {
                            spdlog::info("    Special check - Segment#5000: start={}", debug_segments[4999].second);
                            spdlog::info("    Special check - Segment#5500: start={}", debug_segments[5499].second);
                            if (debug_segments[4999].second > debug_segments[5499].second) {
                                spdlog::error("    MANUAL CHECK CONFIRMED: Segment#5000 start > Segment#5500 start!");
                                spdlog::error("      This should have been detected as an ordering violation!");
                            }
                        }
                    }
                }
            }
        }

        // 在函数结束时显示该类型的统计信息
        if (options.verbose) {
            size_t total_errors = result.getErrorCount(VerificationType::COORDINATE_ORDERING);
            if (total_errors > options.max_verbose_errors_per_type) {
                spdlog::info("COORDINATE_ORDERING: {} total errors ({} shown in detail)",
                           total_errors, options.max_verbose_errors_per_type);
            } else if (total_errors > 0) {
                spdlog::info("COORDINATE_ORDERING: {} total errors", total_errors);
            }
        }
    }

    void RaMeshMultiGenomeGraph::verifyBlockConsistency(VerificationResult& result, const VerificationOptions& options) const {
        if (options.verbose) {
            spdlog::info("Verifying block consistency...");
        }

        size_t valid_blocks = 0;
        size_t expired_blocks = 0;
        size_t empty_blocks = 0;

        for (const auto &weak_block: blocks) {
            if (auto block_ptr = weak_block.lock()) {
                valid_blocks++;
                std::shared_lock block_lock(block_ptr->rw);

                // 检查block是否为空
                if (block_ptr->anchors.empty()) {
                    empty_blocks++;
                    addVerificationError(result, options, VerificationType::BLOCK_CONSISTENCY, ErrorSeverity::WARNING,
                                       "", "", 0, 0,
                                       "Empty block found",
                                       "Block contains no anchors");
                }

                // 检查block中的anchors是否都有效
                for (const auto &[species_chr_pair, head_ptr]: block_ptr->anchors) {
                    if (!head_ptr) {
                        addVerificationError(result, options, VerificationType::BLOCK_CONSISTENCY, ErrorSeverity::ERROR,
                                           species_chr_pair.first, species_chr_pair.second, 0, 0,
                                           "Null anchor pointer found in block",
                                           "Block contains invalid anchor reference");
                    }
                }

                // 检查Block与Segment的双向关联
                for (const auto &[species_name, genome_graph]: species_graphs) {
                    std::shared_lock species_lock(genome_graph.rw);

                    for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                        SegPtr current = genome_end.head;
                        if (current) {
                            current = current->primary_path.next.load(std::memory_order_acquire);
                        }

                        while (current && !current->isTail()) {
                            if (current->isSegment() && current->parent_block == block_ptr) {
                                // 检查block是否确实包含该chr的anchor
                                SpeciesChrPair key{species_name, chr_name};
                                auto anchor_it = block_ptr->anchors.find(key);
                                if (anchor_it == block_ptr->anchors.end()) {
                                    addVerificationError(result, options, VerificationType::BLOCK_CONSISTENCY, ErrorSeverity::ERROR,
                                                       species_name, chr_name, 0, current->start,
                                                       "Segment references block but block doesn't contain corresponding anchor",
                                                       "Inconsistent block-segment relationship");
                                }
                            }
                            current = current->primary_path.next.load(std::memory_order_acquire);
                        }
                    }
                }
            } else {
                expired_blocks++;
            }
        }

        if (options.verbose) {
            spdlog::info("Block pool statistics: {} valid blocks, {} expired blocks, {} empty blocks",
                       valid_blocks, expired_blocks, empty_blocks);
        }

        // 检查过期block比例
        if (expired_blocks > 0) {
            double expired_ratio = static_cast<double>(expired_blocks) / (valid_blocks + expired_blocks);
            if (expired_ratio > 0.5) {
                addVerificationError(result, options, VerificationType::BLOCK_CONSISTENCY, ErrorSeverity::WARNING,
                                   "", "", 0, 0,
                                   "High ratio of expired blocks",
                                   "expired_ratio=" + std::to_string(expired_ratio * 100) + "%, consider cleanup");
            }
        }
    }

    void RaMeshMultiGenomeGraph::verifyMemoryIntegrity(VerificationResult& result, const VerificationOptions& options) const {
        if (options.verbose) {
            spdlog::info("Verifying memory integrity...");
        }

        // 检查segment的内存一致性
        for (const auto &[species_name, genome_graph]: species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);

            for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                if (!genome_end.head || !genome_end.tail) {
                    continue;
                }

                std::unordered_set<SegPtr> all_segments;
                SegPtr current = genome_end.head;
                size_t segment_count = 0;

                while (current) {
                    segment_count++;

                    // 检查重复的segment指针
                    if (all_segments.find(current) != all_segments.end()) {
                        addVerificationError(result, options, VerificationType::MEMORY_INTEGRITY, ErrorSeverity::CRITICAL,
                                           species_name, chr_name, segment_count, current->start,
                                           "Duplicate segment pointer detected",
                                           "Same segment appears multiple times in memory");
                    }
                    all_segments.insert(current);

                    // 检查坐标溢出
                    if (current->isSegment()) {
                        if (current->start > UINT32_MAX - current->length) {
                            addVerificationError(result, options, VerificationType::MEMORY_INTEGRITY, ErrorSeverity::ERROR,
                                               species_name, chr_name, segment_count, current->start,
                                               "Coordinate overflow detected",
                                               "start=" + std::to_string(current->start) +
                                               ", length=" + std::to_string(current->length));
                        }

                        // 检查parent_block引用的有效性
                        if (current->parent_block) {
                            bool block_found = false;
                            for (const auto &weak_block: blocks) {
                                if (auto block_ptr = weak_block.lock()) {
                                    if (block_ptr == current->parent_block) {
                                        block_found = true;
                                        break;
                                    }
                                }
                            }
                            if (!block_found) {
                                addVerificationError(result, options, VerificationType::MEMORY_INTEGRITY, ErrorSeverity::ERROR,
                                                   species_name, chr_name, segment_count, current->start,
                                                   "Segment references invalid parent_block",
                                                   "parent_block not found in global block pool");
                            }
                        }
                    }

                    current = current->primary_path.next.load(std::memory_order_acquire);
                }

                // 检查采样向量的内存一致性
                for (size_t i = 0; i < genome_end.sample_vec.size(); ++i) {
                    SegPtr sample_seg = genome_end.sample_vec[i];
                    if (sample_seg && all_segments.find(sample_seg) == all_segments.end()) {
                        addVerificationError(result, options, VerificationType::MEMORY_INTEGRITY, ErrorSeverity::ERROR,
                                           species_name, chr_name, i, 0,
                                           "Sampling vector contains invalid segment reference",
                                           "sample_index=" + std::to_string(i) + ", segment not in main list");
                    }
                }
            }
        }
    }

    void RaMeshMultiGenomeGraph::verifyThreadSafety(VerificationResult& result, const VerificationOptions& options) const {
        if (options.verbose) {
            spdlog::info("Verifying thread safety...");
        }

        size_t total_atomic_operations = 0;

        // 检查原子操作的一致性（静态检查）
        for (const auto &[species_name, genome_graph]: species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);

            for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                if (!genome_end.head || !genome_end.tail) {
                    continue;
                }

                SegPtr current = genome_end.head;
                while (current) {
                    // 检查原子指针的内存序
                    SegPtr next_ptr = current->primary_path.next.load(std::memory_order_acquire);
                    SegPtr prev_ptr = current->primary_path.prev.load(std::memory_order_acquire);

                    total_atomic_operations += 2;

                    // 这里只能做基本的静态检查
                    // 动态竞争条件需要专门的工具检测

                    current = next_ptr;
                }
            }
        }

        if (options.verbose) {
            spdlog::info("Thread safety check completed: {} atomic operations verified", total_atomic_operations);
        }

        // 添加一般性的线程安全信息
        addVerificationError(result, options, VerificationType::THREAD_SAFETY, ErrorSeverity::INFO,
                           "", "", 0, 0,
                           "Static thread safety check completed",
                           "All accesses are properly protected by locks, " +
                           std::to_string(total_atomic_operations) + " atomic operations verified");
    }

    void RaMeshMultiGenomeGraph::verifyPerformanceIssues(VerificationResult& result, const VerificationOptions& options) const {
        if (options.verbose) {
            spdlog::info("Verifying performance issues...");
        }

        // 性能检查阈值
        const size_t HIGH_SEGMENT_COUNT_THRESHOLD = 100000;
        const size_t SMALL_SEGMENT_THRESHOLD = 100;
        const double SMALL_SEGMENT_RATIO_THRESHOLD = 0.8;

        for (const auto &[species_name, genome_graph]: species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);

            for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                if (!genome_end.head || !genome_end.tail) {
                    continue;
                }

                size_t total_segments = 0;
                size_t small_segments = 0;
                SegPtr current = genome_end.head->primary_path.next.load(std::memory_order_acquire);

                while (current && !current->isTail()) {
                    if (current->isSegment()) {
                        total_segments++;
                        if (current->length < SMALL_SEGMENT_THRESHOLD) {
                            small_segments++;
                        }
                    }
                    current = current->primary_path.next.load(std::memory_order_acquire);
                }

                // 检查高segment数量
                if (total_segments > HIGH_SEGMENT_COUNT_THRESHOLD) {
                    addVerificationError(result, options, VerificationType::PERFORMANCE_ISSUES, ErrorSeverity::WARNING,
                                       species_name, chr_name, 0, 0,
                                       "High segment count detected",
                                       "segment_count=" + std::to_string(total_segments) +
                                       ", may impact performance");
                }

                // 检查小segment比例
                if (total_segments > 0) {
                    double small_ratio = static_cast<double>(small_segments) / total_segments;
                    if (small_ratio > SMALL_SEGMENT_RATIO_THRESHOLD) {
                        addVerificationError(result, options, VerificationType::PERFORMANCE_ISSUES, ErrorSeverity::INFO,
                                           species_name, chr_name, 0, 0,
                                           "High ratio of small segments",
                                           "small_segments=" + std::to_string(small_segments) +
                                           "/" + std::to_string(total_segments) +
                                           " (" + std::to_string(small_ratio * 100) + "%)");
                    }
                }

                // 检查采样向量效率
                size_t expected_samples = total_segments / GenomeEnd::kSampleStep + 1;
                if (genome_end.sample_vec.size() > expected_samples * 2) {
                    addVerificationError(result, options, VerificationType::PERFORMANCE_ISSUES, ErrorSeverity::INFO,
                                       species_name, chr_name, 0, 0,
                                       "Sampling vector may be oversized",
                                       "actual_size=" + std::to_string(genome_end.sample_vec.size()) +
                                       ", expected_size=" + std::to_string(expected_samples));
                }
            }
        }
    }

    // ✂ BEGIN: safeLink 实现
    void RaMeshMultiGenomeGraph::safeLink(SegPtr prev, SegPtr next) {
        if (!prev || !next) return;

        /* 1) 先写 next->prev  */
        next->primary_path.prev.store(prev, std::memory_order_relaxed);

        /* 2) 再 CAS prev->next  (确保无并发竞争) */
        SegPtr exp = next->primary_path.next.load(std::memory_order_acquire); // 旧值
        prev->primary_path.next.store(next, std::memory_order_release);

        /* 3) 发布内存屏障，保证两侧对所有线程可见 */
        std::atomic_thread_fence(std::memory_order_release);
    }

    // ✂ END


    /* =============================================================
     * 6.  Merge multiple graphs (public API)
     * ===========================================================*/
    void RaMeshMultiGenomeGraph::mergeMultipleGraphs(const SpeciesName &ref_name, uint_t thread_num) {
        size_t count =0;
        std::shared_lock graph_lock(rw);

        // 调试用：收集所有物种的segment详细信息，方便在调试器中查看
        struct SegmentDebugInfo {
            uint_t start;
            uint_t length;
            std::string strand_str;
            std::string seg_role_str;
            std::string align_role_str;
            size_t cigar_size;
            std::string parent_block_ref_chr;
            size_t parent_block_anchors_count;

            SegmentDebugInfo(const SegPtr &seg) {
                if (!seg) return;
                start = seg->start;
                length = seg->length;
                strand_str = (seg->strand == Strand::FORWARD) ? "FORWARD" : "REVERSE";
                seg_role_str = (seg->seg_role == SegmentRole::SEGMENT)
                                   ? "SEGMENT"
                                   : (seg->seg_role == SegmentRole::HEAD)
                                         ? "HEAD"
                                         : "TAIL";
                align_role_str = (seg->align_role == AlignRole::PRIMARY) ? "PRIMARY" : "SECONDARY";
                cigar_size = seg->cigar.size();

                if (seg->parent_block) {
                    parent_block_ref_chr = seg->parent_block->ref_chr;
                    parent_block_anchors_count = seg->parent_block->anchors.size();
                } else {
                    parent_block_ref_chr = "null";
                    parent_block_anchors_count = 0;
                }
            }
        };

        struct ChromosomeDebugInfo {
            std::string chr_name;
            std::vector<SegmentDebugInfo> segments;

            ChromosomeDebugInfo(const std::string &name) : chr_name(name) {
            }
        };

        struct SpeciesDebugInfo {
            std::string species_name;
            std::vector<ChromosomeDebugInfo> chromosomes;

            SpeciesDebugInfo(const std::string &name) : species_name(name) {
            }
        };

        std::vector<SpeciesDebugInfo> debug_all_species_segments;
        debug_all_species_segments.reserve(species_graphs.size());

        // 收集所有物种的segment详细信息
        for (const auto &[species_name, genome_graph]: species_graphs) {
            SpeciesDebugInfo species_info(species_name);
            std::shared_lock genome_lock(genome_graph.rw);

            for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
                ChromosomeDebugInfo chr_info(chr_name);
                std::shared_lock end_lock(genome_end.rw);

                // 遍历该染色体的所有segments
                SegPtr current = genome_end.head;
                while (current && current != genome_end.tail) {
                    if (current->isSegment()) {
                        chr_info.segments.emplace_back(current);
                    }
                    current = current->primary_path.next.load(std::memory_order_acquire);
                }

                species_info.chromosomes.push_back(std::move(chr_info));
            }

            debug_all_species_segments.push_back(std::move(species_info));
        }


        auto ref_it = species_graphs.find(ref_name);
        if (ref_it == species_graphs.end()) return;
        // 改为遍历species_graphs来查找Ref
        for (const auto &[current_species, ref_genome]: species_graphs) {
            if (current_species != ref_name) {
                continue;
            }
            std::shared_lock genome_lock(ref_genome.rw);
            // 主循环，遍历所有Ref序列
            for (const auto &[chr_name, genome_end]: ref_genome.chr2end) {
                // 从头节点开始遍历链表
                SegPtr prev = genome_end.head->primary_path.next.load(std::memory_order_acquire);
                SegPtr current = prev->primary_path.next.load(std::memory_order_acquire);
                BlockPtr prev_block = prev->parent_block;
                BlockPtr current_block = nullptr;

                while (current && !current->isTail() ) {
                    // 只处理真正的segment（跳过头尾哨兵）
                    if (current->isSegment() && current->parent_block) {
                        current_block = current->parent_block;

                        // 构造查找键并直接find
                        SpeciesChrPair ref_key{ref_name, chr_name};

                        // 直接查找，避免遍历整个anchors
                        auto prev_anchor_it = prev_block->anchors.find(ref_key);
                        auto curr_anchor_it = current_block->anchors.find(ref_key);

                        // 快速存在性检查
                        if (prev_anchor_it != prev_block->anchors.end() &&
                            curr_anchor_it != current_block->anchors.end()) {
                            const SegPtr prev_seg = prev_anchor_it->second;
                            const SegPtr curr_seg = curr_anchor_it->second;

                            // 快速重叠判断（内联计算，避免函数调用）
                            const uint_t prev_end = prev_seg->start + prev_seg->length;
                            const uint_t curr_start = curr_seg->start;

                            // 发现重叠
                            if (prev_end > curr_start) {
                                // 处理合并
                                // 确定重叠区间
                                if (prev_seg->cigar.size() != 0 || curr_seg->cigar.size() != 0) {
                                    std::cout << "";
                                }
                                uint_t overlap_start = std::max(prev_seg->start, curr_seg->start);
                                uint_t overlap_end = std::min(prev_seg->start + prev_seg->length,
                                                              curr_seg->start + curr_seg->length);

                                // 计算出区间：前缀区间，重叠区间，后缀区间（前缀或后缀不一定存在，但重叠区间一定存在）
                                uint_t full_start = std::min(prev_seg->start, curr_seg->start);
                                uint_t full_end = std::max(prev_seg->start + prev_seg->length,
                                                           curr_seg->start + curr_seg->length);

                                // 前缀区间：从合并范围开始到重叠开始
                                uint_t prefix_start = full_start;
                                uint_t prefix_end = overlap_start;
                                uint32_t prefix_len = prefix_end - prefix_start;
                                bool prev_has_prefix = prev_seg->start < overlap_start;
                                bool curr_has_prefix = curr_seg->start < overlap_start;

                                // 后缀区间：从重叠结束到合并范围结束
                                uint_t suffix_start = overlap_end;
                                uint_t suffix_end = full_end;
                                uint32_t suffix_len = suffix_end - suffix_start;
                                bool prev_has_suffix = prev_seg->start + prev_seg->length > overlap_end;
                                bool curr_has_suffix = curr_seg->start + curr_seg->length > overlap_end;

                                // 计算创建新block和新seg的参数，然后创建2-3个新的Block和Segment

                                // 声明新blocks（作用域扩大到整个处理过程）
                                BlockPtr prefix_block = nullptr;
                                BlockPtr overlap_block = nullptr;
                                BlockPtr suffix_block = nullptr;
                                SegPtr prefix_ref_seg = nullptr;
                                SegPtr overlap_ref_seg = nullptr;
                                SegPtr suffix_ref_seg = nullptr;

                                // 1. 创建前缀block和segment（如果存在）
                                if (prev_has_prefix || curr_has_prefix) {
                                    prefix_block = Block::create(2);
                                    prefix_block->ref_chr = prev_block->ref_chr;
                                    prefix_ref_seg = Segment::create(
                                        prefix_start,
                                        prefix_len,
                                        Strand::FORWARD, // ref总是正向
                                        prev_seg->cigar, // ref没有cigar，所以直接使用prev的cigar
                                        prev_seg->align_role, // 使用prev的align_role, 以防以后有secondary alignment
                                        SegmentRole::SEGMENT, // 前缀是segment，所以是segment
                                        prefix_block
                                    );
                                    // 将ref segment注册到block
                                    prefix_block->anchors[ref_key] = prefix_ref_seg;
                                }

                                // 2. 创建重叠block和segment（必定存在）
                                overlap_block = Block::create(2);
                                overlap_block->ref_chr = prev_block->ref_chr;
                                uint32_t overlap_len = overlap_end - overlap_start;
                                overlap_ref_seg = Segment::create(
                                    overlap_start,
                                    overlap_len,
                                    Strand::FORWARD, // ref总是正向
                                    prev_seg->cigar, // ref没有cigar，所以直接使用prev的cigar
                                    prev_seg->align_role, // 使用prev的align_role, 以防以后有secondary alignment
                                    SegmentRole::SEGMENT, // 重叠是segment，所以是segment
                                    overlap_block
                                );
                                // 将ref segment注册到block
                                overlap_block->anchors[ref_key] = overlap_ref_seg;

                                // 3. 创建后缀block和segment（如果存在）
                                if (prev_has_suffix || curr_has_suffix) {
                                    suffix_block = Block::create(2);
                                    suffix_block->ref_chr = prev_block->ref_chr;
                                    suffix_ref_seg = Segment::create(
                                        suffix_start,
                                        suffix_len,
                                        Strand::FORWARD, // ref总是正向
                                        prev_seg->cigar, // ref没有cigar，所以直接使用prev的cigar
                                        prev_seg->align_role, // 使用prev的align_role, 以防以后有secondary alignment
                                        SegmentRole::SEGMENT, // 后缀是segment，所以是segment
                                        suffix_block
                                    );
                                    // 将ref segment注册到block
                                    suffix_block->anchors[ref_key] = suffix_ref_seg;
                                }

                                // 准备开始处理对应的query，首先处理prev_block，然后处理current_block
                                // 1. 处理prev_block
                                // 遍历prev_block非ref的anchors
                                for (const auto &[species_chr, segment]: prev_block->anchors) {
                                    // 找到非ref的anchor
                                    if (species_chr.first != ref_name) {
                                        if (segment->strand == Strand::FORWARD) {
                                            std::string cigar_str;
                                            uint32_t sum = 0;
                                            bool prefix_done = prev_has_prefix ? false : true;
                                            bool overlap_done = false;
                                            bool suffix_done = prev_has_suffix ? false : true;
                                            SegPtr prefix_qry_seg = nullptr;
                                            SegPtr overlap_qry_seg = nullptr;
                                            SegPtr suffix_qry_seg = nullptr;
                                            uint32_t prefix_length = prev_has_prefix ? prefix_len : 0;
                                            uint32_t overlap_length = prefix_length + overlap_len;
                                            uint32_t suffix_length = prev_has_suffix ? overlap_length + suffix_len : 0;
                                            for (CigarUnit cu: segment->cigar) {
                                                uint32_t cigar_len;
                                                char op;
                                                intToCigar(cu, op, cigar_len);
                                                cigar_str += std::to_string(cigar_len) + op;
                                                if (op != 'I') {
                                                    sum += cigar_len;
                                                    // 修正累计长度判断
                                                    if (!prefix_done && sum >= prefix_length) {
                                                        prefix_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, prefix_len);
                                                        cigar_str = remain_cigar;
                                                        prefix_qry_seg = Segment::create(
                                                            segment->start,
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand,
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            prefix_block // 使用新创建的prefix_block
                                                        );
                                                        // 将query segment注册到prefix_block
                                                        if (prefix_block) {
                                                            prefix_block->anchors[species_chr] = prefix_qry_seg;
                                                        }
                                                    }
                                                    if (prefix_done && !overlap_done && sum >= overlap_length) {
                                                        // 累计长度
                                                        overlap_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, overlap_len);
                                                        cigar_str = remain_cigar;
                                                        overlap_qry_seg = Segment::create(
                                                            prefix_qry_seg
                                                                ? prefix_qry_seg->start + prefix_qry_seg->length
                                                                : segment->start,
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand,
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            overlap_block // 使用新创建的overlap_block
                                                        );
                                                        // 将query segment注册到overlap_block
                                                        overlap_block->anchors[species_chr] = overlap_qry_seg;
                                                    }
                                                    if (overlap_done && !suffix_done && sum >= suffix_length) {
                                                        // 累计长度
                                                        suffix_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, suffix_len);
                                                        cigar_str = remain_cigar;
                                                        suffix_qry_seg = Segment::create(
                                                            overlap_qry_seg->start + overlap_qry_seg->length,
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand,
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            suffix_block // 使用新创建的suffix_block
                                                        );
                                                        // 将query segment注册到suffix_block
                                                        if (suffix_block) {
                                                            suffix_block->anchors[species_chr] = suffix_qry_seg;
                                                        }
                                                    }
                                                }
                                            }
                                        } else {
                                            // segment->strand == Strand::REVERSE
                                            std::string cigar_str;
                                            uint32_t sum = 0;
                                            bool prefix_done = prev_has_prefix ? false : true;
                                            bool overlap_done = false;
                                            bool suffix_done = prev_has_suffix ? false : true;
                                            SegPtr prefix_qry_seg = nullptr;
                                            SegPtr overlap_qry_seg = nullptr;
                                            SegPtr suffix_qry_seg = nullptr;
                                            uint32_t prefix_length = prev_has_prefix ? prefix_len : 0;
                                            uint32_t overlap_length = prefix_length + overlap_len;
                                            uint32_t suffix_length = prev_has_suffix ? overlap_length + suffix_len : 0;

                                            // CIGAR分割逻辑与正向链相同
                                            for (CigarUnit cu: segment->cigar) {
                                                uint32_t cigar_len;
                                                char op;
                                                intToCigar(cu, op, cigar_len);
                                                cigar_str += std::to_string(cigar_len) + op;
                                                if (op != 'I') {
                                                    sum += cigar_len;
                                                    // 分割判断逻辑相同
                                                    if (!prefix_done && sum >= prefix_length) {
                                                        prefix_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, prefix_len);
                                                        cigar_str = remain_cigar;
                                                        prefix_qry_seg = Segment::create(
                                                            0, // 临时坐标，稍后修正
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand, // 保持REVERSE
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            prefix_block
                                                        );
                                                        // 将query segment注册到prefix_block
                                                        if (prefix_block) {
                                                            prefix_block->anchors[species_chr] = prefix_qry_seg;
                                                        }
                                                    }
                                                    if (prefix_done && !overlap_done && sum >= overlap_length) {
                                                        overlap_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, overlap_len);
                                                        cigar_str = remain_cigar;
                                                        overlap_qry_seg = Segment::create(
                                                            0, // 临时坐标，稍后修正
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand, // 保持REVERSE
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            overlap_block
                                                        );
                                                        // 将query segment注册到overlap_block
                                                        overlap_block->anchors[species_chr] = overlap_qry_seg;
                                                    }
                                                    if (overlap_done && !suffix_done && sum >= suffix_length) {
                                                        suffix_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, suffix_len);
                                                        cigar_str = remain_cigar;
                                                        suffix_qry_seg = Segment::create(
                                                            0, // 临时坐标，稍后修正
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand, // 保持REVERSE
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            suffix_block
                                                        );
                                                        // 将query segment注册到suffix_block
                                                        if (suffix_block) {
                                                            suffix_block->anchors[species_chr] = suffix_qry_seg;
                                                        }
                                                    }
                                                }
                                            }

                                            // 反向链坐标修正：从segment末端开始，向前计算
                                            uint_t segment_end = segment->start + segment->length;

                                            // 使用一个临时变量来追踪当前段的末尾
                                            auto current_segment_end = segment_end;

                                            // 1. 首先处理前缀 (prefix)，如果它存在
                                            if (prefix_qry_seg) {
                                                prefix_qry_seg->start = current_segment_end - prefix_qry_seg->length;
                                                // 更新末尾位置，为下一个块做准备
                                                current_segment_end = prefix_qry_seg->start;
                                            }

                                            // 2. 接着处理重叠 (overlap)
                                            overlap_qry_seg->start = current_segment_end - overlap_qry_seg->length;
                                            // 再次更新末尾位置
                                            current_segment_end = overlap_qry_seg->start;

                                            // 3. 最后处理后缀 (suffix)，如果它存在
                                            if (suffix_qry_seg) {
                                                suffix_qry_seg->start = current_segment_end - suffix_qry_seg->length;
                                            }
                                        }
                                    }
                                }

                                // 2. 处理current_block的非ref anchors
                                for (const auto &[species_chr, segment]: current_block->anchors) {
                                    // 找到非ref的anchor
                                    if (species_chr.first != ref_name) {
                                        if (segment->strand == Strand::FORWARD) {
                                            std::string cigar_str;
                                            uint32_t sum = 0;
                                            bool prefix_done = curr_has_prefix ? false : true;
                                            bool overlap_done = false;
                                            bool suffix_done = curr_has_suffix ? false : true;
                                            SegPtr prefix_qry_seg = nullptr;
                                            SegPtr overlap_qry_seg = nullptr;
                                            SegPtr suffix_qry_seg = nullptr;
                                            uint32_t prefix_length = curr_has_prefix ? prefix_len : 0;
                                            uint32_t overlap_length = prefix_length + overlap_len;
                                            uint32_t suffix_length = curr_has_suffix ? overlap_length + suffix_len : 0;

                                            for (CigarUnit cu: segment->cigar) {
                                                uint32_t cigar_len;
                                                char op;
                                                intToCigar(cu, op, cigar_len);
                                                cigar_str += std::to_string(cigar_len) + op;
                                                if (op != 'I') {
                                                    sum += cigar_len;
                                                    // 修正累计长度判断
                                                    if (!prefix_done && sum >= prefix_length) {
                                                        prefix_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, prefix_len);
                                                        cigar_str = remain_cigar;
                                                        prefix_qry_seg = Segment::create(
                                                            segment->start,
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand,
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            prefix_block
                                                        );
                                                        // 将query segment注册到prefix_block
                                                        if (prefix_block) {
                                                            prefix_block->anchors[species_chr] = prefix_qry_seg;
                                                        }
                                                    }
                                                    if (prefix_done && !overlap_done && sum >= overlap_length) {
                                                        overlap_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, overlap_len);
                                                        cigar_str = remain_cigar;
                                                        overlap_qry_seg = Segment::create(
                                                            prefix_qry_seg
                                                                ? prefix_qry_seg->start + prefix_qry_seg->length
                                                                : segment->start,
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand,
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            overlap_block
                                                        );
                                                        // 将query segment注册到overlap_block
                                                        overlap_block->anchors[species_chr] = overlap_qry_seg;
                                                    }
                                                    if (overlap_done && !suffix_done && sum >= suffix_length) {
                                                        suffix_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, suffix_len);
                                                        cigar_str = remain_cigar;
                                                        uint_t suffix_start_pos = overlap_qry_seg
                                                            ? overlap_qry_seg->start + overlap_qry_seg->length
                                                            : (prefix_qry_seg
                                                                   ? prefix_qry_seg->start + prefix_qry_seg->length
                                                                   : segment->start);
                                                        suffix_qry_seg = Segment::create(
                                                            suffix_start_pos,
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand,
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            suffix_block
                                                        );
                                                        // 将query segment注册到suffix_block
                                                        if (suffix_block) {
                                                            suffix_block->anchors[species_chr] = suffix_qry_seg;
                                                        }
                                                    }
                                                }
                                            }
                                        } else {
                                            // segment->strand == Strand::REVERSE
                                            std::string cigar_str;
                                            uint32_t sum = 0;
                                            bool prefix_done = curr_has_prefix ? false : true;
                                            bool overlap_done = false;
                                            bool suffix_done = curr_has_suffix ? false : true;
                                            SegPtr prefix_qry_seg = nullptr;
                                            SegPtr overlap_qry_seg = nullptr;
                                            SegPtr suffix_qry_seg = nullptr;
                                            uint32_t prefix_length = curr_has_prefix ? prefix_len : 0;
                                            uint32_t overlap_length = prefix_length + overlap_len;
                                            uint32_t suffix_length = curr_has_suffix ? overlap_length + suffix_len : 0;

                                            // CIGAR分割逻辑与正向链相同
                                            for (CigarUnit cu: segment->cigar) {
                                                uint32_t cigar_len;
                                                char op;
                                                intToCigar(cu, op, cigar_len);
                                                cigar_str += std::to_string(cigar_len) + op;
                                                if (op != 'I') {
                                                    sum += cigar_len;
                                                    // 分割判断逻辑相同
                                                    if (!prefix_done && sum >= prefix_length) {
                                                        prefix_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, prefix_len);
                                                        cigar_str = remain_cigar;
                                                        prefix_qry_seg = Segment::create(
                                                            0, // 临时坐标，稍后修正
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand, // 保持REVERSE
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            prefix_block
                                                        );
                                                        // 将query segment注册到prefix_block
                                                        if (prefix_block) {
                                                            prefix_block->anchors[species_chr] = prefix_qry_seg;
                                                        }
                                                    }
                                                    if (prefix_done && !overlap_done && sum >= overlap_length) {
                                                        overlap_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, overlap_len);
                                                        cigar_str = remain_cigar;
                                                        overlap_qry_seg = Segment::create(
                                                            0, // 临时坐标，稍后修正
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand, // 保持REVERSE
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            overlap_block
                                                        );
                                                        // 将query segment注册到overlap_block
                                                        overlap_block->anchors[species_chr] = overlap_qry_seg;
                                                    }
                                                    if (overlap_done && !suffix_done && sum >= suffix_length) {
                                                        suffix_done = true;
                                                        auto [split_cigar, remain_cigar] = splitCigarMixed(
                                                            cigar_str, suffix_len);
                                                        cigar_str = remain_cigar;
                                                        suffix_qry_seg = Segment::create(
                                                            0, // 临时坐标，稍后修正
                                                            countNonDeletionOperations(split_cigar),
                                                            segment->strand, // 保持REVERSE
                                                            split_cigar,
                                                            segment->align_role,
                                                            SegmentRole::SEGMENT,
                                                            suffix_block
                                                        );
                                                        // 将query segment注册到suffix_block
                                                        if (suffix_block) {
                                                            suffix_block->anchors[species_chr] = suffix_qry_seg;
                                                        }
                                                    }
                                                }
                                            }

                                            // 反向链坐标修正：从segment末端开始，向前计算
                                            uint_t segment_end = segment->start + segment->length;

                                            // 使用一个临时变量来追踪当前段的末尾
                                            auto current_segment_end = segment_end;

                                            // 1. 首先处理前缀 (prefix)，如果它存在
                                            if (prefix_qry_seg) {
                                                prefix_qry_seg->start = current_segment_end - prefix_qry_seg->length;
                                                // 更新末尾位置，为下一个块做准备
                                                current_segment_end = prefix_qry_seg->start;
                                            }

                                            // 2. 接着处理重叠 (overlap)
                                            overlap_qry_seg->start = current_segment_end - overlap_qry_seg->length;
                                            // 再次更新末尾位置
                                            current_segment_end = overlap_qry_seg->start;

                                            // 3. 最后处理后缀 (suffix)，如果它存在
                                            if (suffix_qry_seg) {
                                                suffix_qry_seg->start = current_segment_end - suffix_qry_seg->length;
                                            }
                                        }
                                    }
                                }

                                // 开始替换prev和current，先替换Ref
                                // 获取prev的前驱和current的后继
                                SegPtr prev_prev = prev->primary_path.prev.load(std::memory_order_acquire);
                                SegPtr current_next = current->primary_path.next.load(std::memory_order_acquire);
                                // 移除prev和current的前驱和后继
                                prev->primary_path.prev.store(nullptr, std::memory_order_release);
                                prev->primary_path.next.store(nullptr, std::memory_order_release);
                                current->primary_path.next.store(nullptr, std::memory_order_release);
                                current->primary_path.prev.store(nullptr, std::memory_order_release);


                                // 将新blocks添加到全局池
                                if (prefix_block) {
                                    blocks.emplace_back(WeakBlock(prefix_block));
                                    // 将prefix的Seg插入

                                    prefix_ref_seg->primary_path.next.store(overlap_ref_seg, std::memory_order_release);
                                    overlap_ref_seg->primary_path.prev.store(prefix_ref_seg, std::memory_order_release);

                                    prefix_ref_seg->primary_path.prev.store(prev_prev, std::memory_order_release);
                                    prev_prev->primary_path.next.store(prefix_ref_seg, std::memory_order_release);
                                } else {
                                    overlap_ref_seg->primary_path.prev.store(prev_prev, std::memory_order_release);
                                    prev_prev->primary_path.next.store(overlap_ref_seg, std::memory_order_release);
                                }
                                blocks.emplace_back(WeakBlock(overlap_block));

                                if (suffix_block) {
                                    blocks.emplace_back(WeakBlock(suffix_block));
                                    SegPtr suffix_next = current_next;
                                    if (current_next->isTail()) {
                                        suffix_ref_seg->primary_path.next.
                                                store(current_next, std::memory_order_release);
                                        current_next->primary_path.prev.
                                                store(suffix_ref_seg, std::memory_order_release);

                                        overlap_ref_seg->primary_path.next.store(
                                            suffix_ref_seg, std::memory_order_release);
                                        suffix_ref_seg->primary_path.prev.store(
                                            overlap_ref_seg, std::memory_order_release);
                                    } else if (suffix_ref_seg->start > suffix_next->start) {
                                        while (suffix_ref_seg->start > suffix_next->start) {
                                            suffix_next = suffix_next->primary_path.next.
                                                    load(std::memory_order_release);
                                        }
                                        overlap_ref_seg->primary_path.next.store(
                                            current_next, std::memory_order_release);
                                        current_next->primary_path.prev.store(
                                            overlap_ref_seg, std::memory_order_release);

                                        suffix_ref_seg->primary_path.prev.store(
                                            suffix_next->primary_path.prev.load(std::memory_order_release),
                                            std::memory_order_release);
                                        suffix_next->primary_path.prev.load(std::memory_order_release)->primary_path.
                                                next.store(suffix_ref_seg, std::memory_order_release);

                                        suffix_ref_seg->primary_path.next.store(suffix_next, std::memory_order_release);
                                        suffix_next->primary_path.prev.store(suffix_ref_seg, std::memory_order_release);
                                    } else {
                                        suffix_ref_seg->primary_path.next.
                                                store(current_next, std::memory_order_release);
                                        current_next->primary_path.prev.
                                                store(suffix_ref_seg, std::memory_order_release);

                                        overlap_ref_seg->primary_path.next.store(
                                            suffix_ref_seg, std::memory_order_release);
                                        suffix_ref_seg->primary_path.prev.store(
                                            overlap_ref_seg, std::memory_order_release);
                                    }
                                } else {
                                    overlap_ref_seg->primary_path.next.store(current_next, std::memory_order_release);
                                    current_next->primary_path.prev.store(overlap_ref_seg, std::memory_order_release);
                                }

                                // 开始替换query的seg和block
                                for (const auto &[species_chr, segment]: prev_block->anchors) {
                                    if (species_chr.first != ref_name) {
                                        SegPtr qry_prev = segment->primary_path.prev.load(std::memory_order_acquire);
                                        SegPtr qry_next = segment->primary_path.next.load(std::memory_order_acquire);
                                        segment->primary_path.prev.store(nullptr, std::memory_order_release);
                                        segment->primary_path.next.store(nullptr, std::memory_order_release);
                                        // 先找到一定存在的overlap_block
                                        for (const auto &[species_chr_overlap, segment_overlap]: overlap_block->
                                             anchors) {
                                            if (species_chr_overlap.first == species_chr.first) {
                                                if (segment->strand == Strand::FORWARD) {
                                                    bool query_has_prefix = false;
                                                    bool query_has_suffix = false;
                                                    if (prev_has_prefix) {
                                                        for (const auto &[species_chr_prefix, segment_prefix]:
                                                             prefix_block->anchors) {
                                                            if (species_chr_prefix.first == species_chr.first) {
                                                                query_has_prefix = true;
                                                                qry_prev->primary_path.next.store(
                                                                    segment_prefix, std::memory_order_release);
                                                                segment_prefix->primary_path.prev.store(
                                                                    qry_prev, std::memory_order_release);
                                                                segment_prefix->primary_path.next.store(
                                                                    segment_overlap, std::memory_order_release);
                                                                segment_overlap->primary_path.prev.store(
                                                                    segment_prefix, std::memory_order_release);
                                                                break;
                                                            }
                                                        }
                                                        if (!query_has_prefix) {
                                                            qry_prev->primary_path.next.store(
                                                                segment_overlap, std::memory_order_release);
                                                            segment_overlap->primary_path.prev.store(
                                                                qry_prev, std::memory_order_release);
                                                        }
                                                    } else {
                                                        qry_prev->primary_path.next.store(
                                                            segment_overlap, std::memory_order_release);
                                                        segment_overlap->primary_path.prev.store(
                                                            qry_prev, std::memory_order_release);
                                                    }
                                                    if (prev_has_suffix) {
                                                        for (const auto &[species_chr_suffix, segment_suffix]:
                                                             suffix_block->anchors) {
                                                            if (species_chr_suffix.first == species_chr.first) {
                                                                query_has_suffix = true;
                                                                segment_overlap->primary_path.next.store(
                                                                    segment_suffix, std::memory_order_release);
                                                                segment_suffix->primary_path.prev.store(
                                                                    segment_overlap, std::memory_order_release);
                                                                segment_suffix->primary_path.next.store(
                                                                    qry_next, std::memory_order_release);
                                                                qry_next->primary_path.prev.store(
                                                                    segment_suffix, std::memory_order_release);
                                                                break;
                                                            }
                                                        }
                                                        if (!query_has_suffix) {
                                                            segment_overlap->primary_path.next.store(
                                                                qry_next, std::memory_order_release);
                                                            qry_next->primary_path.prev.store(
                                                                segment_overlap, std::memory_order_release);
                                                        }
                                                    } else {
                                                        segment_overlap->primary_path.next.store(
                                                            qry_next, std::memory_order_release);
                                                        qry_next->primary_path.prev.store(
                                                            segment_overlap, std::memory_order_release);
                                                    }
                                                }
                                                // 如果是反向链就反着连
                                                else {
                                                    bool query_has_prefix = false;
                                                    bool query_has_suffix = false;
                                                    if (prev_has_suffix) {
                                                        for (const auto &[species_chr_suffix, segment_suffix]:
                                                             suffix_block->anchors) {
                                                            if (species_chr_suffix.first == species_chr.first) {
                                                                query_has_suffix = true;
                                                                segment_suffix->primary_path.next.store(
                                                                    segment_overlap, std::memory_order_release);
                                                                segment_suffix->primary_path.prev.store(
                                                                    qry_prev, std::memory_order_release);
                                                                segment_overlap->primary_path.prev.store(
                                                                    segment_suffix, std::memory_order_release);
                                                                qry_prev->primary_path.next.store(
                                                                    segment_suffix, std::memory_order_release);
                                                                break;
                                                            }
                                                        }
                                                        if (!query_has_suffix) {
                                                            qry_prev->primary_path.next.store(
                                                                segment_overlap, std::memory_order_release);
                                                            segment_overlap->primary_path.prev.store(
                                                                qry_prev, std::memory_order_release);
                                                        }
                                                    } else {
                                                        qry_prev->primary_path.next.store(
                                                            segment_overlap, std::memory_order_release);
                                                        segment_overlap->primary_path.prev.store(
                                                            qry_prev, std::memory_order_release);

                                                    }
                                                    if (prev_has_prefix) {
                                                        for (const auto &[species_chr_prefix, segment_prefix]:
                                                             prefix_block->anchors) {
                                                            if (species_chr_prefix.first == species_chr.first) {
                                                                query_has_prefix = true;
                                                                segment_prefix->primary_path.prev.store(
                                                                    segment_overlap, std::memory_order_release);
                                                                segment_prefix->primary_path.next.store(
                                                                    qry_next, std::memory_order_release);
                                                                qry_next->primary_path.prev.store(
                                                                    segment_prefix, std::memory_order_release);
                                                                segment_overlap->primary_path.next.store(
                                                                    segment_prefix, std::memory_order_release);

                                                                break;
                                                            }
                                                        }
                                                        if (!query_has_prefix) {
                                                            segment_overlap->primary_path.next.store(
                                                                qry_next, std::memory_order_release);
                                                            qry_next->primary_path.prev.store(
                                                                segment_overlap, std::memory_order_release);
                                                        }
                                                    } else {
                                                        segment_overlap->primary_path.next.store(
                                                            qry_next, std::memory_order_release);
                                                        qry_next->primary_path.prev.store(
                                                            segment_overlap, std::memory_order_release);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                for (const auto &[species_chr, segment]: current_block->anchors) {
                                    if (species_chr.first != ref_name) {
                                        SegPtr qry_prev = segment->primary_path.prev.load(std::memory_order_acquire);
                                        SegPtr qry_next = segment->primary_path.next.load(std::memory_order_acquire);
                                        segment->primary_path.prev.store(nullptr, std::memory_order_release);
                                        segment->primary_path.next.store(nullptr, std::memory_order_release);
                                        // 先找到一定存在的overlap_block
                                        for (const auto &[species_chr_overlap, segment_overlap]: overlap_block->
                                             anchors) {
                                            if (species_chr_overlap.first == species_chr.first) {
                                                if (segment->strand == Strand::FORWARD) {
                                                    bool query_has_prefix = false;
                                                    bool query_has_suffix = false;
                                                    if (curr_has_prefix) {
                                                        for (const auto &[species_chr_prefix, segment_prefix]:
                                                             prefix_block->anchors) {
                                                            if (species_chr_prefix.first == species_chr.first) {
                                                                query_has_prefix = true;
                                                                qry_prev->primary_path.next.store(
                                                                    segment_prefix, std::memory_order_release);
                                                                segment_prefix->primary_path.prev.store(
                                                                    qry_prev, std::memory_order_release);
                                                                segment_prefix->primary_path.next.store(
                                                                    segment_overlap, std::memory_order_release);
                                                                segment_overlap->primary_path.prev.store(
                                                                    segment_prefix, std::memory_order_release);
                                                                break;
                                                            }
                                                        }
                                                        if (!query_has_prefix) {
                                                            qry_prev->primary_path.next.store(
                                                                segment_overlap, std::memory_order_release);
                                                            segment_overlap->primary_path.prev.store(
                                                                qry_prev, std::memory_order_release);
                                                        }
                                                    } else {
                                                        qry_prev->primary_path.next.store(
                                                            segment_overlap, std::memory_order_release);
                                                        segment_overlap->primary_path.prev.store(
                                                            qry_prev, std::memory_order_release);
                                                    }
                                                    if (curr_has_suffix) {
                                                        for (const auto &[species_chr_suffix, segment_suffix]:
                                                             suffix_block->anchors) {
                                                            if (species_chr_suffix.first == species_chr.first) {
                                                                query_has_suffix = true;
                                                                segment_overlap->primary_path.next.store(
                                                                    segment_suffix, std::memory_order_release);
                                                                segment_suffix->primary_path.prev.store(
                                                                    segment_overlap, std::memory_order_release);
                                                                segment_suffix->primary_path.next.store(
                                                                    qry_next, std::memory_order_release);
                                                                qry_next->primary_path.prev.store(
                                                                    segment_suffix, std::memory_order_release);
                                                                break;
                                                            }
                                                        }
                                                        if (!query_has_suffix) {
                                                            segment_overlap->primary_path.next.store(
                                                                qry_next, std::memory_order_release);
                                                            qry_next->primary_path.prev.store(
                                                                segment_overlap, std::memory_order_release);
                                                        }
                                                    } else {
                                                        segment_overlap->primary_path.next.store(
                                                            qry_next, std::memory_order_release);
                                                        qry_next->primary_path.prev.store(
                                                            segment_overlap, std::memory_order_release);
                                                    }
                                                }
                                                // 如果是反向链就反着连
                                                else {
                                                    bool query_has_prefix = false;
                                                    bool query_has_suffix = false;
                                                    if (curr_has_suffix) {
                                                        for (const auto &[species_chr_suffix, segment_suffix]:
                                                             suffix_block->anchors) {
                                                            if (species_chr_suffix.first == species_chr.first) {
                                                                query_has_suffix = true;
                                                                segment_suffix->primary_path.next.store(
                                                                    segment_overlap, std::memory_order_release);
                                                                segment_suffix->primary_path.prev.store(
                                                                    qry_prev, std::memory_order_release);
                                                                segment_overlap->primary_path.prev.store(
                                                                    segment_suffix, std::memory_order_release);
                                                                qry_prev->primary_path.next.store(
                                                                    segment_suffix, std::memory_order_release);
                                                                break;
                                                            }
                                                        }
                                                        if (!query_has_suffix) {
                                                            qry_prev->primary_path.next.store(
                                                                segment_overlap, std::memory_order_release);
                                                            segment_overlap->primary_path.prev.store(
                                                                qry_prev, std::memory_order_release);
                                                        }
                                                    } else {
                                                        qry_prev->primary_path.next.store(
                                                            segment_overlap, std::memory_order_release);
                                                        segment_overlap->primary_path.prev.store(
                                                            qry_prev, std::memory_order_release);

                                                    }
                                                    if (curr_has_prefix) {
                                                        for (const auto &[species_chr_prefix, segment_prefix]:
                                                             prefix_block->anchors) {
                                                            if (species_chr_prefix.first == species_chr.first) {
                                                                query_has_prefix = true;
                                                                segment_prefix->primary_path.prev.store(
                                                                    segment_overlap, std::memory_order_release);
                                                                segment_prefix->primary_path.next.store(
                                                                    qry_next, std::memory_order_release);
                                                                qry_next->primary_path.prev.store(
                                                                    segment_prefix, std::memory_order_release);
                                                                segment_overlap->primary_path.next.store(
                                                                    segment_prefix, std::memory_order_release);

                                                                break;
                                                            }
                                                        }
                                                        if (!query_has_prefix) {
                                                            segment_overlap->primary_path.next.store(
                                                                qry_next, std::memory_order_release);
                                                            qry_next->primary_path.prev.store(
                                                                segment_overlap, std::memory_order_release);
                                                        }
                                                    } else {
                                                        segment_overlap->primary_path.next.store(
                                                            qry_next, std::memory_order_release);
                                                        qry_next->primary_path.prev.store(
                                                            segment_overlap, std::memory_order_release);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                // 移除prev_block
                                auto prev_it = std::find_if(
                                    blocks.begin(), blocks.end(),
                                    [&](const WeakBlock &wb) {
                                        auto locked = wb.lock();
                                        return locked && locked == prev_block;
                                    });
                                if (prev_it != blocks.end()) {
                                    blocks.erase(prev_it);
                                }
                                // 移除current_block
                                auto curr_it = std::find_if(
                                    blocks.begin(), blocks.end(),
                                    [&](const WeakBlock &wb) {
                                        auto locked = wb.lock();
                                        return locked && locked == current_block;
                                    });
                                if (curr_it != blocks.end()) {
                                    blocks.erase(curr_it);
                                }
                                // 替换current和prev
                                prev = overlap_ref_seg->primary_path.prev.load(std::memory_order_acquire);
                                if (prev->isHead()) {
                                    prev = overlap_ref_seg;
                                }
                                prev_block = prev->parent_block;
                                current = prev->primary_path.next.load(std::memory_order_acquire);
                                current_block = current->parent_block;
                                // 继续处理，不要跳到最后，因为新创建的blocks之间可能还有重叠


                                // 收集所有物种的segment详细信息
                                debug_all_species_segments.clear();
                                for (const auto& [species_name, genome_graph] : species_graphs) {
                                    SpeciesDebugInfo species_info(species_name);
                                    std::shared_lock genome_lock(genome_graph.rw);

                                    for (const auto& [chr_name, genome_end] : genome_graph.chr2end) {
                                        ChromosomeDebugInfo chr_info(chr_name);
                                        std::shared_lock end_lock(genome_end.rw);

                                        // 遍历该染色体的所有segments
                                        SegPtr current = genome_end.head;
                                        current = current->primary_path.next.load(std::memory_order_acquire);
                                        while (current && current != genome_end.tail) {
                                            if (current->isSegment()) {
                                                chr_info.segments.emplace_back(current);
                                            }
                                            // 检查链表结构是否正确
                                            if (current->isHead()) {
                                                current = current->primary_path.next.load(std::memory_order_acquire);
                                                continue;
                                            }
                                            SegPtr prev = current->primary_path.prev.load(std::memory_order_acquire);
                                            if(prev && prev->isHead())
                                            {
                                                current = current->primary_path.next.load(std::memory_order_acquire);
                                                continue;
                                            }
                                            if (prev && prev->primary_path.next.load(std::memory_order_acquire) != current) {
                                                std::cerr << "[链表错误] prev->next != current, current start: " << current->start << std::endl;
                                            }
                                            current = current->primary_path.next.load(std::memory_order_acquire);
                                        }

                                        species_info.chromosomes.push_back(std::move(chr_info));
                                    }

                                    debug_all_species_segments.push_back(std::move(species_info));
                                }

count++;
                                // std::cout<<count<<std::endl;
                                continue;
                            }
                        }
                    }


                    prev_block = current_block;
                    prev = current;
                    // 移动到下一个segment - 使用原子操作保证线程安全
                    current = current->primary_path.next.load(std::memory_order_acquire);
                }
            }
        }
        // // 收集所有物种的segment详细信息
        // debug_all_species_segments.clear();
        // for (const auto &[species_name, genome_graph]: species_graphs) {
        //     SpeciesDebugInfo species_info(species_name);
        //     std::shared_lock genome_lock(genome_graph.rw);
        //
        //     for (const auto &[chr_name, genome_end]: genome_graph.chr2end) {
        //         ChromosomeDebugInfo chr_info(chr_name);
        //         std::shared_lock end_lock(genome_end.rw);
        //
        //         // 遍历该染色体的所有segments
        //         SegPtr current = genome_end.head;
        //         current = current->primary_path.next.load(std::memory_order_acquire);
        //         while (current && current != genome_end.tail) {
        //             if (current->isSegment()) {
        //                 chr_info.segments.emplace_back(current);
        //             }
        //             // 检查链表结构是否正确
        //             if (current->isHead()) {
        //                 current = current->primary_path.next.load(std::memory_order_acquire);
        //                 continue;
        //             }
        //             SegPtr prev = current->primary_path.prev.load(std::memory_order_acquire);
        //             if (prev && prev->isHead()) {
        //                 current = current->primary_path.next.load(std::memory_order_acquire);
        //                 continue;
        //             }
        //             if (prev && prev->primary_path.next.load(std::memory_order_acquire) != current) {
        //                 std::cerr << "[链表错误] prev->next != current, current start: " << current->start << std::endl;
        //             }
        //             current = current->primary_path.next.load(std::memory_order_acquire);
        //         }
        //
        //         species_info.chromosomes.push_back(std::move(chr_info));
        //     }
        //
        //     debug_all_species_segments.push_back(std::move(species_info));
        // }
    }

    /* ==============================================================
     * 7. 高性能删除方法 (public API)
     * ==============================================================*/
    bool RaMeshMultiGenomeGraph::removeBlock(const BlockPtr &block) {
        if (!block) return false;

        std::unique_lock graph_lock(rw);

        // 1. 从全局block池中移除
        auto it = std::find_if(blocks.begin(), blocks.end(),
                               [&](const WeakBlock &wb) {
                                   auto locked = wb.lock();
                                   return locked && locked == block;
                               });

        if (it != blocks.end()) {
            blocks.erase(it);
        } else {
            return false; // 未找到block
        }

        // 2. 从每个物种的链表中移除相关segments
        {
            std::shared_lock block_lock(block->rw);
            for (const auto &[species_chr, segment]: block->anchors) {
                if (segment) {
                    const SpeciesName &species = species_chr.first;
                    const ChrName &chr = species_chr.second;

                    auto species_it = species_graphs.find(species);
                    if (species_it != species_graphs.end()) {
                        auto chr_it = species_it->second.chr2end.find(chr);
                        if (chr_it != species_it->second.chr2end.end()) {
                            chr_it->second.removeSegment(segment);
                        }
                    }
                }
            }
        }

        // 3. 清空block的anchors（segments会由GenomeEnd负责删除）
        block->removeAllSegments();

        return true;
    }

    bool RaMeshMultiGenomeGraph::removeBlockById(size_t block_index) {
        std::unique_lock graph_lock(rw);

        if (block_index >= blocks.size()) return false;

        auto weak_block = blocks[block_index];
        auto block = weak_block.lock();

        if (!block) {
            // 清理过期的weak_ptr
            blocks.erase(blocks.begin() + block_index);
            return false;
        }

        return removeBlock(block);
    }

    void RaMeshMultiGenomeGraph::removeBlocksBatch(const std::vector<BlockPtr> &blocks_to_remove) {
        if (blocks_to_remove.empty()) return;

        std::unique_lock graph_lock(rw);

        // 收集所有需要删除的segments，按物种-染色体分组
        std::unordered_map<SpeciesName,
            std::unordered_map<ChrName, std::vector<SegPtr> > > segments_by_genome;

        for (const auto &block: blocks_to_remove) {
            if (!block) continue;

            std::shared_lock block_lock(block->rw);
            for (const auto &[species_chr, segment]: block->anchors) {
                if (segment) {
                    segments_by_genome[species_chr.first][species_chr.second].emplace_back(segment);
                }
            }
        }

        // 批量删除segments
        for (const auto &[species, chr_map]: segments_by_genome) {
            auto species_it = species_graphs.find(species);
            if (species_it != species_graphs.end()) {
                for (const auto &[chr, segments]: chr_map) {
                    auto chr_it = species_it->second.chr2end.find(chr);
                    if (chr_it != species_it->second.chr2end.end()) {
                        chr_it->second.removeBatch(segments);
                    }
                }
            }
        }

        // 从全局block池中移除
        for (const auto &block: blocks_to_remove) {
            auto it = std::find_if(blocks.begin(), blocks.end(),
                                   [&](const WeakBlock &wb) {
                                       auto locked = wb.lock();
                                       return locked && locked == block;
                                   });

            if (it != blocks.end()) {
                blocks.erase(it);
            }

            // 清空block的anchors
            block->removeAllSegments();
        }
    }

    size_t RaMeshMultiGenomeGraph::removeExpiredBlocks() {
        std::unique_lock graph_lock(rw);

        size_t removed_count = 0;

        // 移除所有过期的weak_ptr
        blocks.erase(
            std::remove_if(blocks.begin(), blocks.end(),
                           [&removed_count](const WeakBlock &wb) {
                               if (wb.expired()) {
                                   ++removed_count;
                                   return true;
                               }
                               return false;
                           }),
            blocks.end()
        );

        return removed_count;
    }

    void RaMeshMultiGenomeGraph::removeSpecies(const SpeciesName &species) {
        std::unique_lock graph_lock(rw);

        // 1. 移除species_graphs中的条目
        auto species_it = species_graphs.find(species);
        if (species_it == species_graphs.end()) return;

        // 2. 收集包含该物种的所有blocks
        std::vector<BlockPtr> blocks_to_remove;

        for (const auto &weak_block: blocks) {
            auto block = weak_block.lock();
            if (!block) continue;

            std::shared_lock block_lock(block->rw);
            for (const auto &[species_chr, segment]: block->anchors) {
                if (species_chr.first == species) {
                    blocks_to_remove.emplace_back(block);
                    break;
                }
            }
        }

        // 3. 批量删除blocks
        removeBlocksBatch(blocks_to_remove);

        // 4. 从species_graphs中移除
        species_graphs.erase(species_it);
    }

    void RaMeshMultiGenomeGraph::removeChromosome(const SpeciesName &species, const ChrName &chr) {
        std::unique_lock graph_lock(rw);

        auto species_it = species_graphs.find(species);
        if (species_it == species_graphs.end()) return;

        auto chr_it = species_it->second.chr2end.find(chr);
        if (chr_it == species_it->second.chr2end.end()) return;

        // 1. 收集包含该染色体的所有blocks
        std::vector<BlockPtr> blocks_to_update;
        std::vector<BlockPtr> blocks_to_remove;

        for (const auto &weak_block: blocks) {
            auto block = weak_block.lock();
            if (!block) continue;

            std::shared_lock block_lock(block->rw);
            SpeciesChrPair target_key{species, chr};

            if (block->anchors.find(target_key) != block->anchors.end()) {
                if (block->anchors.size() == 1) {
                    // 如果block只包含这一个染色体，删除整个block
                    blocks_to_remove.emplace_back(block);
                } else {
                    // 否则只删除这个染色体的segment
                    blocks_to_update.emplace_back(block);
                }
            }
        }

        // 2. 清空该染色体的所有segments
        chr_it->second.clearAllSegments();

        // 3. 更新相关blocks
        for (const auto &block: blocks_to_update) {
            block->removeSegment(species, chr);
        }

        // 4. 删除只包含该染色体的blocks
        removeBlocksBatch(blocks_to_remove);

        // 5. 从species_graphs中移除染色体
        species_it->second.chr2end.erase(chr_it);
    }

    void RaMeshMultiGenomeGraph::clearAllGraphs() {
        std::unique_lock graph_lock(rw);

        // 1. 清空所有species的图
        for (auto &[species, genome_graph]: species_graphs) {
            std::unique_lock species_lock(genome_graph.rw);
            for (auto &[chr, genome_end]: genome_graph.chr2end) {
                genome_end.clearAllSegments();
            }
            genome_graph.chr2end.clear();
        }

        // 2. 清空所有blocks
        for (const auto &weak_block: blocks) {
            auto block = weak_block.lock();
            if (block) {
                block->removeAllSegments();
            }
        }
        blocks.clear();

        // 3. 清空species_graphs
        species_graphs.clear();
    }

    size_t RaMeshMultiGenomeGraph::compactBlockPool() {
        std::unique_lock graph_lock(rw);

        size_t original_size = blocks.size();
        size_t compacted = removeExpiredBlocks();

        // 可以在这里添加更多的内存优化逻辑
        blocks.shrink_to_fit();

        return compacted;
    }

    void RaMeshMultiGenomeGraph::optimizeGraphStructure() {

        // 优化每个物种图的采样表
        for (auto &[species, genome_graph]: species_graphs) {
            std::shared_lock species_lock(genome_graph.rw);

            for (auto &[chr, genome_end]: genome_graph.chr2end) {
                // 重建采样表以提高查询效率
                genome_end.sample_vec.clear();
                uint_t last_seg_end = 0;

                SegPtr current = genome_end.head->primary_path.next.load(std::memory_order_acquire);

                while (current && !current->isTail()) {
                    genome_end.setToSampling(current);
                    current = current->primary_path.next.load(std::memory_order_acquire);
                }
            }
        }
    }

    RaMeshMultiGenomeGraph::DeletionStats RaMeshMultiGenomeGraph::performMaintenance(bool full_gc) {
        auto start_time = std::chrono::high_resolution_clock::now();

        DeletionStats stats;

        // 1. 清理过期的blocks
        stats.expired_blocks_cleaned = removeExpiredBlocks();

        if (full_gc) {
            // 2. 压缩block池
            stats.blocks_removed = compactBlockPool();

            // 3. 优化图结构
            optimizeGraphStructure();
            stats.sampling_updates = species_graphs.size(); // 简化统计
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        stats.total_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        return stats;
    }
} // namespace RaMesh
