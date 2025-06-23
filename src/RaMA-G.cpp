// RaMAx.cpp: 定义应用程序的入口点。
// 主程序：负责解析命令行参数，初始化配置，执行基因组比对流程

#include "SeqPro.h"
#include "anchor.h"
#include "config.hpp"
#include "data_process.h"
#include "index.h"
#include "anchor.h"
#include "rare_aligner.h"
#include <limits>
#include <algorithm>
#include <variant>

// ------------------------------------------------------------------
// 通用命令行参数结构体（可序列化）
// ------------------------------------------------------------------
struct CommonArgs {
    std::filesystem::path reference_path = "";
    std::filesystem::path query_path = "";
    std::filesystem::path output_path = "";
    std::filesystem::path work_dir_path = "";

    uint_t chunk_size = 10000000;
    uint_t overlap_size = 100000;
    uint_t min_anchor_length = 20;
    uint_t max_anchor_frequency = 50;
    bool restart = false;
    int thread_num = std::thread::hardware_concurrency(); // 默认使用所有 CPU 核心
    PairGenomeOutputFormat output_format = PairGenomeOutputFormat::UNKNOWN;
    bool enable_repeat_masking = false; // 是否启用重复序列遮蔽，默认为 false

    // 支持 cereal 序列化
    template<class Archive>
    void serialize(Archive &ar) {
        ar(
            CEREAL_NVP(reference_path),
            CEREAL_NVP(query_path),
            CEREAL_NVP(output_path),
            CEREAL_NVP(work_dir_path),
            CEREAL_NVP(chunk_size),
            CEREAL_NVP(overlap_size),
            CEREAL_NVP(restart),
            CEREAL_NVP(thread_num),
            CEREAL_NVP(output_format),
            CEREAL_NVP(enable_repeat_masking)
        );
    }
};

// ------------------------------------------------------------------
// 辅助函数：计算并记录reference的序列统计信息
// ------------------------------------------------------------------
template<typename ManagerType>
void recordReferenceSequenceStats(const std::string& species_name, 
                                 const std::unique_ptr<ManagerType>& manager,
                                 SeqPro::Length& reference_min_seq_length) {
    auto seq_count = manager->getSequenceCount();
    
    if (species_name == "reference") {
        // 只对reference计算序列长度统计
        auto seq_names = manager->getSequenceNames();
        SeqPro::Length min_seq_length = std::numeric_limits<SeqPro::Length>::max();
        SeqPro::Length max_seq_length = 0;
        
        for (const auto& seq_name : seq_names) {
            auto seq_length = manager->getSequenceLength(seq_name);
            min_seq_length = std::min(min_seq_length, seq_length);
            max_seq_length = std::max(max_seq_length, seq_length);
        }
        
        spdlog::info("[{}] Loaded {} sequences, min length: {}, max length: {}", 
                     species_name, seq_count, min_seq_length, max_seq_length);
        reference_min_seq_length = min_seq_length;
    } else {
        spdlog::info("[{}] Loaded {} sequences", species_name, seq_count);
    }
}

// ------------------------------------------------------------------
// CLI11 参数注册：配置常用命令行参数（参考 RaMAx 主程序）
// ------------------------------------------------------------------
inline void setupCommonOptions(CLI::App *cmd, CommonArgs &args) {
    auto fmt = std::make_shared<CustomFormatter>();
    fmt->column_width(50);
    cmd->formatter(fmt);

    cmd->set_version_flag("-v,--version", std::string("RaMAx version ") + VERSION);

    auto *ref_opt = cmd->add_option("-r,--reference", args.reference_path,
                                    "Path to the reference genome file (FASTA format).")
            ->group("Input Files")
            ->type_name("<path>")->transform(trim_whitespace);

    auto *qry_opt = cmd->add_option("-q,--query", args.query_path,
                                    "Path to the query genome file (FASTA format).")
            ->group("Input Files")
            ->type_name("<path>")->transform(trim_whitespace);

    auto *output_opt = cmd->add_option("-o,--output", args.output_path,
                                       "Path to save alignment results.")
            ->group("Output")->type_name("<path>")
            ->transform(trim_whitespace);

    auto *workspace_opt = cmd->add_option("-w,--workdir", args.work_dir_path,
                                          "Path to the working directory for temporary files.")
            ->group("Output")->type_name("<path>")
            ->transform(trim_whitespace);

    auto *chunk_size_opt = cmd->add_option("--chunk_size", args.chunk_size,
                                           "Size of each chunk for parallel processing (default: 10000000).")
            ->default_val(10000000)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(1000000, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *overlap_size_opt = cmd->add_option("--overlap_size", args.overlap_size,
                                             "Size of overlap between chunks (default: 100000).")
            ->default_val(100000)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(0, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *min_anchor_length_opt = cmd->add_option(
                "--min_anchor_length", args.min_anchor_length,
                "Minimum anchor length (default: 20).")
            ->default_val(20)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(1, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *max_anchor_frequency_opt = cmd->add_option(
                "--max_anchor_frequency", args.max_anchor_frequency,
                "Maximum anchor frequency filter (default: 50).")
            ->default_val(50)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(0, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);


    auto *threads_opt = cmd->add_option("-t,--threads", args.thread_num,
                                        "Number of threads to use for parallel processing (default: system cores).")
            ->default_val(std::thread::hardware_concurrency())
            ->envname("RAMA_G_THREADS")
            ->capture_default_str()
            ->group("Performance")
            ->check(CLI::Range(1, std::numeric_limits<int>::max()))
            ->type_name("<int>")->transform(trim_whitespace);

    auto *restart_flag = cmd->add_flag("--restart", args.restart,
                                       "Restart the alignment process by skipping the existing index files.")
            ->group("Performance");


    /* auto* mask_repeats_flag = */
    cmd->add_flag("--mask-repeats", args.enable_repeat_masking,
                  "Enable repeat sequence masking.")
            ->group("Software Parameters");

    // Set dependencies and exclusions
    restart_flag->needs(workspace_opt);
    restart_flag->excludes(ref_opt,
                           qry_opt, output_opt, threads_opt,
                           chunk_size_opt, overlap_size_opt,
                           min_anchor_length_opt, max_anchor_frequency_opt);
}


int main(int argc, char **argv) {
    // 初始化异步日志线程池（spdlog）
    spdlog::init_thread_pool(8192, 1); // 日志缓冲区容量 8192，单线程日志写入
    // setupLogger();  // 控制台输出日志（不带文件）

    // 初始化 CLI 命令行应用
    CLI::App app{"RaMA-G: A High-performance Genome Alignment Tool"};

    CommonArgs common_args; // 存储用户输入的参数
    setupCommonOptions(&app, common_args); // 注册参数解析
    CLI11_PARSE(app, argc, argv); // 开始解析命令行参数

    try {
        // ------------------------------
        // 模式 1：重启模式（--restart）
        // ------------------------------
        if (common_args.restart) {
            if (common_args.work_dir_path.empty()) {
                throw CLI::RequiredError("In restart mode, --workdir (-w) is required.");
            }

            // 检查工作目录是否存在
            if (std::filesystem::exists(common_args.work_dir_path)) {
                if (!std::filesystem::is_directory(common_args.work_dir_path)) {
                    throw CLI::ValidationError("Work directory is not valid: " + common_args.work_dir_path.string());
                }
            } else {
                std::filesystem::create_directories(common_args.work_dir_path);
            }

            // 初始化日志器
            setupLoggerWithFile(common_args.work_dir_path);
            spdlog::info("RaMA-G version {}", VERSION);
            spdlog::info("Restart mode enabled.");

            // 加载之前保存的参数配置文件
            FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
            std::ifstream is(config_path);
            if (!is) {
                spdlog::error("Failed to open {} for loading CommonArgs", config_path.string());
                return false;
            }
            cereal::JSONInputArchive archive(is);
            archive(common_args); // 反序列化参数
            spdlog::info("CommonArgs loaded from {}", config_path.string());
        }

        // ------------------------------
        // 模式 2：正常运行模式
        // ------------------------------
        else {
            // 检查必要参数
            if (common_args.reference_path.empty())
                throw CLI::RequiredError("Missing required option: --reference (-r)");
            if (common_args.query_path.empty())
                throw CLI::RequiredError("Missing required option: --query (-q)");
            if (common_args.output_path.empty())
                throw CLI::RequiredError("Missing required option: --output (-o)");
            if (common_args.work_dir_path.empty())
                throw CLI::RequiredError("Missing required option: --workdir (-w)");

#ifndef _DEBUG_
			// 非调试模式下：确保工作目录为空
			if (std::filesystem::exists(common_args.work_dir_path)) {
				if (!std::filesystem::is_directory(common_args.work_dir_path)) {
					throw CLI::ValidationError("Work directory is not valid: " + common_args.work_dir_path.string());
				}
				if (!std::filesystem::is_empty(common_args.work_dir_path)) {
					throw CLI::ValidationError("Work directory is not empty: " + common_args.work_dir_path.string());
				}
			}
			else {
				std::filesystem::create_directories(common_args.work_dir_path);
			}
#endif

            setupLoggerWithFile(common_args.work_dir_path);
            spdlog::info("RaMA-G version {}", VERSION);
            spdlog::info("Alignment mode enabled.");

            // 验证参考文件路径
            std::string ref_str = common_args.reference_path.string();
            if (isUrl(ref_str)) {
                verifyUrlReachable(ref_str);
            } else {
                verifyLocalFile(common_args.reference_path);
            }

            // 验证查询文件路径
            std::string qry_str = common_args.query_path.string();
            if (isUrl(qry_str)) {
                verifyUrlReachable(qry_str);
            } else {
                verifyLocalFile(common_args.query_path);
            }

            // 自动识别输出格式（支持 .sam / .maf / .paf）
            common_args.output_format = detectPairGenomeOutputFormat(common_args.output_path);
            if (common_args.output_format == PairGenomeOutputFormat::UNKNOWN) {
                throw std::runtime_error("Invalid output file extension. Supported: .sam, .maf, .paf");
            }

            // 检查 chunk 和 overlap 参数
            if (common_args.overlap_size >= common_args.chunk_size) {
                throw std::runtime_error("Overlap size must be less than chunk size.");
            }

            // 保存参数配置文件（用于 --restart）
            FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
            std::ofstream os(config_path);
            if (!os) {
                spdlog::error("Failed to open {} for saving CommonArgs", config_path.string());
                return false;
            }
            cereal::JSONOutputArchive archive(os);
            archive(cereal::make_nvp("common_args", common_args));
            spdlog::info("CommonArgs saved to {}", config_path.string());
        }
    } catch (const std::runtime_error &e) {
        spdlog::error("{}", e.what());
        spdlog::error("Use --help for usage information.");
        spdlog::error("Exiting with error code 1.");
        return 1;
    }

    // ------------------------------
    // 主流程开始
    // ------------------------------
    spdlog::info("Command: {}", getCommandLine(argc, argv));

    spdlog::info("Reference: {} (size: {})",
                 common_args.reference_path.string(),
                 getReadableFileSize(common_args.reference_path));

    spdlog::info("Query: {} (size: {})",
                 common_args.query_path.string(),
                 getReadableFileSize(common_args.query_path));

    spdlog::info("Output: {}", common_args.output_path.string());
    spdlog::info("Work directory: {}", common_args.work_dir_path.string());
    spdlog::info("Threads: {}", common_args.thread_num);

    // ------------------------------
    // 数据预处理阶段
    // ------------------------------
    SpeciesPathMap species_path_map;
    species_path_map["reference"] = common_args.reference_path;
    species_path_map["query"] = common_args.query_path;

    // 拷贝或下载原始文件（并行执行）
    copyRawData(common_args.work_dir_path, species_path_map, common_args.thread_num);

    // 声明interval文件映射和SeqPro managers变量
    std::map<SpeciesName, FilePath> interval_files_map;
    std::map<SpeciesName, SeqPro::ManagerVariant> seqpro_managers;
    SeqPro::Length reference_min_seq_length = std::numeric_limits<SeqPro::Length>::max();
    // 清洗 FASTA 文件（统一格式，替换非法字符）
    cleanRawDataset(common_args.work_dir_path, species_path_map,
                    common_args.thread_num);
    // ------------------------------
    // 重复遮蔽
    // ------------------------------
    if (common_args.enable_repeat_masking) {
        spdlog::info("Repeat masking enabled. Generating interval files based on "
            "raw files...");

        // 1. 生成interval文件
        interval_files_map = repeatSeqMasking(
            common_args.work_dir_path, species_path_map, common_args.thread_num);

        if (interval_files_map.empty()) {
            // 没成功生成，中止程序建议用户检查一下
            spdlog::error("No interval files were generated. Please check the input "
                "FASTA files and ensure they are valid.");
            return 1;
        }
        spdlog::info("Interval files generated successfully.");
        spdlog::info("Creating SeqPro managers with repeat masking support...");

        for (const auto &[species_name, cleaned_fasta_path]: species_path_map) {
            if (interval_files_map.contains(species_name)) {
                try {
                    auto original_manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);
                    auto manager = std::make_unique<SeqPro::MaskedSequenceManager>(
                        std::move(original_manager),
                        interval_files_map[species_name]
                    );
                    spdlog::info("[{}] SeqPro Manager created with repeat masking: {}",
                                 species_name, cleaned_fasta_path.string());

                    // 记录序列统计信息
                    recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);
                    
                    // 移动到seqpro_managers中
                    seqpro_managers[species_name] = std::move(manager);
                } catch (const std::exception &e) {
                    spdlog::error("[{}] Error creating SeqPro manager: {}", species_name,
                                  e.what());
                    return 1;
                }
            } else {
                try {
                    auto manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);

                    spdlog::info("[{}] SeqPro Manager created without repeat masking: {}", 
                                 species_name, cleaned_fasta_path.string());
                    
                    // 记录序列统计信息
                    recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);

                     seqpro_managers[species_name] = std::move(manager);
                } catch (const std::exception &e) {
                    spdlog::error("[{}] Error creating SeqPro manager: {}", species_name,
                                  e.what());
                    return 1;
                }
            }
        }
    } else {
        // 不开启重复遮蔽，基于清洗后的文件创建常规SeqPro managers
        spdlog::info(
            "Repeat masking disabled. Creating standard SeqPro managers...");

        for (const auto &[species_name, cleaned_fasta_path]: species_path_map) {
            try {
                auto manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);
                spdlog::info("[{}] SeqPro Manager created: {}", species_name,
                             cleaned_fasta_path.string());

                // 记录序列统计信息
                recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);

                seqpro_managers[species_name] = std::move(manager);
            } catch (const std::exception &e) {
                spdlog::error("[{}] Error creating SeqPro manager: {}", species_name,
                              e.what());
                return 1;
            }
        }
    }


    // 可选：按染色体拆分、按 chunk 切割（目前注释掉）
    // splitRawDataToChr(...)
    // splitChrToChunk(...)

    // ------------------------------
    // 初始化比对器
    // ------------------------------
    PairRareAligner pra(
        common_args.work_dir_path,
        common_args.thread_num,
        common_args.chunk_size,
        common_args.overlap_size,
        common_args.min_anchor_length,
        common_args.max_anchor_frequency
    );


    // ------------------------------
    // 步骤 1：构建索引
    // ------------------------------
    auto t_start_build = std::chrono::steady_clock::now();

    pra.buildIndex("reference", seqpro_managers["reference"], false); // 可切换 CaPS / divsufsort
    auto t_end_build = std::chrono::steady_clock::now();
    std::chrono::duration<double> build_time = t_end_build - t_start_build;
    spdlog::info("Index built in {:.3f} seconds.", build_time.count());

    // ------------------------------
    // 步骤 2：查询序列比对
    // ------------------------------
    auto t_start_align = std::chrono::steady_clock::now();
    sdsl::int_vector<0> ref_global_cache;
    // 初始化ref_global_cache：采样策略避免二分搜索
    auto sampling_interval = std::min(static_cast<SeqPro::Length>(32), reference_min_seq_length);
    auto total_length = std::visit([](auto&& manager_ptr) {
        return manager_ptr->getTotalLength();
    }, seqpro_managers["reference"]);
    auto cache_size = (total_length / sampling_interval) + 1;
    ref_global_cache.resize(cache_size);

    // Fill ref_global_cache: pre-calculate sequence IDs for each sampling point
    spdlog::info("Filling ref_global_cache, sampling_interval={}, cache_size={}", sampling_interval, cache_size);
    auto t_start_cache = std::chrono::steady_clock::now();
    
    std::visit([&](auto&& manager_ptr) {
        // 获取所有序列信息，按 global_start_pos 排序
        auto seq_names = manager_ptr->getSequenceNames();
        std::vector<const SeqPro::SequenceInfo*> seq_infos;
        seq_infos.reserve(seq_names.size());
        
        for (const auto& name : seq_names) {
            if constexpr (std::is_same_v<std::decay_t<decltype(manager_ptr)>, std::unique_ptr<SeqPro::SequenceManager>>) {
                const auto* info = manager_ptr->getIndex().getSequenceInfo(name);
                if (info) seq_infos.push_back(info);
            } else if constexpr (std::is_same_v<std::decay_t<decltype(manager_ptr)>, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                const auto* info = manager_ptr->getOriginalManager().getIndex().getSequenceInfo(name);
                if (info) seq_infos.push_back(info);
            }
        }
        
        // 按全局起始位置排序
        std::sort(seq_infos.begin(), seq_infos.end(), 
                  [](const SeqPro::SequenceInfo* a, const SeqPro::SequenceInfo* b) {
                      return a->global_start_pos < b->global_start_pos;
                  });
        
        // 顺序填充
        size_t current_seq_idx = 0;
        for (SeqPro::Position i = 0; i < cache_size; ++i) {
            SeqPro::Position sample_global_pos = i * sampling_interval;
            
            if (sample_global_pos >= total_length) {
                ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
                continue;
            }
            
            // 向前查找包含当前位置的序列
            while (current_seq_idx < seq_infos.size()) {
                const auto* current_seq = seq_infos[current_seq_idx];
                SeqPro::Position seq_end = current_seq->global_start_pos + current_seq->length;
                
                if (sample_global_pos >= current_seq->global_start_pos && sample_global_pos < seq_end) {
                    // 找到了包含该位置的序列
                    ref_global_cache[i] = current_seq->id;
                    break;
                } else if (sample_global_pos >= seq_end) {
                    // 当前序列已经过了，移动到下一个序列
                    current_seq_idx++;
                } else {
                    // sample_global_pos < current_seq->global_start_pos, shouldn't happen
                    spdlog::warn("Unexpected coordinate order: sample_pos={}, seq_start={}", 
                                sample_global_pos, current_seq->global_start_pos);
                    ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
                    break;
                }
            }
            
            // 如果遍历完所有序列都没找到，标记为无效
            if (current_seq_idx >= seq_infos.size()) {
                ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
            }
        }
    }, seqpro_managers["reference"]);
    sdsl::util::bit_compress(ref_global_cache);
    auto t_end_cache = std::chrono::steady_clock::now();
    std::chrono::duration<double> cache_time = t_end_cache - t_start_cache;
    spdlog::info("ref_global_cache filling completed in {:.3f} seconds", cache_time.count());

    MatchVec3DPtr anchors = pra.alignPairGenome(
        "query", seqpro_managers["query"], FAST_SEARCH, false, ref_global_cache,sampling_interval);
    auto t_end_align = std::chrono::steady_clock::now();
    std::chrono::duration<double> align_time = t_end_align - t_start_align;
    spdlog::info("Query aligned in {:.3f} seconds.", align_time.count());


//    // ------------------------------
//    // 临时验证：检查anchors的序列匹配正确性
//    // ------------------------------
//    spdlog::info("开始验证 anchors 结果的正确性…");
//
//    auto reverseComplement = [](const std::string &seq) -> std::string {
//        std::string result;
//        result.reserve(seq.length());
//        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
//            result.push_back(BASE_COMPLEMENT[static_cast<unsigned char>(*it)]);
//        }
//        return result;
//    };
//
//    uint64_t total_matches = 0;
//    uint64_t correct_matches = 0;
//    uint64_t incorrect_matches = 0;
//
//    uint64_t total_work_items = 0;
//    for (const auto &level1: *anchors) {
//        for (const auto &level2: level1) {
//            total_work_items += level2.size();
//        }
//    }
//
//    if (total_work_items == 0) {
//        spdlog::info("没有需要验证的匹配项。");
//        return 0;
//    }
//
//    spdlog::info("总计需要验证 {} 个匹配项。", total_work_items);
//
//    // --- 新增代码：用于进度的原子计数器 ---
//    std::atomic<uint64_t> processed_items = 0;
//    // 动态计算报告间隔，目标是报告大约100次，避免过于频繁或稀疏
//    const uint64_t report_interval = std::max(1ULL, total_work_items / 100ULL);
//    std::atomic<bool> diagnostic_dump_done = false;
//
//
//#pragma omp parallel reduction(+: total_matches, correct_matches, incorrect_matches)
//    {
//#pragma omp for schedule(dynamic) nowait
//        for (size_t i = 0; i < anchors->size(); ++i) {
//            for (size_t j = 0; j < (*anchors)[i].size(); ++j) {
//                for (size_t k = 0; k < (*anchors)[i][j].size(); ++k) {
//                    uint64_t current_processed = processed_items.fetch_add(1, std::memory_order_relaxed) + 1;
//
//                    if (current_processed % report_interval == 0) {
//                        spdlog::info("验证进度: {} / {} ({:.2f}%)",
//                                     current_processed,
//                                     total_work_items,
//                                     (100.0 * current_processed / total_work_items));
//                    }
//
//                    const Match &match = (*anchors)[i][j][k];
//                    ++total_matches;
//
//                    try {
//                        // 现在anchor返回的坐标是正确的局部坐标，可以直接使用
//                        std::string ref_seq = std::visit([&match](auto &&manager_ptr) -> std::string {
//                            using PtrType = std::decay_t<decltype(manager_ptr)>;
//                            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
//                                return manager_ptr->getSubSequence(match.ref_region.chr_name, match.ref_region.start,
//                                                                   match.ref_region.length);
//                            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
//                                SeqPro::MaskedSequenceManager> >) {
//                                return manager_ptr->getSubSequence(match.ref_region.chr_name, match.ref_region.start,
//                                                                   match.ref_region.length);
//                            } else {
//                                throw std::runtime_error("Unhandled manager type in variant.");
//                            }
//                        }, seqpro_managers["reference"]);
//                        std::string query_seq = std::visit([&match](auto &&manager_ptr) -> std::string {
//                            using PtrType = std::decay_t<decltype(manager_ptr)>;
//                            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager> >) {
//                                return manager_ptr->getSubSequence(match.query_region.chr_name,
//                                                                   match.query_region.start, match.query_region.length);
//                            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<
//                                SeqPro::MaskedSequenceManager> >) {
//                                return manager_ptr->getSubSequence(match.query_region.chr_name,
//                                                                   match.query_region.start, match.query_region.length);
//                            } else {
//                                throw std::runtime_error("Unhandled manager type in variant.");
//                            }
//                        }, seqpro_managers["query"]);
//
//                        if (match.strand == REVERSE) {
//                            query_seq = reverseComplement(query_seq);
//                        }
//
//                        if (ref_seq == query_seq) {
//                            ++correct_matches;
//                        } else {
//                            ++incorrect_matches;
//
//                            // 仅当这是第一个被捕获的错误时，才打印详细信息并退出
//                            bool expected = false;
//                            if (!diagnostic_dump_done.load(std::memory_order_relaxed) &&
//                                diagnostic_dump_done.compare_exchange_strong(expected, true)) {
//                                spdlog::warn(
//                                    "序列不匹配: ref_chr={}, ref_start={}, query_chr={}, "
//                                    "query_start={}, length={}, strand={}\n"
//                                    "  Ref Seq:    {}\n"
//                                    "  Query Seq{}: {}",
//                                    match.ref_region.chr_name, match.ref_region.start, match.query_region.chr_name,
//                                    match.query_region.start, match.ref_region.length,
//                                    (match.strand == FORWARD ? "FORWARD" : "REVERSE"), ref_seq,
//                                    (match.strand == REVERSE ? " (RC)" : ""), query_seq
//                                );
//                                spdlog::error("--- [CAPTURED FIRST MISMATCH] INITIATING DIAGNOSTIC DUMP ---");
//
//                                // 打印错误匹配的详细信息
//                                spdlog::error("Failing Match Details:");
//                                spdlog::error("  - Reference: {}:{} (len:{})", match.ref_region.chr_name,
//                                              match.ref_region.start, match.ref_region.length);
//                                spdlog::error("  - Query:     {}:{} (len:{})", match.query_region.chr_name,
//                                              match.query_region.start, match.query_region.length);
//                                spdlog::error("  - Strand:    {}", (match.strand == FORWARD ? "FORWARD" : "REVERSE"));
//
//                                // // 安全退出
//                                // std::exit(1);
//                            }
//                            // --- END: 线程安全的"首错捕获"逻辑 ---
//                        }
//                    } catch (const std::exception &e) {
//                        ++incorrect_matches;
//                        // --- 修改：异常也直接报告 ---
//                        spdlog::warn("处理匹配项时发生异常: {}", e.what());
//                    }
//                }
//            }
//        } // omp for
//    }
//
//    // 确保最终进度是100%
//    spdlog::info("验证进度: {} / {} (100.00%)", total_work_items, total_work_items);
//
//    spdlog::info("验证完成: 总匹配数={}, 正确匹配数={}, 错误匹配数={}, 正确率={:.2f}%",
//                 total_matches,
//                 correct_matches,
//                 incorrect_matches,
//                 total_matches ? (100.0 * correct_matches / total_matches) : 0.0);

    
    RaMesh::RaMeshMultiGenomeGraph graph(seqpro_managers);
    // ------------------------------
    // 步骤 3：过滤锚点
    // ------------------------------
    spdlog::info("Filtering anchors...");
    auto t_start_filer = std::chrono::steady_clock::now();
    MatchClusterVecPtr cluster_vec_ptr = pra.filterPairSpeciesAnchors("query", anchors, seqpro_managers["query"], graph);
    auto t_end_filer = std::chrono::steady_clock::now();
    std::chrono::duration<double> filter_time = t_end_filer - t_start_filer;
    spdlog::info("Anchors clustered in {:.3f} seconds.", filter_time.count());

    pra.constructGraphByGreedy("query", cluster_vec_ptr, graph, 50);

    // ------------------------------
    // 退出
    // ------------------------------
    spdlog::info("RaMA-G exits!");
    return 0;
}
