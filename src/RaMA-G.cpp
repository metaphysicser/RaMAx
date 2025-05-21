// RaMAx.cpp: 定义应用程序的入口点。
// 主程序：负责解析命令行参数，初始化配置，执行基因组比对流程

#include "data_process.h"
#include "config.hpp"
#include "index.h"
#include "anchor.h"
#include "rare_aligner.h"

// ------------------------------------------------------------------
// CLI11 参数注册：配置常用命令行参数（参考 RaMAx 主程序）
// ------------------------------------------------------------------
inline void setupCommonOptions(CLI::App* cmd, CommonArgs& args) {
	auto fmt = std::make_shared<CustomFormatter>();
	fmt->column_width(50);
	cmd->formatter(fmt);

	cmd->set_version_flag("-v,--version", std::string("RaMAx version ") + VERSION);

	auto* ref_opt = cmd->add_option("-r,--reference", args.reference_path,
		"Path to the reference genome file (FASTA format).")
		->group("Input Files")
		->type_name("<path>")->transform(trim_whitespace);

	auto* qry_opt = cmd->add_option("-q,--query", args.query_path,
		"Path to the query genome file (FASTA format).")
		->group("Input Files")
		->type_name("<path>")->transform(trim_whitespace);

	auto* output_opt = cmd->add_option("-o,--output", args.output_path,
		"Path to save alignment results.")
		->group("Output")->type_name("<path>")
		->transform(trim_whitespace);

	auto* workspace_opt = cmd->add_option("-w,--workdir", args.work_dir_path,
		"Path to the working directory for temporary files.")
		->group("Output")->type_name("<path>")
		->transform(trim_whitespace);

	auto* chunk_size_opt = cmd->add_option("--chunk_size", args.chunk_size,
		"Size of each chunk for parallel processing (default: 10000000).")
		->default_val(10000000)
		->capture_default_str()
		->group("Software Parameters")
		->check(CLI::Range(1000000, std::numeric_limits<int>::max()))
		->type_name("<int>")
		->transform(trim_whitespace);

	auto* overlap_size_opt = cmd->add_option("--overlap_size", args.overlap_size,
		"Size of overlap between chunks (default: 100000).")
		->default_val(100000)
		->capture_default_str()
		->group("Software Parameters")
		->check(CLI::Range(0, std::numeric_limits<int>::max()))
		->type_name("<int>")
		->transform(trim_whitespace);

	auto* min_anchor_length_opt = cmd->add_option(
		"--min_anchor_length", args.min_anchor_length,
		"Minimum anchor length (default: 20).")
		->default_val(20)
		->capture_default_str()
		->group("Software Parameters")
		->check(CLI::Range(1, std::numeric_limits<int>::max()))
		->type_name("<int>")
		->transform(trim_whitespace);

	auto* max_anchor_frequency_opt = cmd->add_option(
		"--max_anchor_frequency", args.max_anchor_frequency,
		"Maximum anchor frequency filter (default: 50).")
		->default_val(50)
		->capture_default_str()
		->group("Software Parameters")
		->check(CLI::Range(0, std::numeric_limits<int>::max()))
		->type_name("<int>")
		->transform(trim_whitespace);


	auto* threads_opt = cmd->add_option("-t,--threads", args.thread_num,
		"Number of threads to use for parallel processing (default: system cores).")
		->default_val(std::thread::hardware_concurrency())
		->envname("RAMA_G_THREADS")
		->capture_default_str()
		->group("Performance")
		->check(CLI::Range(1, std::numeric_limits<int>::max()))
		->type_name("<int>")->transform(trim_whitespace);

	auto* restart_flag = cmd->add_flag("--restart", args.restart,
		"Restart the alignment process by skipping the existing index files.")
		->group("Performance");

	// Set dependencies and exclusions
	restart_flag->needs(workspace_opt);
	restart_flag->excludes(ref_opt,
		qry_opt, output_opt, threads_opt,
		chunk_size_opt, overlap_size_opt,
		min_anchor_length_opt, max_anchor_frequency_opt);
}


int main(int argc, char** argv) {
	// 初始化异步日志线程池（spdlog）
	spdlog::init_thread_pool(8192, 1);  // 日志缓冲区容量 8192，单线程日志写入
	// setupLogger();  // 控制台输出日志（不带文件）

	// 初始化 CLI 命令行应用
	CLI::App app{ "RaMA-G: A High-performance Genome Alignment Tool" };

	CommonArgs common_args;  // 存储用户输入的参数
	setupCommonOptions(&app, common_args);  // 注册参数解析
	CLI11_PARSE(app, argc, argv);  // 开始解析命令行参数

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
			}
			else {
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
			archive(common_args);  // 反序列化参数
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
			}
			else {
				verifyLocalFile(common_args.reference_path);
			}

			// 验证查询文件路径
			std::string qry_str = common_args.query_path.string();
			if (isUrl(qry_str)) {
				verifyUrlReachable(qry_str);
			}
			else {
				verifyLocalFile(common_args.query_path);
			}

			// 自动识别输出格式（支持 .sam / .maf / .paf）
			common_args.output_format = detectOutputFormat(common_args.output_path);
			if (common_args.output_format == OutputFormat::UNKNOWN) {
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
	}
	catch (const std::runtime_error& e) {
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

	// 清洗 FASTA 文件（统一格式，替换非法字符）
	cleanRawDataset(common_args.work_dir_path, species_path_map, common_args.thread_num);

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
	pra.buildIndex("reference", species_path_map["reference"], false);  // 可切换 CaPS / divsufsort
	auto t_end_build = std::chrono::steady_clock::now();
	std::chrono::duration<double> build_time = t_end_build - t_start_build;
	spdlog::info("Index built in {:.3f} seconds.", build_time.count());

	// ------------------------------
	// 步骤 2：查询序列比对
	// ------------------------------
	auto t_start_align = std::chrono::steady_clock::now();
	AnchorPtrListVec anchors = pra.findQueryFileAnchor(
		"query", species_path_map["query"], ACCURATE_SEARCH);
	auto t_end_align = std::chrono::steady_clock::now();
	std::chrono::duration<double> align_time = t_end_align - t_start_align;
	spdlog::info("Query aligned in {:.3f} seconds.", align_time.count());

	// ------------------------------
	// 步骤 3：过滤锚点
	// ------------------------------
	auto t_start_filer = std::chrono::steady_clock::now();
	pra.filterAnchors(anchors);
	auto t_end_filer = std::chrono::steady_clock::now();
	std::chrono::duration<double> filter_time = t_end_filer - t_start_filer;
	spdlog::info("Anchors filtered in {:.3f} seconds.", filter_time.count());

	// ------------------------------
	// 退出
	// ------------------------------
	spdlog::info("RaMA-G exits!");
	return 0;
}
