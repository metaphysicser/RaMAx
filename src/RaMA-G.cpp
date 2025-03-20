// RaMA-G.cpp: 定义应用程序的入口点。
//

#include "raw_data_process.hpp"
#include "config.hpp"

int main(int argc, char** argv) {
	// Initialize logger and thread pool
	spdlog::init_thread_pool(8192, 1);
	// setupLogger();

	CLI::App app{ "RaMA-G: A High-performance Genome Alignment Tool" };

	CommonArgs common_args;
	setupCommonOptions(&app, common_args);
	CLI11_PARSE(app, argc, argv);

	try {
		// === Restart Mode Logic ===
		if (common_args.restart) {
			if (common_args.work_dir_path.empty()) {
				throw CLI::RequiredError("In restart mode, --workdir (-w) is required.");
			}

			// Check or create work directory (must be empty if exists)
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
			setupLoggerWithFile(common_args.work_dir_path);
			spdlog::info("RaMA-G version {}", VERSION);
			spdlog::info("Restart mode enabled.");
		}
		// === Normal Alignment Mode Logic ===
		else {
			// Validate required arguments
			if (common_args.reference_path.empty())
				throw CLI::RequiredError("Missing required option: --reference (-r)");
			if (common_args.query_path.empty())
				throw CLI::RequiredError("Missing required option: --query (-q)");
			if (common_args.output_path.empty())
				throw CLI::RequiredError("Missing required option: --output (-o)");
			if (common_args.work_dir_path.empty())
				throw CLI::RequiredError("Missing required option: --workdir (-w)");
#ifndef _DEBUG_
			// Check or create work directory (must be empty if exists)
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
			// Logging information for normal alignment mode
			spdlog::info("RaMA-G version {}", VERSION);
			spdlog::info("Alignment mode enabled.");

			// 针对参考序列路径 (-r) 的验证
			std::string ref_str = common_args.reference_path.string();
			if (isUrl(ref_str)) {
				// 如果是 URL，检查是否可下载
				verifyUrlReachable(ref_str);
			}
			else {
				// 如果是本地文件，检查文件是否存在且为常规文件
				verifyLocalFile(common_args.reference_path);
			}

			// 针对查询序列路径 (-q) 的验证
			std::string qry_str = common_args.query_path.string();
			if (isUrl(qry_str)) {
				verifyUrlReachable(qry_str);
			}
			else {
				verifyLocalFile(common_args.query_path);
			}

			// Detect output format from file extension
			common_args.output_format = detectOutputFormat(common_args.output_path);
			if (common_args.output_format == OutputFormat::UNKNOWN) {
				throw std::runtime_error("Invalid output file extension. Supported: .sam, .maf, .paf");
			}
		}
	}
	catch (const std::runtime_error& e) {
		// Catch CLI parsing/validation errors and log
		spdlog::error("{}", e.what());
		spdlog::error("Use --help for usage information.");
		spdlog::error("Exiting with error code 1.");
		return 1;
	}

	// TODO: Add alignment execution or restart handling here
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

	std::unordered_map<std::string, std::filesystem::path> species_path_map;
	species_path_map["reference"] = common_args.reference_path;
	species_path_map["query"] = common_args.query_path;

	copyRawData(common_args.work_dir_path, species_path_map, common_args.thread_num);

	std::filesystem::path ref_clean_data_path = cleanRawData(common_args.work_dir_path, species_path_map["reference"]);
	std::filesystem::path qry_clean_data_path = cleanRawData(common_args.work_dir_path, species_path_map["query"]);

	return 0;
}