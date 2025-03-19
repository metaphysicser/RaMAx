// RaMA-G.cpp: 定义应用程序的入口点。
//

#include <iostream>

#include "config.hpp"

int main(int argc, char** argv) {
	// Initialize logger and thread pool
	spdlog::init_thread_pool(8192, 1);
	setupLogger();

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

			spdlog::info("RaMA-G version {}", VERSION);
			spdlog::info("Command: {}", getCommandLine(argc, argv));
			spdlog::info("Restart mode enabled.");
			spdlog::info("Work directory: {}", common_args.work_dir_path.string());
			// Skip input file checks, assume checkpoint handling
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
				spdlog::info("Work directory does not exist, created: {}", common_args.work_dir_path.string());
			}

			// Detect output format from file extension
			common_args.output_format = detectOutputFormat(common_args.output_path);
			if (common_args.output_format == OutputFormat::UNKNOWN) {
				throw CLI::ValidationError("Invalid output file extension. Supported: .sam, .maf, .paf");
			}

			// Logging information for normal alignment mode
			spdlog::info("RaMA-G version {}", VERSION);
			spdlog::info("Command: {}", getCommandLine(argc, argv));
			spdlog::info("Alignment mode enabled.");

			spdlog::info("Reference: {} (size: {})",
				common_args.reference_path.string(),
				getReadableFileSize(common_args.reference_path));

			spdlog::info("Query: {} (size: {})",
				common_args.query_path.string(),
				getReadableFileSize(common_args.query_path));

			spdlog::info("Output: {}", common_args.output_path.string());
			spdlog::info("Work directory: {}", common_args.work_dir_path.string());
			spdlog::info("Threads: {}", common_args.thread_num);
		}
	}
	catch (const CLI::Error& e) {
		// Catch CLI parsing/validation errors and log
		spdlog::error("{}", e.what());
		return 1;
	}

	// TODO: Add alignment execution or restart handling here

	return 0;
}
