// RaMA-G.cpp: 定义应用程序的入口点。
//

#include "data_process.h"
#include "config.hpp"
#include "index.h"

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
			}
			else {
				std::filesystem::create_directories(common_args.work_dir_path);
			}
			setupLoggerWithFile(common_args.work_dir_path);
			spdlog::info("RaMA-G version {}", VERSION);
			spdlog::info("Restart mode enabled.");

			FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
			std::ifstream is(config_path);
			if (!is) {
				spdlog::error("Failed to open {} for loading CommonArgs", config_path.string());
				return false;
			}
			cereal::JSONInputArchive archive(is);
			archive(common_args);
			spdlog::info("CommonArgs loaded from {}", config_path.string());
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

			// overlap should be less than chunk_size
			if (common_args.overlap_size >= common_args.chunk_size) {
				throw std::runtime_error("Overlap size must be less than chunk size.");
			}


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

	SpeciesPathMap species_path_map;
	//SpeciesChrPathMap species_chr_path_map;
	//SpeciesChunkInfoMap species_chunk_info_map;
	species_path_map["reference"] = common_args.reference_path;
	species_path_map["query"] = common_args.query_path;

	copyRawData(common_args.work_dir_path, species_path_map, common_args.thread_num);
	cleanRawDataset(common_args.work_dir_path, species_path_map, common_args.thread_num);

	FastaManager ref_fasta_manager(species_path_map["reference"], getFaiIndexPath(species_path_map["reference"]));
	//splitRawDataToChr(common_args.work_dir_path, species_path_map, species_chr_path_map, common_args.thread_num);

	//SpeciesChrPathMap tmp_species_chr_path_map = species_chr_path_map;
	//// remove reference
	//tmp_species_chr_path_map.erase("reference");
	//splitChrToChunk(common_args.work_dir_path, tmp_species_chr_path_map, species_chunk_info_map, common_args.chunk_size, common_args.overlap_size, common_args.thread_num);
	
	IndexManager index_manager(common_args.work_dir_path, common_args.thread_num);
	IndexPathMap index_path_map;
	index_path_map["reference"] = index_manager.buildIndex("reference", ref_fasta_manager, IndexType::RIndexType);
	
	spdlog::info("RaMA-G exits!");

	return 0;
}