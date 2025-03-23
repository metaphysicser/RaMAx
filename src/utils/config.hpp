#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "CLI/CLI.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/async.h"
#include <cereal/archives/json.hpp> 
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <memory>
#include <filesystem>

#define VERSION "1.0.0"

#define LOGGER_NAME "logger"
#define LOGGER_FILE "RaMA-G.log"
#define CONFIG_FILE "config.json"

#define DATA_DIR "data"	
#define RAW_DATA_DIR "raw_data"
#define CLEAN_DATA_DIR "clean_data"
#define SPLIT_CHR_DIR "split_chr"
#define CHUNK_DIR "chunk"
#define CHUNK_MAP_FILE "chunk_map.json"

namespace cereal {
	template <class Archive>
	void save(Archive& ar, const std::filesystem::path& p)
	{
		std::string path_str = p.string();
		ar(path_str);
	}

	template <class Archive>
	void load(Archive& ar, std::filesystem::path& p)
	{
		std::string path_str;
		ar(path_str);
		p = std::filesystem::path(path_str);
	}
} // namespace cereal

// Custom formatter for CLI11, unify option display style
class CustomFormatter : public CLI::Formatter {
public:
	CustomFormatter() : Formatter() {}

	// Display option format as [VALUE], remove "TEXT REQUIRED" style
	std::string make_option_opts(const CLI::Option* opt) const override {
		if (opt->get_type_size() == 0) return "";  // No display for flag options
		std::ostringstream out;
		out << " " << opt->get_type_name();  // Automatically display <path> or <int>
		if (!opt->get_default_str().empty())
			out << " (default: " << opt->get_default_str() << ")";
		return out.str();
	}

	// Customize usage example for better readability
	std::string make_usage(const CLI::App* app, std::string name) const override {
		std::ostringstream out;
		out << "Usage:\n"
			<< "  ./RaMA-G -r <ref.fa> -q <query.fa> -o <out_dir> [options]\n\n"
			<< "Example:\n"
			<< "  ./RaMA-G -r ref.fa -q query.fa -o output/ -t 8\n\n";
		return out.str();
	}
};

// Whitespace trimming validator for CLI11
inline CLI::Validator trim_whitespace = CLI::Validator(
	[](std::string& s) {
		auto start = s.find_first_not_of(" \t\n\r");
		auto end = s.find_last_not_of(" \t\n\r");
		s = (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
		return std::string();  // empty means valid
	}, ""
);

// Enum for output format types
enum class OutputFormat {
	SAM,
	MAF,
	PAF,
	UNKNOWN
};

// Auto detect output format based on file extension
inline OutputFormat detectOutputFormat(const std::filesystem::path& output_path) {
	std::string ext = output_path.extension().string();
	if (ext == ".sam") return OutputFormat::SAM;
	if (ext == ".maf") return OutputFormat::MAF;
	if (ext == ".paf") return OutputFormat::PAF;
	return OutputFormat::UNKNOWN;
}

// Structure to store common command-line arguments
struct CommonArgs {
	std::filesystem::path reference_path = "";
	std::filesystem::path query_path = "";
	std::filesystem::path output_path = "";
	std::filesystem::path work_dir_path = "";

	size_t chunk_size = 10000000;
	size_t overlap_size = 100000;
	bool restart = false;
	int thread_num = std::thread::hardware_concurrency();  // Default to hardware concurrency
	OutputFormat output_format = OutputFormat::UNKNOWN;

	template<class Archive>
	void serialize(Archive& ar) {
		ar(
			CEREAL_NVP(reference_path),
			CEREAL_NVP(query_path),
			CEREAL_NVP(output_path),
			CEREAL_NVP(work_dir_path),
			CEREAL_NVP(chunk_size),
			CEREAL_NVP(overlap_size),
			CEREAL_NVP(restart),
			CEREAL_NVP(thread_num),
			CEREAL_NVP(output_format)
		);
	}
};
//bool saveCommonArgs(const CommonArgs& args, const std::string& filename) {
//	std::ofstream os(filename);
//	if (!os) {
//		spdlog::error("Failed to open {} for saving CommonArgs", filename);
//		return false;
//	}
//	cereal::JSONOutputArchive archive(os);
//	archive(args);
//	spdlog::info("CommonArgs saved to {}", filename);
//	return true;
//}
//
//bool loadCommonArgs(CommonArgs& args, const std::string& filename) {
//	std::ifstream is(filename);
//	if (!is) {
//		spdlog::error("Failed to open {} for loading CommonArgs", filename);
//		return false;
//	}
//	cereal::JSONInputArchive archive(is);
//	archive(args);
//	spdlog::info("CommonArgs loaded from {}", filename);
//	return true;
//}

// Setup CLI11 common options
inline void setupCommonOptions(CLI::App* cmd, CommonArgs& args) {
	auto fmt = std::make_shared<CustomFormatter>();
	fmt->column_width(50);
	cmd->formatter(fmt);

	cmd->set_version_flag("-v,--version", std::string("RaMA-G version ") + VERSION);

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

	auto* threads_opt = cmd->add_option("-t,--threads", args.thread_num,
		"Number of threads to use for parallel processing (default: system cores).")
		->default_val(std::thread::hardware_concurrency())
		->envname("RAMA_G_THREADS")
		->capture_default_str()
		->group("Performance")
		->check(CLI::Range(2, std::numeric_limits<int>::max()))
		->type_name("<int>")->transform(trim_whitespace);

	auto* restart_flag = cmd->add_flag("--restart", args.restart,
		"Restart the alignment process by skipping the existing index files.")
		->group("Performance");

	// Set dependencies and exclusions
	restart_flag->needs(workspace_opt);
	restart_flag->excludes(ref_opt, qry_opt, output_opt, threads_opt);
}

// Setup logger with optional file output
inline void setupLoggerWithFile(std::filesystem::path log_dir) {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");

	std::filesystem::path log_file = log_dir / LOGGER_FILE;
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file.string(), true);
	file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");

	spdlog::sinks_init_list sinks = { console_sink, file_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
	spdlog::set_level(spdlog::level::trace);
	spdlog::flush_every(std::chrono::seconds(3));
}

// Setup console-only logger
inline void setupLogger() {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");

	spdlog::sinks_init_list sinks = { console_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
	spdlog::set_level(spdlog::level::trace);
	spdlog::flush_every(std::chrono::seconds(3));
}

// Get full command-line string for logging and reproducibility
inline std::string getCommandLine(int argc, char** argv) {
	std::ostringstream cmd;
	for (int i = 0; i < argc; ++i) {
		cmd << argv[i];
		if (i != argc - 1) cmd << " ";
	}
	return cmd.str();
}

#endif // CONFIG_HPP
