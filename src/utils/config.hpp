#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "CLI/CLI.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/async.h"
#include <iostream>
#include <string>
#include <thread>
#include <memory>
#include <filesystem>

#define LOGGER_NAME "logger"
#define LOGGER_FILE "RaMA-G.log"

#define VERSION "1.0.0"

// 自定义 Formatter，让 option 格式统一为 [OPTION]，更美观
class CustomFormatter : public CLI::Formatter {
public:
	CustomFormatter() : Formatter() {}

	// 统一每个 option 参数格式为 [VALUE]，去掉 CLI11 自带的 "TEXT REQUIRED"
	std::string make_option_opts(const CLI::Option* opt) const override {
		if (opt->get_type_size() == 0) return "";  // flag 类型不显示

		std::ostringstream out;
		out << " " << opt->get_type_name();  // 自动读取 <path> 或 <int>

		if (!opt->get_default_str().empty())
			out << " (default: " << opt->get_default_str() << ")";
		//if (!opt->get_envname().empty())
		//	out << " (env: " << opt->get_envname() << ")";
		return out.str();
	}


	// 去掉默认的 executable 路径 usage 输出（美观）
	std::string make_usage(const CLI::App* app, std::string name) const override {
		std::ostringstream out;
		out << "Usage:\n"
			<< "  ./RaMA-G -r <ref.fa> -q <query.fa> -o <out_dir> [options]\n\n"
			<< "Example:\n"
			<< "  ./RaMA-G -r ref.fa -q query.fa -o output/ -t 8\n\n";
		return out.str();
	}

};
// Structure for common command-line arguments
struct CommonArgs {
	std::filesystem::path reference_path = "";
	std::filesystem::path query_path = "";
	std::filesystem::path output_path = "";

	int thread_num = std::thread::hardware_concurrency();  // Default to system hardware concurrency
};

void addCommonOptions(CLI::App* cmd, CommonArgs& args) {
	auto fmt = std::make_shared<CustomFormatter>();
	fmt->column_width(50);
	cmd->formatter(fmt);

	cmd->set_version_flag("-v,--version", std::string("RaMA-G version ") + VERSION);

	cmd->add_option("-r,--reference", args.reference_path,
		"Path to the reference genome file (FASTA format).")
		->required()->group("Input Files")
		->check(CLI::ExistingFile)
		->type_name("<path>");

	cmd->add_option("-q,--query", args.query_path,
		"Path to the query genome file (FASTA format).")
		->required()->group("Input Files")
		->check(CLI::ExistingFile)
		->type_name("<path>");

	cmd->add_option("-o,--output", args.output_path,
		"Path to save alignment results.")
		->required()->group("Output")
		->type_name("<path>");

	cmd->add_option("-t,--threads", args.thread_num,
		"Number of threads to use for parallel processing (default: system cores).")
		->default_val(std::thread::hardware_concurrency())
		->envname("RAMA_G_THREADS")
		->capture_default_str()
		->group("Performance")
		->check(CLI::Range(2, std::numeric_limits<int>::max()))
		->type_name("<int>");

}


// Function to set up the logger
void setupLoggerWithFile(std::filesystem::path log_dir) {
	// Create console log sink (supports multi-threading + colored output)
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

	// Set the format of console log, including line color
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [thread %t] [%l] %v%$");

	// Define the path for the log file
	std::filesystem::path log_file = log_dir / LOGGER_FILE;

	// Create file log sink (without color)
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file.string(), true);
	file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [thread %t] [%l] %v");

	// Create multi-target logger
	spdlog::sinks_init_list sinks = { console_sink, file_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	// Set the logger as default
	spdlog::set_default_logger(logger);

	// Set global log level
	// TODO change to info
	spdlog::set_level(spdlog::level::trace);

	// Periodically flush logs every 3 seconds
	spdlog::flush_every(std::chrono::seconds(3));
}

// Function to set up the logger (console only)
void setupLogger() {
	// Create console log sink (multi-thread + colored output)
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

	// Set the format of console log
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [thread %t] [%l] %v%$");

	// Create async logger with only console sink
	spdlog::sinks_init_list sinks = { console_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	// Set as default logger
	spdlog::set_default_logger(logger);

	// Set global log level (你可以改成 info）
	spdlog::set_level(spdlog::level::trace);

	// Flush every 3 seconds
	spdlog::flush_every(std::chrono::seconds(3));
}

// Function to trim whitespace and newline characters from the start and end of a string
std::string trim(const std::string& str) {
	auto start = str.find_first_not_of(" \t\r\n");
	auto end = str.find_last_not_of(" \t\r\n");
	return (start == std::string::npos || end == std::string::npos) ? "" : str.substr(start, end - start + 1);
}



#endif // !CONFIF_HPP