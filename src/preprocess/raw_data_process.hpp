#ifndef RAW_DATA_PROCESS_HPP
#define RAW_DATA_PROCESS_HPP

#include <regex>
#include <filesystem>
#include <unordered_map>
#include <string>
#include <fstream>
#include <stdexcept>
#include <spdlog/spdlog.h>
#include <curl/curl.h>
#include <zlib.h>
#include <cctype>
#include "config.hpp"
#include "threadpool.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

// Check if a string is a URL
bool isUrl(const std::string& path_str) {
	static const std::regex url_pattern(R"(^(https?|ftp)://)");
	return std::regex_search(path_str, url_pattern);
}

// Verify if a URL is reachable
bool verifyUrlReachable(const std::string& url) {
	if (url.empty()) {
		throw std::runtime_error("URL is empty and cannot be reached: " + url);
	}

	CURL* curl = curl_easy_init();
	if (!curl) {
		throw std::runtime_error("Failed to initialize CURL for URL verification" + url);
	}

	CURLcode res;
	long response_code = 0;

	// Set CURL options for URL verification
	curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
	curl_easy_setopt(curl, CURLOPT_NOBODY, 1L); // Perform a HEAD request
	curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L); // Timeout after 10 seconds
	curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);

	// Perform the request
	res = curl_easy_perform(curl);
	if (res != CURLE_OK) {
		curl_easy_cleanup(curl);
		throw std::runtime_error("Failed to verify URL: " + url + std::string(curl_easy_strerror(res)));
	}

	// Get the HTTP response code
	curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);
	curl_easy_cleanup(curl);

	if (response_code >= 200 && response_code < 400) {
		spdlog::info("Verified URL is reachable: {} (HTTP {})", url, response_code);
		return true;
	}
	else {
		throw std::runtime_error("URL verification failed with HTTP code: " + std::to_string(response_code));
	}
}

// Verify if a local file exists and is valid
void verifyLocalFile(const std::filesystem::path& file_path) {
	if (!std::filesystem::exists(file_path)) {
		throw std::runtime_error("Local file does not exist: " + file_path.string());
	}
	if (!std::filesystem::is_regular_file(file_path)) {
		throw std::runtime_error("Path is not a regular file: " + file_path.string());
	}
	spdlog::info("Verified local file exists: {}", file_path.string());
}

// Callback function for CURL write operation
size_t writeData(void* ptr, size_t size, size_t nmemb, void* stream) {
	std::ofstream* output_stream = static_cast<std::ofstream*>(stream);
	size_t written = size * nmemb;
	output_stream->write(static_cast<const char*>(ptr), written);
	return written;
}

// Download a file from a URL using libcurl
void downloadFile(const std::string& url, const std::filesystem::path& destination) {
	spdlog::info("Downloading {} to {}", url, destination.string());
	CURL* curl;
	CURLcode result;
	std::ofstream output_stream(destination, std::ios::binary);

	if (!output_stream) {
		throw std::runtime_error("Failed to create file at: " + destination.string());
	}

	curl = curl_easy_init();
	if (!curl) {
		throw std::runtime_error("Failed to initialize CURL");
	}

	curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeData);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, &output_stream);
	curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);

	result = curl_easy_perform(curl);
	if (result != CURLE_OK) {
		curl_easy_cleanup(curl);
		throw std::runtime_error("CURL failed: " + std::string(curl_easy_strerror(result)));
	}

	curl_easy_cleanup(curl);
	output_stream.close();
	spdlog::info("Download completed: {}", destination.string());
}

// Copy a local file to the destination
void copyLocalFile(const std::filesystem::path& source, const std::filesystem::path& destination) {
	try {
		std::filesystem::copy_file(source, destination, std::filesystem::copy_options::overwrite_existing);
	}
	catch (const std::filesystem::filesystem_error& error) {
		throw std::runtime_error("Failed to copy file from " + source.string() + " to " + destination.string() + " : " + error.what());
	}
	spdlog::info("Copy local file from {} to {}", source.string(), destination.string());
}

// Get the file extension from a path
std::string getFileExtension(const std::filesystem::path& file_path) {
	return file_path.extension().string();
}

// Main function to copy or download raw data
bool copyRawData(const std::filesystem::path workdir_path, std::unordered_map<std::string, std::filesystem::path>& species_path_map, int thread_num) {
	try {
		// Validate paths
		for (const auto& [key, path] : species_path_map) {
			if (isUrl(path.string())) {
				verifyUrlReachable(path.string());
			}
			else {
				verifyLocalFile(path);
			}
		}

		// Create directories
		std::filesystem::path raw_data_dir = workdir_path / DATA_DIR / RAW_DATA_DIR;
		std::filesystem::create_directories(raw_data_dir);
		spdlog::info("Created directory: {}", raw_data_dir.string());

		// Use ThreadPool for parallel processing
		ThreadPool pool(thread_num);
		for (auto it = species_path_map.begin(); it != species_path_map.end(); ++it) {
			const std::string& key = it->first;
			const std::filesystem::path& path = it->second;
			std::string extension = getFileExtension(path);
			std::string final_name = key + extension;
			std::filesystem::path final_dest = raw_data_dir / final_name;
			pool.enqueue([key, path, extension, raw_data_dir, final_dest]() {
				if (std::filesystem::exists(final_dest)) {
					spdlog::warn("File already exists, skipping: {}", final_dest.string());
					return;
				}

				if (isUrl(path.string())) {
					std::filesystem::path temp_dest = raw_data_dir / (key + "_in_download" + extension);
					downloadFile(path.string(), temp_dest);
					std::filesystem::rename(temp_dest, final_dest);
				}
				else {
					copyLocalFile(path, final_dest);
				}

				spdlog::info("Successfully processed species {} -> {}", key, final_dest.string());
				});
			species_path_map[key] = final_dest;
		}

		pool.waitAllTasksDone();
		spdlog::info("All raw sequence data copy to work directory successfully.");
		return true;
	}
	catch (const std::exception& error) {
		spdlog::error("Error occurred: {}", error.what());
		throw std::runtime_error("Failed to copy raw data to work directory.");
		return false;
	}
}

std::filesystem::path cleanRawData(const std::filesystem::path workdir_path, const std::filesystem::path raw_data_path) {
	// 构造输出目录： workdir_path / DATA_DIR / CLEAN_DATA
	std::filesystem::path out_dir = workdir_path / DATA_DIR / CLEAN_DATA;
	if (!std::filesystem::exists(out_dir)) {
		std::filesystem::create_directories(out_dir);
		spdlog::info("Created directory: {}", out_dir.string());
	}

	// 构造目标文件路径：保持原有文件名
	std::filesystem::path target_file = out_dir / raw_data_path.filename();
	if (std::filesystem::exists(target_file)) {
		spdlog::info("Cleaned file already exists, skipping: {}", target_file.string());
		return target_file;
	}

	// 构造临时文件路径（在文件名后加 _in_process）
	std::filesystem::path temp_file;
	if (target_file.has_extension()) {
		std::string stem = target_file.stem().string();
		std::string extension = target_file.extension().string();
		temp_file = target_file.parent_path() / (stem + "_in_process" + extension);
	}
	else {
		temp_file = target_file.parent_path() / (target_file.filename().string() + "_in_process");
	}

	// 如果临时文件已存在，则删除（重新开始清理）
	if (std::filesystem::exists(temp_file)) {
		spdlog::warn("Temporary file exists, removing: {}", temp_file.string());
		std::filesystem::remove(temp_file);
	}

	// 使用 kseq.h 进行流式读取（支持 .gz 格式）
	gzFile fp = gzopen(raw_data_path.string().c_str(), "r");
	if (!fp) {
		throw std::runtime_error("Failed to open raw data file: " + raw_data_path.string());
	}
	kseq_t* seq = kseq_init(fp);
	int64_t l;

	// 打开临时文件进行输出
	std::ofstream ofs(temp_file, std::ios::out);
	if (!ofs) {
		kseq_destroy(seq);
		gzclose(fp);
		throw std::runtime_error("Failed to open temporary file for writing: " + temp_file.string());
	}

	// 流式读取每条 FASTA 记录，逐条处理
	while ((l = kseq_read(seq)) >= 0) {
		// 清理序列：将所有字符转换为大写，若不属于 A、G、C、T、N 则替换为 N
		std::string cleaned;
		cleaned.reserve(seq->seq.l);
		for (int i = 0; i < seq->seq.l; i++) {
			char c = seq->seq.s[i];
			char uc = std::toupper(static_cast<unsigned char>(c));
			if (uc == 'A' || uc == 'G' || uc == 'C' || uc == 'T' || uc == 'N') {
				cleaned.push_back(uc);
			}
			else {
				cleaned.push_back('N');
			}
		}
		// 写入 FASTA 格式记录：标题行和清理后的序列
		ofs << ">" << seq->name.s;
		ofs << "\n" << cleaned << "\n";
	}

	// 关闭输出流和 kseq
	ofs.close();
	kseq_destroy(seq);
	gzclose(fp);

	// 清理完成后，将临时文件重命名为目标文件（满足断点续传）
	std::filesystem::rename(temp_file, target_file);
	spdlog::info("Cleaned raw data saved to: {}", target_file.string());
	return target_file;
}

//
//// Function to split genome data into chromosome-level files
//bool splitRawData(const std::filesystem::path workdir_path,
//	std::unordered_map<std::string, std::filesystem::path>& species_path_map,
//	std::unordered_map<std::string, std::unordered_map<std::string, std::filesystem::path>>& species_chr_path_map,
//	int thread_num) {
//	try {
//		// Create "split_chr" directory if it doesn't exist
//		std::filesystem::path split_chr_dir = workdir_path / DATA_DIR / SPLIT_CHR_DIR;
//		if (!std::filesystem::exists(split_chr_dir)) {
//			std::filesystem::create_directories(split_chr_dir);
//			spdlog::info("Created split_chr directory: {}", split_chr_dir.string());
//		}
//		else {
//			spdlog::warn("split_chr directory already exists: {}", split_chr_dir.string());
//		}
//
//		// Create a thread pool for parallel processing of multiple species
//		ThreadPool pool(thread_num);
//
//		// Loop through each species in the species_path_map
//		for (auto it = species_path_map.begin(); it != species_path_map.end(); ++it) {
//			const std::string& species = it->first;  // Species name
//			const std::filesystem::path& fasta_path = it->second; // Path to the genome file (FASTA format)
//
//			// Enqueue a task to process each species
//			pool.enqueue([&species, &fasta_path, &split_chr_dir, &species_chr_path_map]() {
//				try {
//					// Create a directory for the current species
//					std::filesystem::path species_dir = split_chr_dir / species;
//					if (!std::filesystem::exists(species_dir)) {
//						std::filesystem::create_directories(species_dir);
//						spdlog::info("Created directory for species: {}", species_dir.string());
//					}
//					else {
//						spdlog::warn("Directory for species {} already exists: {}", species, species_dir.string());
//					}
//
//					// Map to store the chromosome paths for the species
//					std::unordered_map<std::string, std::filesystem::path> chr_path_map;
//
//					// Open the FASTA file (supporting both .gz and .fasta formats)
//					gzFile fp = gzopen(fasta_path.string().c_str(), "r");
//					if (!fp) {
//						throw std::runtime_error("Failed to open FASTA file: " + fasta_path.string());
//					}
//
//					// Initialize kseq to read the FASTA file
//					kseq_t* seq = kseq_init(fp);
//					int ret;
//					int count = 0;
//
//					// Read sequences from the FASTA file
//					while ((ret = kseq_read(seq)) >= 0) {
//						std::string header(seq->name.s);  // Chromosome header (e.g., "chr1")
//						std::string sequence(seq->seq.s); // Sequence data
//
//						// Path to save the chromosome sequence file
//						std::filesystem::path chr_file = species_dir / (header + ".fasta");
//						count++; // Count the number of chromosomes processed
//
//						// Record the chromosome file path in the map
//						chr_path_map[header] = chr_file;
//
//						// Skip if the chromosome file already exists
//						if (std::filesystem::exists(chr_file)) {
//							spdlog::warn("File already exists, skipping: {}", chr_file.string());
//							continue;
//						}
//
//						// Write the chromosome sequence to the file
//						std::ofstream chr_ofs(chr_file, std::ios::out);
//						chr_ofs << ">" << header << std::endl;
//						chr_ofs << sequence << std::endl;
//						chr_ofs.close();
//					}
//
//					// Store the chromosome paths for the species in the species_chr_path_map
//					species_chr_path_map[species] = chr_path_map;
//					spdlog::info("Processed {} chromosomes for species {}, saved to: {}", count, species, species_dir.string());
//				}
//				catch (const std::exception& e) {
//					spdlog::error("Error processing species {}: {}", species, e.what());
//					throw std::runtime_error("Failed to split genome data for species: " + species);
//				}
//				});
//		}
//
//		// Wait for all tasks in the thread pool to finish
//		pool.waitAllTasksDone();
//		spdlog::info("Successfully split genomes into chromosomes!");
//		return true;
//	}
//	catch (const std::exception& e) {
//		spdlog::error("Error in splitRawData: {}", e.what());
//		throw std::runtime_error("Failed to split genome data into chromosomes.");
//		return false;
//	}
//}

// Get human-readable file size (auto convert to KB/MB/GB)
inline std::string getReadableFileSize(const std::filesystem::path& filePath) {
	std::string pathStr = filePath.string();
	// 检查路径是否为 URL
	if (pathStr.rfind("http://", 0) == 0 || pathStr.rfind("https://", 0) == 0) {
		CURL* curl = curl_easy_init();
		if (!curl) {
			spdlog::error("curl_easy_init failed for URL: {}", pathStr);
			return "0 B";
		}
		curl_easy_setopt(curl, CURLOPT_URL, pathStr.c_str());
		curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);            // 使用 HEAD 请求
		curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
		curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L);

		CURLcode res = curl_easy_perform(curl);
		if (res != CURLE_OK) {
			spdlog::error("curl_easy_perform failed for URL {}: {}", pathStr, curl_easy_strerror(res));
			curl_easy_cleanup(curl);
			return "0 B";
		}

		double content_length = 0;
		res = curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD, &content_length);
		curl_easy_cleanup(curl);
		if (res != CURLE_OK || content_length < 0) {
			spdlog::error("Failed to get Content-Length for URL: {}", pathStr);
			return "0 B";
		}
		size_t size = static_cast<size_t>(content_length);

		const char* units[] = { "B", "KB", "MB", "GB", "TB" };
		int unit_index = 0;
		double display_size = static_cast<double>(size);
		while (display_size >= 1024 && unit_index < 4) {
			display_size /= 1024;
			++unit_index;
		}
		char buf[64];
		snprintf(buf, sizeof(buf), "%.2f %s", display_size, units[unit_index]);
		return std::string(buf);
	}
	else {
		try {
			auto size = std::filesystem::file_size(filePath);
			const char* units[] = { "B", "KB", "MB", "GB", "TB" };
			int unit_index = 0;
			double display_size = static_cast<double>(size);
			while (display_size >= 1024 && unit_index < 4) {
				display_size /= 1024;
				++unit_index;
			}
			char buf[64];
			snprintf(buf, sizeof(buf), "%.2f %s", display_size, units[unit_index]);
			return std::string(buf);
		}
		catch (const std::filesystem::filesystem_error& e) {
			spdlog::error("Failed to get file size for {}: {}", filePath.string(), e.what());
			return "0 B";
		}
	}
}

#endif