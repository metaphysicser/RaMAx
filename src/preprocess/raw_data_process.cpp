#ifndef RAW_DATA_PROCESS_HPP
#define RAW_DATA_PROCESS_HPP

#include "data_process.h"


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
void verifyLocalFile(const FilePath& file_path) {
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
void downloadFile(const std::string& url, const FilePath& destination) {
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
void copyLocalFile(const FilePath& source, const FilePath& destination) {
	try {
		std::filesystem::copy_file(source, destination, std::filesystem::copy_options::overwrite_existing);
	}
	catch (const std::filesystem::filesystem_error& error) {
		throw std::runtime_error("Failed to copy file from " + source.string() + " to " + destination.string() + " : " + error.what());
	}
	spdlog::info("Copy local file from {} to {}", source.string(), destination.string());
}

// Get the file extension from a path
std::string getFileExtension(const FilePath& file_path) {
	return file_path.extension().string();
}

// Main function to copy or download raw data
bool copyRawData(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num) {
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
		FilePath raw_data_dir = workdir_path / DATA_DIR / RAW_DATA_DIR;
		std::filesystem::create_directories(raw_data_dir);
		spdlog::info("Created directory: {}", raw_data_dir.string());

		// Use ThreadPool for parallel processing
		ThreadPool pool(thread_num);
		for (auto it = species_path_map.begin(); it != species_path_map.end(); ++it) {
			const std::string& key = it->first;
			const FilePath& path = it->second;
			std::string extension = getFileExtension(path);
			std::string final_name = key + extension;
			FilePath final_dest = raw_data_dir / final_name;
			pool.enqueue([key, path, extension, raw_data_dir, final_dest]() {
				if (std::filesystem::exists(final_dest)) {
					spdlog::warn("File already exists, skipping: {}", final_dest.string());
					return;
				}

				if (isUrl(path.string())) {
					FilePath temp_dest = raw_data_dir / (key + "_in_download" + extension);
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

bool cleanRawDataset(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num) {
	ThreadPool pool(thread_num);

	// Iterate over the species_path_map using iterator form
	for (auto it = species_path_map.begin(); it != species_path_map.end(); ++it) {
		// Capture the species name and raw file path
		std::string species = it->first;
		std::filesystem::path raw_path = it->second;

		// Enqueue a cleaning task for each species file
		pool.enqueue([species, raw_path, &workdir_path, &species_path_map]() {
			// Clean the raw data file for this species
			// std::filesystem::path cleaned_path = cleanRawData(workdir_path, raw_path);
			FastaManager fasta_manager(raw_path);
			FilePath out_dir = workdir_path / DATA_DIR / CLEAN_DATA_DIR;
			FilePath cleaned_path = fasta_manager.cleanAndIndexFasta(out_dir, species);

			// Update the species_path_map safely
			species_path_map[species] = cleaned_path;

			spdlog::info("Species {} cleaned file path updated to: {}", species, cleaned_path.string());
			});
	}

	// Wait for all tasks to finish before returning
	pool.waitAllTasksDone();
	return true;
}

// Get human-readable file size (auto convert to KB/MB/GB)
std::string getReadableFileSize(const FilePath& filePath) {
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

		curl_off_t content_length = 0;
		res = curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &content_length);

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

bool isFileSmallerThan(const FilePath& filePath, size_t maxSizeMB = 1024) {
	const std::string pathStr = filePath.string();
	// 将阈值转换为字节
	const curl_off_t maxBytes = static_cast<curl_off_t>(maxSizeMB) * 1024 * 1024;

	// 1) 如果是 HTTP/HTTPS URL
	if (pathStr.rfind("http://", 0) == 0 || pathStr.rfind("https://", 0) == 0) {
		CURL* curl = curl_easy_init();
		if (!curl) {
			spdlog::error("curl_easy_init failed for URL: {}", pathStr);
			return false;
		}
		curl_easy_setopt(curl, CURLOPT_URL, pathStr.c_str());
		curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
		curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
		curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L);

		if (curl_easy_perform(curl) != CURLE_OK) {
			spdlog::error("curl_easy_perform failed for URL: {}", pathStr);
			curl_easy_cleanup(curl);
			return false;
		}

		curl_off_t content_length = 0;
		if (curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &content_length) != CURLE_OK
			|| content_length < 0) {
			spdlog::error("Failed to get Content-Length for URL: {}", pathStr);
			curl_easy_cleanup(curl);
			return false;
		}
		curl_easy_cleanup(curl);

		return content_length < maxBytes;
	}
	// 2) 本地文件
	else {
		try {
			auto size = static_cast<curl_off_t>(std::filesystem::file_size(filePath));
			return size < maxBytes;
		}
		catch (const std::filesystem::filesystem_error& e) {
			spdlog::error("Failed to get file size for {}: {}", pathStr, e.what());
			return false;
		}
	}
}

FilePath getTempFilePath(const FilePath& input_path) {
	// Construct a temporary file path by appending '_in_process' to the file name
	FilePath temp_file;
	if (input_path.has_extension()) {
		std::string stem = input_path.stem().string();
		std::string extension = input_path.extension().string();
		temp_file = input_path.parent_path() / (stem + "_in_process" + extension);
	}
	else {
		temp_file = input_path.parent_path() / (input_path.filename().string() + "_in_process");
	}
	return temp_file;
}

FilePath getFaiIndexPath(const FilePath& fasta_path) {
	// 输入FilePath, 返回这个路径加上.fai的FilePath类型
	FilePath fai_path = fasta_path;
	fai_path += ".fai";
	return fai_path;
	
}

#endif