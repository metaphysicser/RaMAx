#include "data_process.h"

// -----------------------------
// 工具函数：判断字符串是否为 URL
// -----------------------------
bool isUrl(const std::string& path_str) {
    static const std::regex url_pattern(R"(^(https?|ftp)://)");
    return std::regex_search(path_str, url_pattern);
}

// -----------------------------
// 验证 URL 是否可达（通过 CURL 发送 HEAD 请求）
// -----------------------------
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

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L);
    curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);

    res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        curl_easy_cleanup(curl);
        throw std::runtime_error("Failed to verify URL: " + url + std::string(curl_easy_strerror(res)));
    }

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

// -----------------------------
// 验证本地文件是否存在且为常规文件
// -----------------------------
void verifyLocalFile(const FilePath& file_path) {
    if (!std::filesystem::exists(file_path)) {
        throw std::runtime_error("Local file does not exist: " + file_path.string());
    }
    if (!std::filesystem::is_regular_file(file_path)) {
        throw std::runtime_error("Path is not a regular file: " + file_path.string());
    }
    spdlog::info("Verified local file exists: {}", file_path.string());
}

// -----------------------------
// CURL 写入回调函数，用于下载时写入文件流
// -----------------------------
size_t writeData(void* ptr, size_t size, size_t nmemb, void* stream) {
    std::ofstream* output_stream = static_cast<std::ofstream*>(stream);
    size_t written = size * nmemb;
    output_stream->write(static_cast<const char*>(ptr), written);
    return written;
}

// -----------------------------
// 使用 CURL 下载文件到指定路径
// -----------------------------
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

// -----------------------------
// 拷贝本地文件到目标路径
// -----------------------------
void copyLocalFile(const FilePath& source, const FilePath& destination) {
    try {
        std::filesystem::copy_file(source, destination, std::filesystem::copy_options::overwrite_existing);
    }
    catch (const std::filesystem::filesystem_error& error) {
        throw std::runtime_error("Failed to copy file from " + source.string() + " to " + destination.string() + " : " + error.what());
    }
    spdlog::info("Copy local file from {} to {}", source.string(), destination.string());
}

// -----------------------------
// 获取文件扩展名字符串
// -----------------------------
std::string getFileExtension(const FilePath& file_path) {
    return file_path.extension().string();
}

// -----------------------------
// 主流程函数：复制或下载所有原始数据文件
// -----------------------------
bool copyRawData(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num) {
    try {
        // 预验证文件路径或 URL
        for (const auto& [key, path] : species_path_map) {
            if (isUrl(path.string())) {
                verifyUrlReachable(path.string());
            }
            else {
                verifyLocalFile(path);
            }
        }

        // 创建原始数据文件夹
        FilePath raw_data_dir = workdir_path / DATA_DIR / RAW_DATA_DIR;
        std::filesystem::create_directories(raw_data_dir);
        spdlog::info("Created directory: {}", raw_data_dir.string());

        // 使用线程池并发执行下载或复制
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
    }
}

// -----------------------------
// 清洗原始数据：调用 FastaManager 清理并索引
// -----------------------------
bool cleanRawDataset(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num) {
    ThreadPool pool(thread_num);
    for (auto it = species_path_map.begin(); it != species_path_map.end(); ++it) {
        std::string species = it->first;
        std::filesystem::path raw_path = it->second;

        pool.enqueue([species, raw_path, &workdir_path, &species_path_map]() {
            FastaManager fasta_manager(raw_path);
            FilePath out_dir = workdir_path / DATA_DIR / CLEAN_DATA_DIR;
            FilePath cleaned_path = fasta_manager.cleanAndIndexFasta(out_dir, species);
            species_path_map[species] = cleaned_path;
            spdlog::info("Species {} cleaned file path updated to: {}", species, cleaned_path.string());
            });
    }
    pool.waitAllTasksDone();
    return true;
}


// -----------------------------
// 获取人类可读的文件大小字符串（自动转化为 KB / MB / GB）
// -----------------------------
std::string getReadableFileSize(const FilePath& filePath) {
    std::string pathStr = filePath.string();

    // 处理远程 URL 文件
    if (pathStr.rfind("http://", 0) == 0 || pathStr.rfind("https://", 0) == 0) {
        CURL* curl = curl_easy_init();
        if (!curl) {
            spdlog::error("curl_easy_init failed for URL: {}", pathStr);
            return "0 B";
        }
        curl_easy_setopt(curl, CURLOPT_URL, pathStr.c_str());
        curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
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

        // 转换为可读单位
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

        // 处理本地文件
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

// -----------------------------
// 判断文件是否小于指定大小（单位：MB）
// -----------------------------
bool isFileSmallerThan(const FilePath& filePath, size_t maxSizeMB) {
    const std::string pathStr = filePath.string();
    const curl_off_t maxBytes = static_cast<curl_off_t>(maxSizeMB) * 1024 * 1024;

    // 对 URL 文件使用 CURL 获取 Content-Length
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

        // 本地文件直接获取大小
    }
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

// -----------------------------
// 生成临时文件路径（用于处理中间文件）
// 例如 a.fasta -> a_in_process.fasta
// -----------------------------
FilePath getTempFilePath(const FilePath& input_path) {
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

// -----------------------------
// 获取对应的 .fai 索引文件路径
// 例如 test.fasta -> test.fasta.fai
// -----------------------------
FilePath getFaiIndexPath(const FilePath& fasta_path) {
    FilePath fai_path = fasta_path;
    fai_path += ".fai";
    return fai_path;
}
