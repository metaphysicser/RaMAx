#include "SeqPro.h"
#include "data_process.h"

// -----------------------------
// 查找 windowmasker 可执行文件路径
// -----------------------------
std::string findWindowmaskerPath() {
    // 1. 首先尝试在系统 PATH 中查找
    if (std::system("which windowmasker > /dev/null 2>&1") == 0) {
        spdlog::info("Found windowmasker in system PATH");
        return "windowmasker";
    }

    // 2. 尝试相对于当前可执行文件的路径（安装后的情况）
    try {
        std::filesystem::path exe_path = std::filesystem::canonical("/proc/self/exe");
        std::filesystem::path exe_dir = exe_path.parent_path();
        std::filesystem::path windowmasker_path = exe_dir / "windowmasker";

        if (std::filesystem::exists(windowmasker_path) &&
            std::filesystem::is_regular_file(windowmasker_path)) {
            spdlog::info("Found windowmasker relative to executable: {}", windowmasker_path.string());
            return windowmasker_path.string();
        }
    }
    catch (const std::exception& e) {
        spdlog::warn("Failed to determine executable path: {}", e.what());
    }

    // 3. 尝试项目开发目录的 bin 路径
    if (std::filesystem::exists("bin/windowmasker") &&
        std::filesystem::is_regular_file("bin/windowmasker")) {
        spdlog::info("Found windowmasker in development bin directory");
        return "bin/windowmasker";
    }

    // 4. 如果都找不到，抛出异常
    throw std::runtime_error("windowmasker executable not found. Please ensure it is installed or available in PATH.");
}

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
        throw std::runtime_error("Failed to initialize CURL for URL verification" +
            url);
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
        throw std::runtime_error("Failed to verify URL: " + url +
            std::string(curl_easy_strerror(res)));
    }

    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);
    curl_easy_cleanup(curl);

    if (response_code >= 200 && response_code < 400) {
        spdlog::info("Verified URL is reachable: {} (HTTP {})", url, response_code);
        return true;
    }
    else {
        throw std::runtime_error("URL verification failed with HTTP code: " +
            std::to_string(response_code));
    }
}

// -----------------------------
// 验证本地文件是否存在且为常规文件
// -----------------------------
void verifyLocalFile(const FilePath& file_path) {
    if (!std::filesystem::exists(file_path)) {
        throw std::runtime_error("Local file does not exist: " +
            file_path.string());
    }
    if (!std::filesystem::is_regular_file(file_path)) {
        throw std::runtime_error("Path is not a regular file: " +
            file_path.string());
    }
    spdlog::info("Verified local file exists: {}", file_path.string());
}

// ------------------------------
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
        throw std::runtime_error("Failed to create file at: " +
            destination.string());
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
        throw std::runtime_error("CURL failed: " +
            std::string(curl_easy_strerror(result)));
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
        std::filesystem::copy_file(
            source, destination, std::filesystem::copy_options::overwrite_existing);
    }
    catch (const std::filesystem::filesystem_error& error) {
        throw std::runtime_error("Failed to copy file from " + source.string() +
            " to " + destination.string() + " : " +
            error.what());
    }
    spdlog::info("Copy local file from {} to {}", source.string(),
        destination.string());
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
bool copyRawData(const FilePath workdir_path, SpeciesPathMap& species_path_map,
    int thread_num) {
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
        for (auto it = species_path_map.begin(); it != species_path_map.end();
            ++it) {
            const std::string& key = it->first;
            const FilePath& path = it->second;
            std::string extension = getFileExtension(path);
            std::string final_name = key + extension;
            FilePath final_dest = raw_data_dir / final_name;

            pool.enqueue([key, path, extension, raw_data_dir, final_dest]() {
                if (std::filesystem::exists(final_dest)) {
                    spdlog::warn("File already exists, skipping: {}",
                        final_dest.string());
                    return;
                }

                if (isUrl(path.string())) {
                    FilePath temp_dest =
                        raw_data_dir / (key + "_in_download" + extension);
                    downloadFile(path.string(), temp_dest);
                    std::filesystem::rename(temp_dest, final_dest);
                }
                else {
                    copyLocalFile(path, final_dest);
                }
                spdlog::info("Successfully processed species {} -> {}", key,
                    final_dest.string());
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
// 清洗原始数据：调用 SeqPro 清理并索引
// -----------------------------
bool cleanRawDataset(const FilePath workdir_path,
    SpeciesPathMap& species_path_map, int thread_num) {
    try {
        // 确保清洗数据输出目录存在
        FilePath clean_data_dir = workdir_path / DATA_DIR / CLEAN_DATA_DIR;
        if (!std::filesystem::exists(clean_data_dir)) {
            std::filesystem::create_directories(clean_data_dir);
            spdlog::info("Created directory for cleaned data: {}", clean_data_dir.string());
        }
        else {
            spdlog::info("Cleaned data directory already exists: {}", clean_data_dir.string());
        }

        ThreadPool pool(thread_num);
        for (auto it = species_path_map.begin(); it != species_path_map.end(); ++it) {
            std::string species = it->first;
            std::filesystem::path raw_path = it->second;

            pool.enqueue([species, raw_path, &workdir_path, &species_path_map]() {
                FilePath out_dir = workdir_path / DATA_DIR / CLEAN_DATA_DIR;
                FilePath out_fasta = out_dir / (species + ".fasta");

                // 检查清洗后的文件是否已存在
                if (std::filesystem::exists(out_fasta)) {
                    spdlog::warn("Species {} cleaned file already exists, skipping: {}",
                        species, out_fasta.string());
                    species_path_map[species] = out_fasta;
                    return;
                }

                SeqPro::utils::cleanFastaFile(raw_path, out_fasta, 60);
                species_path_map[species] = out_fasta;
                spdlog::info("Species {} cleaned file path updated to: {}", species,
                    out_fasta.string());
                });
        }
        pool.waitAllTasksDone();
        spdlog::info("All raw sequence data cleaned successfully.");
        return true;
    }
    catch (const std::exception& error) {
        spdlog::error("Error occurred during data cleaning: {}", error.what());
        throw std::runtime_error("Failed to clean raw dataset.");
    }
}

// -----------------------------
// 运行 WindowMasker 对 FASTA 文件进行处理
// -----------------------------
std::map<SpeciesName, FilePath>
repeatSeqMasking(const FilePath workdir_path,
    const SpeciesPathMap& species_path_map_const, int thread_num) {
    SpeciesPathMap species_path_map =
        species_path_map_const; // Create a mutable copy
    std::map<SpeciesName, FilePath> interval_files_map;
    const std::string DATA_DIR_NAME = "data";
    const std::string MASKED_DATA_DIR_NAME = "masked_data";

    try {
        FilePath masked_data_dir =
            workdir_path / DATA_DIR_NAME / MASKED_DATA_DIR_NAME;
        if (!std::filesystem::exists(masked_data_dir)) {
            std::filesystem::create_directories(masked_data_dir);
            spdlog::info("Created directory for windowmasker output: {}",
                masked_data_dir.string());
        }
        else {
            spdlog::info("Windowmasker output directory already exists: {}",
                masked_data_dir.string());
        }

        // 获取 windowmasker 路径
        std::string windowmasker_path;
        try {
            windowmasker_path = findWindowmaskerPath();
        }
        catch (const std::exception& e) {
            spdlog::error("Failed to locate windowmasker: {}", e.what());
            throw std::runtime_error("Cannot proceed with repeat masking: " + std::string(e.what()));
        }

        ThreadPool pool(thread_num);

        for (auto it = species_path_map.begin(); it != species_path_map.end();
            ++it) {
            std::string species_key = it->first;
            FilePath input_fasta_path = it->second; // 待处理的FASTA文件

            pool.enqueue([species_key, input_fasta_path, masked_data_dir,
                &species_path_map, &interval_files_map, windowmasker_path]() {
                    FilePath counts_file = masked_data_dir / (species_key + ".counts");
                    FilePath interval_file = masked_data_dir / (species_key + ".interval");

                    // 检查最终的interval文件是否已存在
                    if (std::filesystem::exists(interval_file)) {
                        spdlog::warn("[{}] WindowMasker output interval file already exists, "
                            "skipping: {}",
                            species_key, interval_file.string());
                        interval_files_map[species_key] =
                            interval_file; // Store in the new map
                        return;
                    }

                    // 命令1: windowmasker -mk_counts
                    std::string cmd1_str = "\"" + windowmasker_path + "\" -mk_counts -in \"" +
                        input_fasta_path.string() + "\" -out \"" +
                        counts_file.string() + "\" -sformat obinary";
                    spdlog::info("[{}] Executing: {}", species_key, cmd1_str);
                    int ret1 = std::system(cmd1_str.c_str());

                    if (ret1 != 0) {
                        spdlog::error("[{}] windowmasker -mk_counts failed with exit code "
                            "{}. Command: {}",
                            species_key, ret1, cmd1_str);
                        if (std::filesystem::exists(counts_file)) {
                            std::filesystem::remove(counts_file);
                        }
                        return; // 该物种任务失败
                    }
                    spdlog::info("[{}] windowmasker -mk_counts successful.", species_key);

                    // 命令2: windowmasker -ustat
                    std::string cmd2_str =
                        "\"" + windowmasker_path + "\" -ustat \"" + counts_file.string() + "\" -in \"" +
                        input_fasta_path.string() + "\" -out \"" + interval_file.string() +
                        "\" -outfmt interval -dust true";
                    spdlog::info("[{}] Executing: {}", species_key, cmd2_str);
                    int ret2 = std::system(cmd2_str.c_str());

                    if (ret2 != 0) {
                        spdlog::error(
                            "[{}] windowmasker -ustat failed with exit code {}. Command: {}",
                            species_key, ret2, cmd2_str);
                        if (std::filesystem::exists(interval_file)) {
                            std::filesystem::remove(interval_file);
                        }
                        return; // 该物种任务失败
                    }
                    spdlog::info("[{}] windowmasker -ustat successful. Output: {}",
                        species_key, interval_file.string());

                    interval_files_map[species_key] = interval_file;
                    spdlog::info("[{}] Interval file generated at: {}", species_key,
                        interval_file.string());
                });
        }

        pool.waitAllTasksDone();

        spdlog::info("All windowmasker tasks submitted and processing attempted. "
            "Check logs for individual species status.");
        return interval_files_map; // 返回包含interval文件路径的映射

    }
    catch (const std::exception& error) {
        spdlog::error("Critical error during windowmasker processing setup or "
            "thread pool execution: {}",
            error.what());
        throw std::runtime_error(
            std::string(
                "Failed to run windowmasker processing due to critical error: ") +
            error.what());
    }
    return {}; // 发生严重错误时返回空映射
}

// 辅助函数：解析区间字符串，格式如 "start - end"
std::pair<size_t, size_t> parseIntervalLine(const std::string& line) {
    std::stringstream ss(line);
    size_t start, end;
    char dash;
    ss >> start >> dash >> end;
    if (ss.fail() || dash != '-') {
        throw std::runtime_error("Invalid interval format: " + line);
    }
    return { start, end };
}
//
// // 内部辅助函数：根据区间文件对单个FASTA文件进行掩码处理
// FilePath applyMaskToFastaInternal(const FilePath& original_fasta_path,
//     const FilePath& interval_file_path,
//     const FilePath& output_masked_fasta_path,
//     const std::string& species_key) {
//
//     spdlog::info("[{}] Applying mask using {} to {} -> {}", species_key,
//         interval_file_path.string(), original_fasta_path.string(),
//         output_masked_fasta_path.string());
//
//     // 1. 解析区间文件
//     std::map<std::string, std::vector<std::pair<size_t, size_t>>>
//         intervals_by_seq;
//     std::ifstream interval_stream(interval_file_path);
//     if (!interval_stream.is_open()) {
//         throw std::runtime_error("Failed to open interval file: " +
//             interval_file_path.string());
//     }
//     std::string line;
//     std::string current_seq_name;
//     while (std::getline(interval_stream, line)) {
//         if (line.empty())
//             continue;
//         if (line[0] == '>') {
//             current_seq_name = line.substr(1);
//             current_seq_name.erase(0,
//                 current_seq_name.find_first_not_of(" \t\n\r\f\v"));
//             current_seq_name.erase(current_seq_name.find_last_not_of(" \t\n\r\f\v") +
//                 1);
//         }
//         else if (!current_seq_name.empty()) {
//             try {
//                 intervals_by_seq[current_seq_name].push_back(parseIntervalLine(line));
//             }
//             catch (const std::runtime_error& e) {
//                 spdlog::warn("[{}] Skipping invalid interval line in {}: {} - {}",
//                     species_key, interval_file_path.string(), line, e.what());
//             }
//         }
//     }
//     interval_stream.close();
//
//     if (intervals_by_seq.empty()) {
//         spdlog::warn("[{}] No intervals found or all intervals were invalid in {}. "
//             "Copying original file.",
//             species_key, interval_file_path.string());
//         copyLocalFile(original_fasta_path, output_masked_fasta_path);
//         return output_masked_fasta_path;
//     }
//
//     // 2. 读取原始FASTA，应用掩码，写入新FASTA
//     FastaManager reader(original_fasta_path);
//     std::ofstream out_stream(output_masked_fasta_path);
//     if (!out_stream.is_open()) {
//         throw std::runtime_error("Failed to open output masked FASTA file: " +
//             output_masked_fasta_path.string());
//     }
//
//     std::string header, sequence;
//     const size_t line_width = 60; // 标准FASTA行宽
//
//     while (reader.nextRecord(header, sequence)) {
//         std::string clean_header = header;
//         size_t first_space = header.find_first_of(" \t");
//         if (first_space != std::string::npos) {
//             clean_header = header.substr(0, first_space);
//         }
//
//         if (intervals_by_seq.count(clean_header)) {
//             auto& intervals = intervals_by_seq[clean_header];
//             std::sort(intervals.begin(), intervals.end());
//
//             for (const auto& interval : intervals) {
//                 size_t start_1based = interval.first;
//                 size_t end_1based = interval.second;
//
//                 if (start_1based == 0) {
//                     spdlog::warn("[{}] Invalid 0-based start interval {} for sequence {} "
//                         "in {}. Skipping.",
//                         species_key, start_1based, clean_header,
//                         interval_file_path.string());
//                     continue;
//                 }
//                 size_t start_0based = start_1based - 1;
//                 size_t end_0based = end_1based - 1;
//
//                 if (start_0based < sequence.length() &&
//                     end_0based < sequence.length() && start_0based <= end_0based) {
//                     for (size_t i = start_0based; i <= end_0based; ++i) {
//                         sequence[i] = std::tolower(sequence[i]);
//                     }
//                 }
//                 else {
//                     spdlog::warn("[{}] Interval {}-{} out of bounds for sequence {} (len "
//                         "{}) in {}. Skipping.",
//                         species_key, start_1based, end_1based, clean_header,
//                         sequence.length(), interval_file_path.string());
//                 }
//             }
//         }
//
//         out_stream << ">" << header << "\n";
//         for (size_t i = 0; i < sequence.length(); i += line_width) {
//             out_stream << sequence.substr(i,
//                 std::min(line_width, sequence.length() - i))
//                 << "\n";
//         }
//     }
//     out_stream.close();
//     spdlog::info("[{}] Successfully applied mask and wrote to {}", species_key,
//         output_masked_fasta_path.string());
//     return output_masked_fasta_path;
// }
//
// // 根据区间文件应用掩码并更新物种路径映射
// bool applyMaskingAndUpdatePaths(
//     const FilePath workdir_path,
//     SpeciesPathMap& species_path_map, // 将被更新为指向掩码后的文件
//     const std::map<SpeciesName, FilePath>& interval_files_map, int thread_num) {
//
//     const std::string MASKED_FINAL_DIR_NAME = "masked_final_sequences";
//     FilePath masked_final_dir = workdir_path / DATA_DIR / MASKED_FINAL_DIR_NAME;
//
//     try {
//         if (!std::filesystem::exists(masked_final_dir)) {
//             std::filesystem::create_directories(masked_final_dir);
//             spdlog::info("Created directory for final masked FASTA files: {}",
//                 masked_final_dir.string());
//         }
//         else {
//             spdlog::info("Final masked FASTA output directory already exists: {}",
//                 masked_final_dir.string());
//         }
//
//         ThreadPool pool(thread_num);
//         std::vector<std::future<std::pair<SpeciesName, FilePath>>> task_futures;
//
//         for (auto const& [species_key, interval_file_path] : interval_files_map) {
//             if (!species_path_map.count(species_key)) {
//                 spdlog::error("[{}] Original FASTA path (cleaned) not found in "
//                     "species_path_map. Skipping masking for this species.",
//                     species_key);
//                 continue;
//             }
//             const FilePath& original_cleaned_fasta_path =
//                 species_path_map.at(species_key);
//
//             FilePath output_masked_fasta_path =
//                 masked_final_dir / (species_key + ".masked.fasta");
//
//             task_futures.emplace_back(
//                 pool.enqueue([original_cleaned_fasta_path, interval_file_path,
//                     output_masked_fasta_path,
//                     species_key]() -> std::pair<SpeciesName, FilePath> {
//                         try {
//                             FilePath final_masked_path = applyMaskToFastaInternal(
//                                 original_cleaned_fasta_path, interval_file_path,
//                                 output_masked_fasta_path, species_key);
//                             return { species_key, final_masked_path };
//                         }
//                         catch (const std::exception& e) {
//                             spdlog::error(
//                                 "[{}] Error applying mask: {}. Original path will be kept.",
//                                 species_key, e.what());
//                             return { species_key, FilePath() }; // 失败时返回空路径
//                         }
//                     }));
//         }
//
//         for (auto& fut : task_futures) {
//             std::pair<SpeciesName, FilePath> result = fut.get();
//             if (!result.second.empty()) {
//                 species_path_map[result.first] = result.second;
//                 spdlog::info("[{}] Path updated to masked file: {}", result.first,
//                     result.second.string());
//             }
//             else {
//                 spdlog::warn("[{}] Masking failed or skipped. Path not updated in "
//                     "species_path_map.",
//                     result.first);
//             }
//         }
//
//         spdlog::info(
//             "All masking tasks for apply_masking_and_update_paths completed.");
//         return true;
//
//     }
//     catch (const std::exception& error) {
//         spdlog::error("Critical error during masked FASTA generation setup: {}",
//             error.what());
//         return false;
//     }
// }

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
            spdlog::error("curl_easy_perform failed for URL {}: {}", pathStr,
                curl_easy_strerror(res));
            curl_easy_cleanup(curl);
            return "0 B";
        }

        curl_off_t content_length = 0;
        res = curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T,
            &content_length);
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
            spdlog::error("Failed to get file size for {}: {}", filePath.string(),
                e.what());
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
        if (curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T,
            &content_length) != CURLE_OK ||
            content_length < 0) {
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
        temp_file = input_path.parent_path() /
            (input_path.filename().string() + "_in_process");
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

// -----------------------------
// parseSeqfile：解析 seqfile，同时使用 NewickParser 对第一行 Newick 树进行解析
// -----------------------------
bool parseSeqfile(const FilePath& seqfile_path,
    NewickParser& newick_tree,
    SpeciesPathMap& species_map,
    const std::string& root /* = "" */) {
    if (!std::filesystem::exists(seqfile_path)) {
        throw std::runtime_error("parseSeqfile: cannot locate file: " +
            seqfile_path.string());
    }

    std::ifstream ifs(seqfile_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("parseSeqfile: failed to open: " +
            seqfile_path.string());
    }

    std::string line;
    bool got_tree = false;

    species_map.clear();
    newick_tree.clear();

    // 用于子树过滤的叶子名集合
    std::unordered_set<std::string> allowed_leaves;

    while (std::getline(ifs, line)) {
        auto l = line.find_first_not_of(" \t\r\n");
        if (l == std::string::npos) continue;
        auto r = line.find_last_not_of(" \t\r\n");
        std::string trimmed = line.substr(l, r - l + 1);

        if (!got_tree) {
            // 第一条非空行为 Newick
            try {
                newick_tree.parse(trimmed);
            }
            catch (const std::runtime_error& e) {
                throw std::runtime_error(
                    std::string("parseSeqfile: invalid Newick tree: ") + e.what());
            }
            got_tree = true;

            // 如果指定了 root，则裁剪到该子树
            if (!root.empty()) {
                int rootId = newick_tree.findNodeIdByName(root);
                if (rootId == -1) {
                    throw std::runtime_error("parseSeqfile: root name `" + root + "` not found in Newick tree");
                }
                newick_tree.restrictToSubtreeByRootId(rootId);
            }

            // 构建允许的叶子名集合（用于过滤 species）
            {
                auto leafNames = newick_tree.getLeafNames();
                allowed_leaves.clear();
                for (auto& nm : leafNames) {
                    std::string tmp = nm;
                    // 与 trimString 语义一致，防御性去空白
                    newick_tree.trimString(tmp);
                    if (!tmp.empty()) allowed_leaves.insert(std::move(tmp));
                }
            }

            continue;
        }

        // 其余行：<speciesName> <path>
        std::istringstream iss(trimmed);
        std::string speciesName;
        std::string filePath;
        if (!(iss >> speciesName)) {
            continue;
        }
        if (!(iss >> std::ws) || !std::getline(iss, filePath)) {
            throw std::runtime_error("parseSeqfile: missing file path for species `" +
                speciesName + "`");
        }
        auto p_l = filePath.find_first_not_of(" \t\r\n");
        auto p_r = filePath.find_last_not_of(" \t\r\n");
        if (p_l == std::string::npos) {
            throw std::runtime_error("parseSeqfile: empty file path for species `" +
                speciesName + "`");
        }
        filePath = filePath.substr(p_l, p_r - p_l + 1);

        // 如果指定了 root，只接收子树中的叶子；否则全部接收
        if (!root.empty()) {
            // speciesName 必须在 allowed_leaves 里（和叶节点名匹配）
            if (allowed_leaves.find(speciesName) == allowed_leaves.end()) {
                // 不属于该子树，跳过
                continue;
            }
        }

        species_map.emplace(speciesName, FilePath(filePath));
    }

    ifs.close();

    if (!got_tree) {
        throw std::runtime_error(
            "parseSeqfile: no non-empty line found as Newick tree");
    }
    if (species_map.empty()) {
        // 如果有 root，空也可能意味着 root 子树没有任何叶子/或映射里没给这些叶子
        throw std::runtime_error("parseSeqfile: no species=>path mappings found"
            + std::string(root.empty() ? "" : " under root `" + root + "`"));
    }

    return true;
}

//
// void repeatMaskRawData(const FilePath& work_dir, int thread_num,
//     SpeciesPathMap& species_path_map) {
//     spdlog::info("Repeat masking enabled. Generating interval files...");
//
//     // 1. 运行 WindowMasker 生成 interval 文件
//     std::map<SpeciesName, FilePath> interval_files_map =
//         repeatSeqMasking(work_dir, species_path_map, thread_num);
//
//     if (interval_files_map.empty()) {
//         spdlog::warn("No interval files were generated by repeatSeqMasking. "
//             "Skipping application of masks.");
//         return;
//     }
//
//     spdlog::info(
//         "Applying masks to FASTA files based on generated interval files...");
//     // 2. 应用掩码并更新 species_path_map 指向新的掩码后文件
//     bool masking_applied_successfully = applyMaskingAndUpdatePaths(
//         work_dir, species_path_map, interval_files_map, thread_num);
//
//     if (!masking_applied_successfully) {
//         spdlog::error("Failed to apply masks or update species paths. Subsequent "
//             "steps will use unmasked data.");
//         return;
//     }
//
//     spdlog::info("Successfully applied masks and updated species paths.");
//
//     // 3. 隐藏遮蔽区间并生成新的 FASTA 和 FAI 文件
//     spdlog::info("Hiding masked intervals and generating new FASTA files with "
//         "FAI indices...");
//     FilePath hidden_data_dir = work_dir / "data" / "hidden_intervals";
//     if (!std::filesystem::exists(hidden_data_dir)) {
//         std::filesystem::create_directories(hidden_data_dir);
//         spdlog::info("Created directory for hidden intervals data: {}",
//             hidden_data_dir.string());
//     }
//
//     for (const auto& [species_name, interval_file_path] : interval_files_map) {
//         // 仅处理在 species_path_map 中存在的物种
//         if (!species_path_map.count(species_name)) {
//             continue;
//         }
//
//         FilePath input_fasta = species_path_map[species_name];
//         FilePath output_fasta = hidden_data_dir / (species_name + ".hidden.fasta");
//         FilePath output_fai =
//             hidden_data_dir / (species_name + ".hidden.fasta.fai");
//
//         try {
//             FastaManager fasta_manager(input_fasta);
//             bool hide_success = fasta_manager.hideIntervalsFromFileAndGenerateFai(
//                 interval_file_path, output_fasta, output_fai,
//                 60 // 每行最大字符数
//             );
//
//             if (hide_success) {
//                 spdlog::info("[{}] Successfully generated hidden intervals FASTA: {}",
//                     species_name, output_fasta.string());
//                 // 更新映射，使后续步骤使用屏蔽后序列
//                 species_path_map[species_name] = output_fasta;
//             }
//             else {
//                 spdlog::error("[{}] Failed to generate hidden intervals FASTA",
//                     species_name);
//             }
//         }
//         catch (const std::exception& e) {
//             spdlog::error("[{}] Error processing hidden intervals: {}", species_name,
//                 e.what());
//         }
//     }
// }