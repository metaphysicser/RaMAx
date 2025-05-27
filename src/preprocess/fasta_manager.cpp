#include "data_process.h"

// -----------------------------
// FastaManager: 管理 FASTA 文件的读取、清洗、索引等操作
// -----------------------------
//FastaManager::FastaManager()
//{
//
//}

//FastaManager::FastaManager(const FilePath& fasta_path, const FilePath& fai_path)
//    : fasta_path_(fasta_path)
//    , fai_path_(fai_path)
//{
//    // If fai_path is not empty, load the FAI records into an unordered_map
//    if (!fai_path.empty()) {
//        loadFaiRecords(fai_path);
//    }
//    open(); // Open FASTA file
//}

//FastaManager::~FastaManager() {
//    if (seq_) {
//        kseq_destroy(seq_);
//        seq_ = nullptr;
//    }
//    if (fp_) {
//        gzclose(fp_);
//        fp_ = nullptr;
//    }
//}
// 加载 .fai 索引文件，填充到 fai_records 中
void FastaManager::loadFaiRecords(const FilePath& fai_path)
{
    std::ifstream in(fai_path);
    if (!in) {
        spdlog::error("Failed to open fai file: {}", fai_path.string());
        throw std::runtime_error("Cannot open " + fai_path.string());
    }

    std::string line;

    // 读取第一行：是否包含 'N' 字符（YES 或 NO）
    std::getline(in, line);
    if (line == "YES") {
        has_n_in_fasta = true;
    }
    else if (line == "NO") {
        has_n_in_fasta = false;
    }
    else {
        spdlog::error("Invalid FAI file format: {}", fai_path.string());
        throw std::runtime_error("Invalid FAI file format: " + fai_path.string());
    }

    fai_records.clear();

    // 每一行一个序列记录
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);

        FaiRecord rec;
        iss >> rec.seq_name >> rec.global_start_pos >> rec.length >> rec.offset >> rec.line_bases >> rec.line_bytes;

        if (iss.fail()) {
            spdlog::warn("Skipping malformed line in {}: {}", fai_path.string(), line);
            continue;
        }
        fai_records.push_back(rec);
    }
    in.close();
}

// 读取下一个序列记录（通过 kseq），返回 false 表示读完
bool FastaProcessor::nextRecord(std::string& header, std::string& sequence) {
    int ret = kseq_read(seq_.get());
    if (ret < 0) return false;

    header.assign(seq_->name.s);
    sequence.assign(seq_->seq.s, seq_->seq.l);

    if (clean_data_) {
        cleanSequence(sequence);
    }
    return true;
}

// 关闭并重新打开 fasta 文件流，重置读取状态
void FastaProcessor::reset() {
    gz_file_wrapper_.reset();
    seq_.reset();
    fasta_open();
}

// 读取若干序列片段并拼接成一个长序列，设置终止符 terminator
std::string FastaProcessor::concatRecords(char terminator, size_t limit) {
    reset();
    std::ostringstream oss;
    std::string hdr, seq;
    size_t count = 0;

    while (nextRecord(hdr, seq) && count < limit) {
        oss << seq;
        ++count;
    }
    std::string result = oss.str();
    if (!result.empty()) {
        result.back() = terminator;
    }
    return result;
}

// 统计所有序列的数量、长度、最小值、最大值、平均值
FastaProcessor::Stats FastaProcessor::getStats() {
    reset();
    Stats s;
    std::string hdr, seq;
    while (nextRecord(hdr, seq)) {
        size_t len = seq.size();
        s.record_count++;
        s.total_bases += len;
        s.min_len = std::min(s.min_len, len);
        s.max_len = std::max(s.max_len, len);
    }
    if (s.record_count > 0) {
        s.average_len = s.total_bases / s.record_count;
    }
    else {
        s.min_len = 0;
    }
    return s;
}

// 将所有碱基字符转为大写，非法字符统一设为 N
void FastaProcessor::cleanSequence(std::string& seq) {
    for (char& c : seq) {
        unsigned char uc = static_cast<unsigned char>(c);
        uc = std::toupper(uc);
        if (!has_n_in_fasta && uc == 'N') {
            has_n_in_fasta = true;
        }
        if (uc != 'A' && uc != 'C' && uc != 'G' && uc != 'T' && uc != 'N') {
            uc = 'N';
        }
        c = static_cast<char>(uc);
    }
}

// 打开 gzip 格式的 fasta 文件，创建 kseq 结构体
void FastaProcessor::fasta_open() {
    gz_file_wrapper_ = std::make_unique<GzFileWrapper>(fasta_path_.string());
    seq_.reset(kseq_init(gz_file_wrapper_->get()));
}

// 如果 FAI 文件不存在，重新扫描 fasta 生成 FAI 文件（写入 YES/NO + 多行记录）
bool FastaManager::reScanAndWriteFai(const FilePath& fa_path,
    const FilePath& fai_path,
    size_t line_width) const {
    if (std::filesystem::exists(fai_path)) {
        spdlog::warn("Fai file already exists: {}", fai_path.string());
        return true;
    }

    std::ifstream in(fa_path);
    if (!in) {
        spdlog::error("Failed to open {} for reading", fa_path.string());
        return false;
    }

    FilePath tmp_fai_path = getTempFilePath(fai_path);
    std::ofstream out(tmp_fai_path);
    if (!out) {
        spdlog::error("Failed to open {} for writing", tmp_fai_path.string());
        return false;
    }

    out << (has_n_in_fasta ? "YES" : "NO") << "\n";

    uint_t global_start_pos = 0;
    size_t global_offset = 0;
    std::string line;
    std::string seq_name;
    size_t seq_len = 0;
    size_t seq_start = 0;
    bool reading_seq = false;
    size_t line_bytes_in_fasta = line_width + 1; // 包括换行符

    while (std::getline(in, line)) {
        size_t this_line_bytes = line.size() + 1;

        if (!line.empty() && line[0] == '>') {
            if (reading_seq) {
                out << seq_name << "\t" << global_start_pos << "\t" << seq_len << "\t"
                    << seq_start << "\t" << line_width << "\t" << line_bytes_in_fasta << "\n";
                global_start_pos += seq_len;
            }
            seq_name = line.substr(1);
            seq_len = 0;
            reading_seq = true;
            seq_start = global_offset + this_line_bytes;
        }
        else {
            seq_len += line.size();
        }
        global_offset += this_line_bytes;
    }

    if (reading_seq) {
        out << seq_name << "\t" << global_start_pos << "\t" << seq_len << "\t"
            << seq_start << "\t" << line_width << "\t" << line_bytes_in_fasta << "\n";
    }

    in.close();
    out.close();
    std::filesystem::rename(tmp_fai_path, fai_path);
    spdlog::info("FAI index created: {}", fai_path.string());
    return true;
}


FilePath FastaProcessor::writeCleanedFasta(const FilePath& output_file, uint64_t line_width) {
    if (std::filesystem::exists(output_file)) {
        spdlog::warn("Output FASTA already exists: {}", output_file.string());
        return output_file;
    }

    FilePath tmp_path = getTempFilePath(output_file);
    std::ofstream ofs_fasta(tmp_path);
    if (!ofs_fasta) {
        throw std::runtime_error("Failed to create " + tmp_path.string());
    }

    reset();

    std::string header, sequence;
    while (nextRecord(header, sequence)) {
        ofs_fasta << ">" << header << "\n";
        cleanSequence(sequence);
        size_t seq_len = sequence.size();
        size_t start = 0;
        while (start < seq_len) {
            size_t chunk_end = std::min(start + line_width, (uint64_t)seq_len);
            ofs_fasta.write(&sequence[start], chunk_end - start);
            ofs_fasta.put('\n');
            start = chunk_end;
        }
    }

    ofs_fasta.close();
    std::filesystem::rename(tmp_path, output_file);
    spdlog::info("Wrote cleaned FASTA: {}", output_file.string());

    return output_file;
}

FilePath FastaManager::cleanAndIndexFasta(const FilePath& output_dir,
    const std::string& prefix,
    uint64_t line_width)
{
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
        spdlog::info("Created directory: {}", output_dir.string());
    }

    FilePath out_fasta = output_dir / (prefix + ".fasta");
    FilePath out_fai = output_dir / (prefix + ".fasta.fai");

    writeCleanedFasta(out_fasta, line_width);
    reScanAndWriteFai(out_fasta, out_fai, line_width);

    this->fai_path_ = out_fai;
    return out_fasta;
}

std::string FastaManager::getSubConcatSequence(size_t start, size_t length)
{
    // 如果请求长度为 0，直接返回空
    if (length == 0) {
        return "";
    }

    // 将结果存储在 result 中
    std::string result;
    result.reserve(length); // 预留长度, 避免反复分配

    // current_offset 表示当前遍历到的序列在“全局坐标”中的起始位置
    size_t current_offset = 0;

    // 剩余需要读取的长度
    size_t remain = length;

    // 遍历 FAI 中所有记录（假设顺序与原始文件一致）
    for (const auto& rec : fai_records) {
        // 每条序列在全局坐标中的范围是 [current_offset, current_offset + rec.length - 1]
        size_t seq_start_global = current_offset;
        size_t seq_end_global = current_offset + rec.length - 1;

        // 判断请求的 start 是否落在这条序列上（或者部分区间与这条序列重叠）
        if (start <= seq_end_global && (start + remain - 1) >= seq_start_global) {
            // 计算在本序列中的局部起始位置 local_start
            // 例如，如果 start=1500 而这条序列 global 区间 [1000..2999]，则 local_start=1500-1000=500
            size_t local_start = (start > seq_start_global) ? (start - seq_start_global) : 0;

            // 当前序列从 local_start 到末端还剩多少碱基
            size_t can_read_in_this_seq = rec.length - local_start;

            // 我们实际要从本序列读出的长度
            size_t to_read = (remain < can_read_in_this_seq) ? remain : can_read_in_this_seq;

            // 调用已有的函数，从当前序列中读出区间
            // 注意该函数是基于 “序列本地坐标”(start_in_seq, length_in_seq)
            std::string part = getSubSequence(rec.seq_name, local_start, to_read);
            result += part;

            // 更新剩余需要读取的长度
            remain -= to_read;
            // 更新全局坐标的起点（因为我们已经读完这部分）
            start += to_read;

            // 如果剩余长度为 0，说明已完成读取
            if (remain == 0) {
                break;
            }
        }

        // 更新下一条序列的 global offset
        current_offset += rec.length;
    }

    // 如果循环结束后 remain 还不为 0，说明请求区间超出了所有序列总长度
    // 你可以选择抛出异常或仅返回已能读取的部分
    if (remain > 0) {
        throw std::runtime_error("Requested range exceeds total length of all sequences.");
    }

    return result;
}

std::string FastaManager::getSubSequence(const std::string& seq_name, size_t start, size_t length)
{
    // 使用 std::vector 查找目标序列
    FaiRecord rec;
    bool find = false;
    for (const auto& fai_record : fai_records) {
        if (fai_record.seq_name == seq_name) {
            rec = fai_record;
            find = true;
            break;
        }
    }

    if (!find) {
        throw std::runtime_error("Cannot find sequence " + seq_name + " in FAI records.");
    }

    if ( start+length > rec.length) {
        throw std::runtime_error("Invalid range [" + std::to_string(start) + ", " + std::to_string(start+length) +
            "] for sequence " + seq_name);
    }

    // 使用 mmap 映射 FASTA 文件到内存
    int fd = open(fasta_path_.c_str(), O_RDONLY);
    if (fd == -1) {
        throw std::runtime_error("Failed to open " + fasta_path_.string() + " for reading sub-sequence.");
    }

    // 获取文件大小
    off_t file_size = lseek(fd, 0, SEEK_END);
    if (file_size == -1) {
        close(fd);
        throw std::runtime_error("Failed to determine the size of " + fasta_path_.string());
    }

    // 将文件映射到内存中
    char* file_data = (char*)mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
    if (file_data == MAP_FAILED) {
        close(fd);
        throw std::runtime_error("Failed to mmap file " + fasta_path_.string());
    }

    // 关闭文件描述符，mmap 后不再需要它
    close(fd);

    // 计算请求的子序列的起始位置
    size_t line_bases = rec.line_bases;
    size_t line_bytes = rec.line_bytes;
    size_t seq_offset = rec.offset;

    size_t rowIndex = start / line_bases;
    size_t colIndex = start % line_bases;
    size_t filePos = seq_offset + rowIndex * line_bytes + colIndex;

    // 计算需要读取的字符数
    size_t req_len = length;
    std::string result;
    result.reserve(req_len);

    // 使用内存映射的数据进行读取
    size_t to_read = req_len;
    size_t current_pos = filePos;

    while (to_read > 0) {
        char c = file_data[current_pos];
        if (c != '\n' && c != '\r') {  // 忽略换行符
            result.push_back(c);
            --to_read;
        }
        ++current_pos;
        if (current_pos >= file_size) {
            break; // 如果超出文件范围
        }
    }

    // 释放内存映射
    munmap(file_data, file_size);

    if (to_read > 0) {
        throw std::runtime_error("Reached EOF unexpectedly while reading sub-sequence for " + seq_name);
    }

    return result;
}

uint_t FastaManager::getConcatSeqLength() {
    uint_t total_length = 0;
	for (const auto& rec : fai_records) {
		total_length += rec.length;
	}
    return total_length;
}

ChrName FastaManager::getChrName(uint_t global_start, uint_t length) {
    // fai_records 必须按 global_start_pos 升序排列
    size_t l = 0, r = fai_records.size();
    // 二分查找：找到第一个 global_start_pos > global_start
    while (l < r) {
        size_t m = l + (r - l) / 2;
        if (fai_records[m].global_start_pos > global_start)
            r = m;
        else
            l = m + 1;
    }
    if (l == 0) {
        throw std::out_of_range("请求的 global_start 在第一个染色体之前");
    }
    const FaiRecord& rec = fai_records[l - 1];

    // 检查这段区间是否越界到下一个染色体
    uint_t chr_end_global = rec.global_start_pos + rec.length;
    if (global_start + length > chr_end_global) {
        return "";
    }

    // 计算在该染色体上的局部坐标
    return rec.seq_name;
}

// data_process.cpp

RegionVec FastaManager::preAllocateChunks(uint_t chunk_size, uint_t overlap_size)
{
    RegionVec chunks;
    // 预估总 chunk 数以减少 realloc（可选）
    size_t total_len = getConcatSeqLength();
    size_t est_chunks = (total_len + chunk_size - 1) / chunk_size;
    chunks.reserve(est_chunks);

    // 对每条染色体分别切分
    for (const auto& rec : fai_records) {
        const auto& chr = rec.seq_name;
        size_t chr_len = rec.length;

        // 如果染色体长度小于等于 chunk_size，则只生成一个不重叠 chunk
        if (chr_len <= chunk_size) {
            chunks.push_back({ chr, 0, chr_len });
            continue;
        }

        // 多个 chunk，需要在它们之间保留 overlap_size
        size_t start = 0;
        while (start < chr_len) {
            // 本 chunk 的实际长度（最后一块可能不足 chunk_size）
            size_t this_len = std::min<size_t>(chunk_size,
                chr_len - start);
            chunks.push_back({ chr, start, this_len });

            // 计算下一个 chunk 的起始位置：前进 chunk_size，然后回退 overlap_size
            if (start + this_len >= chr_len) {
                break;
            }
            start += chunk_size;
            // 保证不回退超过已经走过的距离
            start = (start >= overlap_size ? start - overlap_size : 0);
        }
    }
    return chunks;
}


// 通用的隐藏区间函数，根据提供的区间数据隐藏指定区间，生成新的FASTA文件和对应的FAI索引
// intervals_by_seq: 按序列名称组织的区间映射，格式为 {seq_name: [(start1, end1), (start2, end2), ...]}
// output_fasta_path: 输出的隐藏区间后的FASTA文件路径
// output_fai_path: 输出的FAI索引文件路径
// line_width: FASTA文件的行宽
// 返回值: 成功返回true，失败返回false
bool FastaManager::hideIntervalsAndGenerateFai(
    const std::map<std::string, std::vector<std::pair<size_t, size_t>>>& intervals_by_seq,
    const FilePath& output_fasta_path,
    const FilePath& output_fai_path,
    size_t line_width) {
    
    try {
        // 1. 复制区间数据并进行处理
        std::map<std::string, std::vector<std::pair<size_t, size_t>>> processed_intervals = intervals_by_seq;
        
        // 2. 对每个序列的区间进行排序和合并
        for (auto& [seq_name, intervals] : processed_intervals) {
            if (intervals.empty()) continue;
            
            // 按起始位置排序
            std::sort(intervals.begin(), intervals.end());
            
            // 合并重叠的区间
            std::vector<std::pair<size_t, size_t>> merged_intervals;
            merged_intervals.push_back(intervals[0]);
            
            for (size_t i = 1; i < intervals.size(); ++i) {
                auto& last = merged_intervals.back();
                auto& current = intervals[i];
                
                if (current.first <= last.second + 1) {
                    // 区间重叠或相邻，合并
                    last.second = std::max(last.second, current.second);
                } else {
                    // 不重叠，添加新区间
                    merged_intervals.push_back(current);
                }
            }
            
            processed_intervals[seq_name] = std::move(merged_intervals);
        }
        
        // 3. 创建输出文件
        std::ofstream out_fasta(output_fasta_path);
        if (!out_fasta.is_open()) {
            spdlog::error("Failed to create output FASTA file: {}", output_fasta_path.string());
            return false;
        }
        
        std::ofstream out_fai(output_fai_path);
        if (!out_fai.is_open()) {
            spdlog::error("Failed to create output FAI file: {}", output_fai_path.string());
            return false;
        }
        
        // 写入FAI文件头部信息（是否包含N字符）
        out_fai << (has_n_in_fasta ? "YES" : "NO") << "\n";
        
        // 4. 处理每个序列
        reset(); // 重置FASTA读取器
        std::string header, sequence;
        uint_t global_start_pos = 0;
        size_t global_offset = 0;
        
        while (nextRecord(header, sequence)) {
            // 获取清理后的序列名称
            std::string clean_header = header;
            size_t first_space = header.find_first_of(" \t");
            if (first_space != std::string::npos) {
                clean_header = header.substr(0, first_space);
            }
            
            // 获取该序列需要隐藏的区间
            std::vector<std::pair<size_t, size_t>> hide_intervals;
            if (processed_intervals.count(clean_header)) {
                hide_intervals = processed_intervals[clean_header];
            }
            
            // 5. 生成保留的片段
            std::vector<std::string> segments;
            std::vector<std::pair<size_t, size_t>> segment_coords; // 记录每个片段在原序列中的坐标
            
            if (hide_intervals.empty()) {
                // 没有需要隐藏的区间，保留整个序列
                segments.push_back(sequence);
                segment_coords.push_back({1, sequence.length()}); // 1-based坐标
            } else {
                size_t current_pos = 1; // 1-based坐标
                
                for (const auto& interval : hide_intervals) {
                    size_t hide_start = interval.first;
                    size_t hide_end = interval.second;
                    
                    // 添加隐藏区间之前的片段
                    if (current_pos < hide_start) {
                        size_t segment_start = current_pos - 1; // 转换为0-based
                        size_t segment_end = hide_start - 2;    // 转换为0-based
                        
                        if (segment_start < sequence.length() && segment_end < sequence.length()) {
                            std::string segment = sequence.substr(segment_start, segment_end - segment_start + 1);
                            segments.push_back(segment);
                            segment_coords.push_back({current_pos, hide_start - 1});
                        }
                    }
                    
                    current_pos = hide_end + 1;
                }
                
                // 添加最后一个隐藏区间之后的片段
                if (current_pos <= sequence.length()) {
                    size_t segment_start = current_pos - 1; // 转换为0-based
                    std::string segment = sequence.substr(segment_start);
                    segments.push_back(segment);
                    segment_coords.push_back({current_pos, sequence.length()});
                }
            }
            
            // 6. 写入FASTA文件和FAI记录
            for (size_t i = 0; i < segments.size(); ++i) {
                const std::string& segment = segments[i];
                const auto& coords = segment_coords[i];
                
                if (segment.empty()) continue;
                
                // 生成片段名称
                std::string segment_name;
                if (segments.size() == 1) {
                    segment_name = clean_header;
                } else {
                    segment_name = clean_header + "_segment_" + std::to_string(i + 1) + 
                                  "_" + std::to_string(coords.first) + "-" + std::to_string(coords.second);
                }
                
                // 写入FASTA序列
                out_fasta << ">" << segment_name << "\n";
                for (size_t pos = 0; pos < segment.length(); pos += line_width) {
                    size_t chunk_size = std::min(line_width, segment.length() - pos);
                    out_fasta << segment.substr(pos, chunk_size) << "\n";
                }
                
                // 计算FAI记录信息
                size_t seq_offset = global_offset;
                size_t line_bases = line_width;
                size_t line_bytes = line_width + 1; // 包括换行符
                
                // 写入FAI记录
                out_fai << segment_name << "\t" 
                       << global_start_pos << "\t" 
                       << segment.length() << "\t" 
                       << seq_offset << "\t" 
                       << line_bases << "\t" 
                       << line_bytes << "\n";
                
                // 更新全局位置
                global_start_pos += segment.length();
                
                // 更新文件偏移量
                global_offset += segment_name.length() + 2; // ">" + name + "\n"
                size_t lines_count = (segment.length() + line_width - 1) / line_width;
                global_offset += segment.length() + lines_count; // 序列内容 + 换行符
            }
        }
        
        out_fasta.close();
        out_fai.close();
        
        spdlog::info("Successfully generated hidden intervals FASTA: {}", output_fasta_path.string());
        spdlog::info("Successfully generated corresponding FAI index: {}", output_fai_path.string());
        
        return true;
        
    } catch (const std::exception& e) {
        spdlog::error("Error in hideIntervalsAndGenerateFai: {}", e.what());
        return false;
    }
}

// 根据interval文件隐藏指定区间，生成新的FASTA文件和对应的FAI索引
// interval_file_path: WindowMasker生成的interval文件路径
// output_fasta_path: 输出的隐藏区间后的FASTA文件路径
// output_fai_path: 输出的FAI索引文件路径
// line_width: FASTA文件的行宽
// 返回值: 成功返回true，失败返回false
bool FastaManager::hideIntervalsFromFileAndGenerateFai(
    const FilePath& interval_file_path,
    const FilePath& output_fasta_path,
    const FilePath& output_fai_path,
    size_t line_width) {
    
    try {
        // 1. 解析interval文件，获取需要隐藏的区间
        std::map<std::string, std::vector<std::pair<size_t, size_t>>> intervals_by_seq;
        std::ifstream interval_stream(interval_file_path);
        if (!interval_stream.is_open()) {
            spdlog::error("Failed to open interval file: {}", interval_file_path.string());
            return false;
        }
        
        std::string line;
        std::string current_seq_name;
        while (std::getline(interval_stream, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                current_seq_name = line.substr(1);
                // 清理序列名称中的空白字符
                current_seq_name.erase(0, current_seq_name.find_first_not_of(" \t\n\r\f\v"));
                current_seq_name.erase(current_seq_name.find_last_not_of(" \t\n\r\f\v") + 1);
            } else if (!current_seq_name.empty()) {
                try {
                    std::stringstream ss(line);
                    size_t start, end;
                    char dash;
                    ss >> start >> dash >> end;
                    if (ss.fail() || dash != '-') {
                        spdlog::warn("Invalid interval format: {}", line);
                        continue;
                    }
                    intervals_by_seq[current_seq_name].push_back({start, end});
                } catch (const std::exception& e) {
                    spdlog::warn("Skipping invalid interval line: {} - {}", line, e.what());
                }
            }
        }
        interval_stream.close();
        
        // 2. 调用通用的隐藏区间函数
        return hideIntervalsAndGenerateFai(intervals_by_seq, output_fasta_path, output_fai_path, line_width);
        
    } catch (const std::exception& e) {
        spdlog::error("Error in hideIntervalsFromFileAndGenerateFai: {}", e.what());
        return false;
    }
}

