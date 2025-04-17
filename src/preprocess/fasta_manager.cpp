#include "data_process.h"

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

void FastaManager::loadFaiRecords(const FilePath& fai_path)
{
    // Open the .fai file
    std::ifstream in(fai_path);
    if (!in) {
        spdlog::error("Failed to open fai file: {}", fai_path.string());
        throw std::runtime_error("Cannot open " + fai_path.string());
    }

    std::string line;
    
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
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // Parse the FAI line format: name length offset line_bases line_bytes
        std::istringstream iss(line);

        FaiRecord rec;
        iss >> rec.seq_name >> rec.length >> rec.offset >> rec.line_bases >> rec.line_bytes;
        if (iss.fail()) {
            spdlog::warn("Skipping malformed line in {}: {}", fai_path.string(), line);
            continue;
        }

        // Insert the record into the unordered_map, using seq_name as the key
        fai_records.push_back(rec);
    }

    in.close();
}

bool FastaManager::nextRecord(std::string& header, std::string& sequence) {
    int ret = kseq_read(seq_.get());
    if (ret < 0) {
        return false; // 无更多记录
    }
    header.assign(seq_->name.s);
    sequence.assign(seq_->seq.s, seq_->seq.l);

    if (clean_data_) {
        cleanSequence(sequence);
    }
    return true;
}

void FastaManager::reset() {
    // 通过 unique_ptr 的 reset 释放旧资源，然后重新打开
    gz_file_wrapper_.reset();
    seq_.reset();
    fasta_open();
}

std::string FastaManager::concatRecords(char separator, char terminator, size_t limit) {
    reset();
    std::ostringstream oss;
    std::string hdr, seq;
    size_t count = 0;

    while (nextRecord(hdr, seq) && count < limit) {
        oss << seq << separator;  // 使用 ostringstream 拼接，减少不必要的内存拷贝
        ++count;
    }
    std::string result = oss.str();
    if (!result.empty()) {
        result.back() = terminator; // 替换最后一个分隔符为终止符
    }
    return result;
}

typename FastaManager::Stats FastaManager::getStats() {
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

void FastaManager::cleanSequence(std::string& seq) {
    for (char& c : seq) {
        unsigned char uc = static_cast<unsigned char>(c);
        uc = std::toupper(uc);
        if (!has_n_in_fasta && uc == 'N') {
            has_n_in_fasta = true;
        }
        // 可加入 OpenMP 并行化优化（需要确定线程安全）
        if (uc != 'A' && uc != 'C' && uc != 'G' && uc != 'T' && uc != 'N') {
            uc = 'N';
        }
        c = static_cast<char>(uc);
    }
}

void FastaManager::fasta_open() {
    // 使用 RAII 封装自动管理 gzFile 和 kseq_t 的生命周期
    gz_file_wrapper_ = std::make_unique<GzFileWrapper>(fasta_path_.string());
    seq_.reset(kseq_init(gz_file_wrapper_->get()));
}

bool FastaManager::reScanAndWriteFai(const FilePath& fa_path,
    const FilePath& fai_path,
    size_t line_width) const
{
    if (std::filesystem::exists(fai_path)) {
        spdlog::warn("Fai file already exists: {}", fai_path.string());
        return true;
    }
    std::ifstream in(fa_path);
    if (!in) {
        spdlog::error("Failed to open {} for reading in reScanAndWriteFai", fa_path.string());
        return false;
    }

    FilePath tmp_fai_path = getTempFilePath(fai_path);
    std::ofstream out(tmp_fai_path);
    if (!out) {
        spdlog::error("Failed to open {} for writing .fai", fai_path.string());
        return false;
    }

    out << (has_n_in_fasta ? "YES" : "NO") << "\n";

    size_t global_offset = 0;
    std::string line;
    std::string seq_name;
    size_t seq_len = 0;
    size_t seq_start = 0;
    bool reading_seq = false;

    size_t line_bytes_in_fasta = line_width + 1;

    while (std::getline(in, line)) {
        size_t this_line_bytes = line.size() + 1;

        if (!line.empty() && line[0] == '>') {
            if (reading_seq) {
                out << seq_name << "\t" << seq_len << "\t"
                    << seq_start << "\t"
                    << line_width << "\t"
                    << line_bytes_in_fasta << "\n";
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
        out << seq_name << "\t" << seq_len << "\t"
            << seq_start << "\t"
            << line_width << "\t"
            << line_bytes_in_fasta << "\n";
    }

    in.close();
    out.close();

    std::filesystem::rename(tmp_fai_path, fai_path);
    spdlog::info("FAI index created: {}", fai_path.string());
    return true;
}

FilePath FastaManager::writeCleanedFasta(const FilePath& output_file, uint64_t line_width) {
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

    // 将文件映射到内存
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
    return total_length + 1;
}
