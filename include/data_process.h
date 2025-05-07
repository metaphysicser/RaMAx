#ifndef DATA_PROCESS_H
#define DATA_PROCESS_H

#include <regex>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cctype>
#include <curl/curl.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <zlib.h>
#include "config.hpp"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

struct ChunkInfo {
    ChrName chr_name;  // 染色体／序列名称
    uint_t start;          // 在该染色体上的起始偏移（0-based）
    uint_t length;         // 本 chunk 的长度（最后一个可能不足 chunk_size）
};
using ChunkInfoVec = std::vector<ChunkInfo>;

class GzFileWrapper {
public:
    explicit GzFileWrapper(const std::string& filename, const std::string& mode = "r") {
        fp_ = gzopen(filename.c_str(), mode.c_str());
        if (!fp_) {
            throw std::runtime_error("Failed to open FASTA file: " + filename);
        }
    }
    ~GzFileWrapper() {
        if (fp_) {
            gzclose(fp_);
        }
    }
    gzFile get() const { return fp_; }
    // 禁止复制
    GzFileWrapper(const GzFileWrapper&) = delete;
    GzFileWrapper& operator=(const GzFileWrapper&) = delete;
private:
    gzFile fp_{ nullptr };
};

// 封装 kseq_t*，使用自定义删除器自动调用 kseq_destroy
struct KseqDeleter {
    void operator()(kseq_t* seq) const {
        if (seq) kseq_destroy(seq);
    }
};
using KseqPtr = std::unique_ptr<kseq_t, KseqDeleter>;

class FastaManager {
public:
    FilePath fasta_path_;  // FASTA 文件路径
    FilePath fai_path_;    // 可选 .fai 索引文件路径
    bool has_n_in_fasta = false;

    // 构造函数与析构函数：采用 RAII，不再需要在析构函数中手动释放资源
    FastaManager() = default;
    FastaManager(const FilePath& fasta_path, const FilePath& fai_path = FilePath())
        : fasta_path_(fasta_path), fai_path_(fai_path)
    {
        if (!fai_path.empty()) {
            loadFaiRecords(fai_path);
        }
        fasta_open(); // 打开 FASTA 文件并初始化解析器
    }
    ~FastaManager() = default; // unique_ptr 自动释放资源

    FastaManager(FastaManager&&) = default;
    FastaManager& operator=(FastaManager&&) = default;

    // 读取下一个 FASTA 记录（header 与序列）
    bool nextRecord(std::string& header, std::string& sequence);

    // 重置读取到文件开头（采用 RAII 重建解析器）
    void reset();

    // 将所有记录拼接成一个字符串（使用 ostringstream 优化字符串拼接）
    std::string concatRecords(char separator = '0', char terminator = '1',
        size_t limit = std::numeric_limits<size_t>::max());

    // FASTA 文件统计信息
    struct Stats {
        size_t record_count = 0;
        size_t total_bases = 0;
        size_t min_len = std::numeric_limits<size_t>::max();
        size_t max_len = 0;
        size_t average_len = 0;
    };
    Stats getStats();

    // 清洗 FASTA 序列并写入新文件
    FilePath writeCleanedFasta(const FilePath& output_file, uint64_t line_width = 60);

    // 清洗并生成 FASTA 与 FAI 文件
    FilePath cleanAndIndexFasta(const FilePath& output_dir,
        const std::string& prefix,
        uint64_t line_width = 60);

    // FAI 记录数据结构
    struct FaiRecord {
        std::string seq_name;    // 序列名称（例如 "chr1"）
        uint_t global_start_pos;
        uint_t length;           // 序列长度
        uint_t offset;           // 序列在文件中的字节偏移量
        uint_t line_bases;       // 每行碱基数（例如 60）
        uint_t line_bytes;       // 每行字节数（包括换行符）
    };

    // 使用 unordered_map 保存 FAI 记录
    std::vector<FaiRecord> fai_records;

    // 从指定序列中根据范围 [start, end] 提取子序列
    std::string getSubSequence(const std::string& seq_name, size_t start, size_t length);
    std::string getSubConcatSequence(size_t start, size_t length);

    uint_t getConcatSeqLength();
    ChrName getChrName(uint_t global_start, uint_t length);

    ChunkInfoVec preAllocateChunks(uint_t chunk_size, uint_t overlap_size);

private:
    // RAII 封装后，不再直接使用原始指针
    std::unique_ptr<GzFileWrapper> gz_file_wrapper_;
    KseqPtr seq_;

    bool clean_data_{ false };

    // 序列清洗：将非 ACGTN 的字符替换为 N，同时记录是否存在 N
    void cleanSequence(std::string& seq);

    // 打开 FASTA 文件和初始化 kseq 解析器
    void fasta_open();

    // 重新扫描 FASTA 并写入 .fai 文件
    bool reScanAndWriteFai(const FilePath& fa_path,
        const FilePath& fai_path,
        size_t line_width) const;

    // 加载 .fai 文件中的记录到 fai_records
    void loadFaiRecords(const FilePath& fai_path);
};


using SpeciesPathMap = std::unordered_map<SpeciesName, FilePath>;
//using ChrPathMap = std::unordered_map<Chr, FilePath>;
//using SpeciesChrPathMap = std::unordered_map<Chr, ChrPathMap>;
//using ChunkInfoVec = std::vector<ChunkInfo>;
//using ChrChunkInfoMap = std::unordered_map<Chr, ChunkInfoVec>;
//using SpeciesChunkInfoMap = std::unordered_map<Species, ChrChunkInfoMap>;

//bool saveSpeciesChunkInfoMap(const SpeciesChunkInfoMap& species_chunk_info_map, const std::string& filename);
//
//bool loadSpeciesChunkInfoMap(SpeciesChunkInfoMap& species_chunk_info_map, const std::string& filename);


// Check if a string is a URL
bool isUrl(const std::string& path_str);

// Verify if a URL is reachable
bool verifyUrlReachable(const std::string& url);

// Verify if a local file exists and is valid
void verifyLocalFile(const FilePath& file_path);

// Callback function for CURL write operation
size_t writeData(void* ptr, size_t size, size_t nmemb, void* stream);

// Download a file from a URL using libcurl
void downloadFile(const std::string& url, const FilePath& destination);

// Copy a local file to the destination
void copyLocalFile(const FilePath& source, const FilePath& destination);

// Get the file extension from a path
std::string getFileExtension(const FilePath& file_path);

// Main function to copy or download raw data
bool copyRawData(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num);

bool cleanRawDataset(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num);

// Clean the raw data file by converting to uppercase and replacing invalid characters.
// The cleaned file is saved to workdir_path / DATA_DIR / CLEAN_DATA, preserving the original filename.
// If the target file exists, the process is skipped. A temporary file is used for in-process data.
// FilePath cleanRawData(const FilePath workdir_path, const FilePath raw_data_path);

// Get human-readable file size (auto convert to KB/MB/GB). Supports both local files and URLs.
std::string getReadableFileSize(const FilePath& filePath);

//bool splitRawDataToChr(const FilePath workdir_path, SpeciesPathMap& species_path_map, SpeciesChrPathMap& species_chr_path_map, int thread_num);
//
//bool splitChrToChunk(FilePath work_dir, SpeciesChrPathMap& species_chr_path_map, SpeciesChunkInfoMap& species_chunk_info_map, uint_t chunk_length, uint_t overlap_length, int thread_num);

FilePath getTempFilePath(const FilePath& input_path);

FilePath getFaiIndexPath(const FilePath& fasta_path);
#endif // RAW_DATA_PROCESS_HPP
