#ifndef DATA_PROCESS_H
#define DATA_PROCESS_H

// -----------------------------
// 包含依赖头文件
// -----------------------------
#include <regex>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cctype>
#include <curl/curl.h>      // 用于 URL 请求
#include <fcntl.h>
#include <sys/mman.h>       // 用于内存映射读取文件
#include <zlib.h>           // 用于处理 gzip 压缩格式
#include "config.hpp"
#include "anchor.h"
#include "kseq.h"           // 用于解析 FASTA 格式

// 初始化 kseq 使用 gzFile 类型
KSEQ_INIT(gzFile, gzread)

// 物种名到文件路径的映射
using SpeciesPathMap = std::unordered_map<SpeciesName, FilePath>;

// -----------------------------
// GzFileWrapper：RAII 封装 gzopen/gzclose
// -----------------------------
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

    // 禁止复制操作（只允许移动或指针管理）
    GzFileWrapper(const GzFileWrapper&) = delete;
    GzFileWrapper& operator=(const GzFileWrapper&) = delete;

private:
    gzFile fp_{ nullptr };
};

// -----------------------------
// kseq_t 的智能指针删除器
// -----------------------------
struct KseqDeleter {
    void operator()(kseq_t* seq) const {
        if (seq) kseq_destroy(seq);
    }
};

using KseqPtr = std::unique_ptr<kseq_t, KseqDeleter>;

// -----------------------------
// FastaManager 类：用于读取、清洗、索引 fasta 序列数据
// -----------------------------
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

    // 获取下一个序列记录（返回 false 表示读完）
    bool nextRecord(std::string& header, std::string& sequence);

    // 重置读取状态，从头开始读取
    void reset();

    // 将多个记录拼接为一个序列（用于压缩或合并）
    std::string concatRecords(char terminator = '\0',
        size_t limit = std::numeric_limits<size_t>::max());

    // 获取序列文件的基本统计信息
    struct Stats {
        size_t record_count = 0;
        size_t total_bases = 0;
        size_t min_len = std::numeric_limits<size_t>::max();
        size_t max_len = 0;
        size_t average_len = 0;
    };
    Stats getStats();

    // 将清洗后的序列写入新的 fasta 文件
    FilePath writeCleanedFasta(const FilePath& output_file, uint64_t line_width = 60);

    // 清洗并生成 .fai 索引，返回新 fasta 文件路径
    FilePath cleanAndIndexFasta(const FilePath& output_dir,
        const std::string& prefix,
        uint64_t line_width = 60);

    // -----------------------------
    // FAI 索引记录结构
    // -----------------------------
    struct FaiRecord {
        std::string seq_name;        // 序列名（如 chr1）
        uint_t global_start_pos;    // 该序列在拼接序列中的全局起始坐标
        uint_t length;              // 碱基总数
        uint_t offset;              // 在文件中的偏移量（跳过注释）
        uint_t line_bases;          // 每行的碱基数
        uint_t line_bytes;          // 每行占用的字节数（含换行）
    };

    // 所有 FAI 记录
    std::vector<FaiRecord> fai_records;

    // 从指定序列中获取区间子串（基于名称 + 坐标）
    std::string getSubSequence(const std::string& seq_name, size_t start, size_t length);

    // 从拼接序列的全局坐标中提取连续片段
    std::string getSubConcatSequence(size_t start, size_t length);

    // 获取拼接后的总序列长度
    uint_t getConcatSeqLength();

    // 从全局坐标范围反查染色体名称
    ChrName getChrName(uint_t global_start, uint_t length);

    // 按 chunk_size 和 overlap_size 对所有序列预分段
    RegionVec preAllocateChunks(uint_t chunk_size, uint_t overlap_size);

private:
    // 内部成员变量：gzip 文件包装和 kseq_t 管理器
    std::unique_ptr<GzFileWrapper> gz_file_wrapper_;
    KseqPtr seq_;

    bool clean_data_{ false }; // 是否启用清洗

    // 清洗序列（全部转大写 + 过滤非法字符）
    void cleanSequence(std::string& seq);

    // 打开 fasta 文件并初始化 kseq_t
    void fasta_open();

    // 若不存在 .fai 索引则重新扫描并生成
    bool reScanAndWriteFai(const FilePath& fa_path,
        const FilePath& fai_path,
        size_t line_width) const;

    // 加载已有 .fai 文件
    void loadFaiRecords(const FilePath& fai_path);
};

// -----------------------------
// 工具函数声明
// -----------------------------

// 判断字符串是否是 URL 格式（http/https/ftp）
bool isUrl(const std::string& path_str);

// 验证 URL 是否可以访问
bool verifyUrlReachable(const std::string& url);

// 验证本地文件路径是否存在且合法
void verifyLocalFile(const FilePath& file_path);

// CURL 下载文件时的写入回调函数
size_t writeData(void* ptr, size_t size, size_t nmemb, void* stream);

// 使用 curl 下载远程文件到本地
void downloadFile(const std::string& url, const FilePath& destination);

// 将本地文件复制到目标路径
void copyLocalFile(const FilePath& source, const FilePath& destination);

// 提取文件扩展名字符串
std::string getFileExtension(const FilePath& file_path);

// 拷贝或下载原始数据至工作目录
bool copyRawData(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num);

// 清洗原始数据集，并更新路径映射
bool cleanRawDataset(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num);

// 获取文件大小（支持本地和 URL），输出带单位字符串
std::string getReadableFileSize(const FilePath& filePath);

// 判断文件是否小于指定大小（单位 MB）
bool isFileSmallerThan(const FilePath& filePath, size_t maxSizeMB = 1024);

// 生成一个临时中间处理路径（如 _in_process.fasta）
FilePath getTempFilePath(const FilePath& input_path);

// 获取 .fai 索引的标准路径（原路径加 .fai 后缀）
FilePath getFaiIndexPath(const FilePath& fasta_path);

#endif // DATA_PROCESS_H
