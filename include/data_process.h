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
#include "anchor.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
using SpeciesPathMap = std::unordered_map<SpeciesName, FilePath>;
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
    // ��ֹ����
    GzFileWrapper(const GzFileWrapper&) = delete;
    GzFileWrapper& operator=(const GzFileWrapper&) = delete;
private:
    gzFile fp_{ nullptr };
};

// ��װ kseq_t*��ʹ���Զ���ɾ�����Զ����� kseq_destroy
struct KseqDeleter {
    void operator()(kseq_t* seq) const {
        if (seq) kseq_destroy(seq);
    }
};
using KseqPtr = std::unique_ptr<kseq_t, KseqDeleter>;


class FastaManager {
public:
    FilePath fasta_path_;  // FASTA �ļ�·��
    FilePath fai_path_;    // ��ѡ .fai �����ļ�·��
    bool has_n_in_fasta = false;

    // ���캯������������������ RAII��������Ҫ�������������ֶ��ͷ���Դ
    FastaManager() = default;
    FastaManager(const FilePath& fasta_path, const FilePath& fai_path = FilePath())
        : fasta_path_(fasta_path), fai_path_(fai_path)
    {
        if (!fai_path.empty()) {
            loadFaiRecords(fai_path);
        }
        fasta_open(); // �� FASTA �ļ�����ʼ��������
    }
    ~FastaManager() = default; // unique_ptr �Զ��ͷ���Դ

    FastaManager(FastaManager&&) = default;
    FastaManager& operator=(FastaManager&&) = default;

    // ��ȡ��һ�� FASTA ��¼��header �����У�
    bool nextRecord(std::string& header, std::string& sequence);

    // ���ö�ȡ���ļ���ͷ������ RAII �ؽ���������
    void reset();

    // �����м�¼ƴ�ӳ�һ���ַ�����ʹ�� ostringstream �Ż��ַ���ƴ�ӣ�
    std::string concatRecords(char terminator = '\0',
        size_t limit = std::numeric_limits<size_t>::max());

    // FASTA �ļ�ͳ����Ϣ
    struct Stats {
        size_t record_count = 0;
        size_t total_bases = 0;
        size_t min_len = std::numeric_limits<size_t>::max();
        size_t max_len = 0;
        size_t average_len = 0;
    };
    Stats getStats();

    // ��ϴ FASTA ���в�д�����ļ�
    FilePath writeCleanedFasta(const FilePath& output_file, uint64_t line_width = 60);

    // ��ϴ������ FASTA �� FAI �ļ�
    FilePath cleanAndIndexFasta(const FilePath& output_dir,
        const std::string& prefix,
        uint64_t line_width = 60);

    // FAI ��¼���ݽṹ
    struct FaiRecord {
        std::string seq_name;    // �������ƣ����� "chr1"��
        uint_t global_start_pos;
        uint_t length;           // ���г���
        uint_t offset;           // �������ļ��е��ֽ�ƫ����
        uint_t line_bases;       // ÿ�м���������� 60��
        uint_t line_bytes;       // ÿ���ֽ������������з���
    };

    // ʹ�� unordered_map ���� FAI ��¼
    std::vector<FaiRecord> fai_records;

    // ��ָ�������и��ݷ�Χ [start, end] ��ȡ������
    std::string getSubSequence(const std::string& seq_name, size_t start, size_t length);
    std::string getSubConcatSequence(size_t start, size_t length);

    uint_t getConcatSeqLength();
    ChrName getChrName(uint_t global_start, uint_t length);

    RegionVec preAllocateChunks(uint_t chunk_size, uint_t overlap_size);

private:
    // RAII ��װ�󣬲���ֱ��ʹ��ԭʼָ��
    std::unique_ptr<GzFileWrapper> gz_file_wrapper_;
    KseqPtr seq_;

    bool clean_data_{ false };

    // ������ϴ������ ACGTN ���ַ��滻Ϊ N��ͬʱ��¼�Ƿ���� N
    void cleanSequence(std::string& seq);

    // �� FASTA �ļ��ͳ�ʼ�� kseq ������
    void fasta_open();

    // ����ɨ�� FASTA ��д�� .fai �ļ�
    bool reScanAndWriteFai(const FilePath& fa_path,
        const FilePath& fai_path,
        size_t line_width) const;

    // ���� .fai �ļ��еļ�¼�� fai_records
    void loadFaiRecords(const FilePath& fai_path);
};

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

bool isFileSmallerThan(const FilePath& filePath, size_t maxSizeMB);

//bool splitRawDataToChr(const FilePath workdir_path, SpeciesPathMap& species_path_map, SpeciesChrPathMap& species_chr_path_map, int thread_num);
//
//bool splitChrToChunk(FilePath work_dir, SpeciesChrPathMap& species_chr_path_map, SpeciesChunkInfoMap& species_chunk_info_map, uint_t chunk_length, uint_t overlap_length, int thread_num);

FilePath getTempFilePath(const FilePath& input_path);

FilePath getFaiIndexPath(const FilePath& fasta_path);
#endif // RAW_DATA_PROCESS_HPP
