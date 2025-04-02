#ifndef DATA_PROCESS_H
#define DATA_PROCESS_H

#include <regex>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cctype>
#include <curl/curl.h>
#include <zlib.h>
#include "config.hpp"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

// Type alias for mapping species names to file paths

//struct ChunkInfo {
//	FilePath file_path;
//	uint_t chunk_index;
//	std::string species;
//	std::string chr;
//	uint_t start;
//	uint_t end;
//	uint_t length;
//
//	template <class Archive>
//	void serialize(Archive& ar) {
//		ar(
//			CEREAL_NVP(file_path),
//			CEREAL_NVP(chunk_index),
//			CEREAL_NVP(species),
//			CEREAL_NVP(chr),
//			CEREAL_NVP(start),
//			CEREAL_NVP(end),
//			CEREAL_NVP(length)
//		);
//	}
//};

class FastaManager {
public:
    FilePath fasta_path_;  // Path to .fasta file
    FilePath fai_path_;    // Optional path to .fai file
    /// Constructor with optional FAI index path
    /// If fai_path is not provided, only sequence reading is enabled
    FastaManager(const FilePath& fasta_path, const FilePath& fai_path = FilePath());
    ~FastaManager();

    /// Read the next FASTA record (header and sequence)
    bool nextRecord(std::string& header, std::string& sequence);

    /// Reset reading to the beginning of the file
    void reset();

    /// Concatenate all records into a single sequence
    /// with optional separator and terminator characters
    std::string concatRecords(char separator = '0', char terminator = '1',
        size_t limit = std::numeric_limits<size_t>::max());

    /// Statistics of the FASTA file (record count, lengths, etc.)
    struct Stats {
        size_t record_count = 0;
        size_t total_bases = 0;
        size_t min_len = SIZE_MAX;
        size_t max_len = 0;
        size_t average_len = 0;
    };
    Stats getStats();

    /// Write cleaned sequences to a new FASTA file (no index)
    FilePath writeCleanedFasta(const FilePath& output_dir, uint64_t line_width = 60);

    /// Clean sequences and generate both FASTA and FAI files
    FilePath cleanAndIndexFasta(const FilePath& output_dir,
        const std::string& prefix,
        uint64_t line_width = 60);

    /// FAI record (sequence index metadata)
    struct FaiRecord {
        std::string seq_name;    // Sequence name (e.g., "chr1")
        size_t length;           // Total number of bases
        size_t offset;           // Byte offset of the first base in the file
        size_t line_bases;       // Number of bases per line (e.g., 60)
        size_t line_bytes;       // Number of bytes per line (including newline)
    };

    /// Loaded FAI records from .fai file
    std::unordered_map<std::string, FaiRecord> fai_records;

    /// Extract a sub-sequence from a given sequence by name and range [start, end]
    std::string getSubSequence(const std::string& seq_name, size_t start, size_t end);

private:
    gzFile fp_{ nullptr };  // File pointer (used for compressed or uncompressed FASTA)
    kseq_t* seq_{ nullptr };  // KSEQ parser handle

    bool clean_data_{ false };

    void cleanSequence(std::string& seq);  // Optional cleaning logic (replace non-ACGTN)

    void open();  // (Re)open fasta file and initialize parser

    /// Scan fasta and write .fai index (used after writing cleaned fasta)
    bool reScanAndWriteFai(const FilePath& fa_path,
        const FilePath& fai_path,
        size_t line_width) const;

    /// Load .fai file into fai_records
    void loadFaiRecords(const FilePath& fai_path);
};

using SpeciesPathMap = std::unordered_map<Species, FilePath>;
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
