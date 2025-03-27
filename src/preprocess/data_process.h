#ifndef DATA_PROCESS_H
#define DATA_PROCESS_HP

#include <regex>
#include <filesystem>
#include <unordered_map>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cctype>
#include <curl/curl.h>
#include <zlib.h>
#include "config.hpp"
#include "threadpool.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


// Type alias for mapping species names to file paths
using FilePath = std::filesystem::path;
using Species = std::string;
using Chr = std::string;

struct ChunkInfo {
	FilePath file_path;
	uint_t chunk_index;
	std::string species;
	std::string chr;
	uint_t start;
	uint_t end;
	uint_t length;

	template <class Archive>
	void serialize(Archive& ar) {
		ar(
			CEREAL_NVP(file_path),
			CEREAL_NVP(chunk_index),
			CEREAL_NVP(species),
			CEREAL_NVP(chr),
			CEREAL_NVP(start),
			CEREAL_NVP(end),
			CEREAL_NVP(length)
		);
	}
};



using SpeciesPathMap = std::unordered_map<Species, FilePath>;
using ChrPathMap = std::unordered_map<Chr, FilePath>;
using SpeciesChrPathMap = std::unordered_map<Chr, ChrPathMap>;
using ChunkInfoVec = std::vector<ChunkInfo>;
using ChrChunkInfoMap = std::unordered_map<Chr, ChunkInfoVec>;
using SpeciesChunkInfoMap = std::unordered_map<Species, ChrChunkInfoMap>;

bool saveSpeciesChunkInfoMap(const SpeciesChunkInfoMap& species_chunk_info_map, const std::string& filename);

bool loadSpeciesChunkInfoMap(SpeciesChunkInfoMap& species_chunk_info_map, const std::string& filename);




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
FilePath cleanRawData(const FilePath workdir_path, const FilePath raw_data_path);

// Get human-readable file size (auto convert to KB/MB/GB). Supports both local files and URLs.
std::string getReadableFileSize(const FilePath& filePath);

bool splitRawDataToChr(const FilePath workdir_path, SpeciesPathMap& species_path_map, SpeciesChrPathMap& species_chr_path_map, int thread_num);


bool splitChrToChunk(FilePath work_dir, SpeciesChrPathMap& species_chr_path_map, SpeciesChunkInfoMap& species_chunk_info_map, uint_t chunk_length, uint_t overlap_length, int thread_num);
#endif // RAW_DATA_PROCESS_HPP
