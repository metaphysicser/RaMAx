#ifndef SPLIT_DATA_HPP
#define SPLIT_DATA_HPP

#include "data_process.h"


// Function to split genome data into chromosome-level files
bool splitRawDataToChr(const FilePath workdir_path,
	SpeciesPathMap& species_path_map,
	SpeciesChrPathMap& species_chr_path_map,
	int thread_num) {
	try {
		// Create "split_chr" directory if it doesn't exist
		FilePath split_chr_dir = workdir_path / DATA_DIR / SPLIT_CHR_DIR;
		if (!std::filesystem::exists(split_chr_dir)) {
			std::filesystem::create_directories(split_chr_dir);
			spdlog::info("Created split_chr directory: {}", split_chr_dir.string());
		}
		else {
			spdlog::warn("split_chr directory already exists: {}", split_chr_dir.string());
		}

		// Create a thread pool for parallel processing of multiple species
		ThreadPool pool(thread_num);

		// Loop through each species in the species_path_map
		for (auto it = species_path_map.begin(); it != species_path_map.end(); ++it) {
			const std::string& species = it->first;  // Species name
			const FilePath& fasta_path = it->second; // Path to the genome file (FASTA format)

			// Enqueue a task to process each species
			pool.enqueue([&species, &fasta_path, &split_chr_dir, &species_chr_path_map]() {
				try {
					// Create a directory for the current species
					FilePath species_dir = split_chr_dir / species;
					if (!std::filesystem::exists(species_dir)) {
						std::filesystem::create_directories(species_dir);
						spdlog::info("Created directory for species: {}", species_dir.string());
					}
					else {
						spdlog::warn("Directory for species {} already exists: {}", species, species_dir.string());
					}

					// Map to store the chromosome paths for the species
					std::unordered_map<std::string, FilePath> chr_path_map;

					// Open the FASTA file (supporting both .gz and .fasta formats)
					gzFile fp = gzopen(fasta_path.string().c_str(), "r");
					if (!fp) {
						throw std::runtime_error("Failed to open FASTA file: " + fasta_path.string());
					}

					// Initialize kseq to read the FASTA file
					kseq_t* seq = kseq_init(fp);
					int ret;
					int count = 0;

					// Read sequences from the FASTA file
					while ((ret = kseq_read(seq)) >= 0) {
						std::string header(seq->name.s);  // Chromosome header (e.g., "chr1")
						std::string sequence(seq->seq.s); // Sequence data

						// Path to save the chromosome sequence file
						FilePath chr_file = species_dir / (header + ".fasta");
						count++; // Count the number of chromosomes processed

						// Record the chromosome file path in the map
						chr_path_map[header] = chr_file;

						// Skip if the chromosome file already exists
						if (std::filesystem::exists(chr_file)) {
							spdlog::warn("File already exists, skipping: {}", chr_file.string());
							continue;
						}

						// Write the chromosome sequence to the file
						std::ofstream chr_ofs(chr_file, std::ios::out);
						chr_ofs << ">" << header << std::endl;
						chr_ofs << sequence << std::endl;
						chr_ofs.close();
					}

					// Store the chromosome paths for the species in the species_chr_path_map
					species_chr_path_map[species] = chr_path_map;
					spdlog::info("Processed {} chromosomes for species {}, saved to: {}", count, species, species_dir.string());
				}
				catch (const std::exception& e) {
					spdlog::error("Error processing species {}: {}", species, e.what());
					throw std::runtime_error("Failed to split genome data for species: " + species);
				}
				});
		}

		// Wait for all tasks in the thread pool to finish
		pool.waitAllTasksDone();
		spdlog::info("Successfully split genomes into chromosomes!");
		return true;
	}
	catch (const std::exception& e) {
		spdlog::error("Error in splitRawData: {}", e.what());
		throw std::runtime_error("Failed to split genome data into chromosomes.");
		return false;
	}
}

bool splitChrToChunk(FilePath work_dir,
	SpeciesChrPathMap& species_chr_path_map,
	SpeciesChunkInfoMap& species_chunk_info_map,
	uint_t chunk_length,
	uint_t overlap_length,
	int thread_num) {

	FilePath chunk_map_filename = work_dir / DATA_DIR / CHUNK_DIR / CHUNK_MAP_FILE;

	if (std::filesystem::exists(chunk_map_filename)) {
		spdlog::info("Found existing chunk info map file: {}. Loading...", chunk_map_filename.string());
		if (!loadSpeciesChunkInfoMap(species_chunk_info_map, chunk_map_filename)) {
			spdlog::error("Failed to load species_chunk_info_map from file: {}", chunk_map_filename.string());
			return false;
		}
		spdlog::info("Successfully loaded chunk info map. Skipping chunk splitting.");
		return true;
	}

	ThreadPool pool(thread_num);

	// Iterate over each species in the species_chr_path_map
	// (Assuming SpeciesChrPathMap maps from Species to ChrPathMap)
	for (auto species_it = species_chr_path_map.begin(); species_it != species_chr_path_map.end(); ++species_it) {
		Species species = species_it->first;
		ChrPathMap& chr_map = species_it->second;

		// Iterate over each chromosome in the ChrPathMap
		for (auto chr_it = chr_map.begin(); chr_it != chr_map.end(); ++chr_it) {
			Chr chr = chr_it->first;
			FilePath chr_file = chr_it->second;

			// Enqueue a task to split the chromosome file into chunks
			pool.enqueue([species, chr, chr_file, &species_chunk_info_map, work_dir, chunk_length, overlap_length]() {
				// Create output directory for chunks: work_dir/CHUNK_DIR/<species>/<chr>
				FilePath out_dir = work_dir / DATA_DIR / CHUNK_DIR / species / chr;
				if (!std::filesystem::exists(out_dir)) {
					std::filesystem::create_directories(out_dir);
				}

				// Open the chromosome file using kseq (supports gzipped FASTA)
				gzFile fp = gzopen(chr_file.string().c_str(), "r");
				if (!fp) {
					spdlog::error("Failed to open chromosome file: {}", chr_file.string());
					return;
				}
				kseq_t* seq = kseq_init(fp);
				int_t l = kseq_read(seq);
				if (l < 0) {
					spdlog::error("Failed to read sequence from file: {}", chr_file.string());
					kseq_destroy(seq);
					gzclose(fp);
					return;
				}
				// Assume each chromosome file contains a single FASTA record
				std::string sequence(seq->seq.s, seq->seq.l);
				kseq_destroy(seq);
				gzclose(fp);

				uint_t seq_length = sequence.size();
				ChunkInfoVec chunks;
				uint_t chunk_index = 0;

				// Slide through the sequence to generate chunks
				for (uint_t start = 0; start < seq_length; ) {
					uint_t end = start + chunk_length;
					if (end > seq_length)
						end = seq_length;
					uint_t length = end - start;

					std::string filename = chr + "_chunk_" + std::to_string(chunk_index) + ".fasta";
					FilePath chunk_file = out_dir / filename;

					// If target chunk file does not exist, create it via a temporary file
					if (!std::filesystem::exists(chunk_file)) {
						// Construct temporary file path: add "_in_process" before extension
						FilePath temp_chunk_file;
						if (chunk_file.has_extension()) {
							std::string stem = chunk_file.stem().string();
							std::string extension = chunk_file.extension().string();
							temp_chunk_file = chunk_file.parent_path() / (stem + "_in_process" + extension);
						}
						else {
							temp_chunk_file = chunk_file.parent_path() / (chunk_file.filename().string() + "_in_process");
						}

						// If temporary file exists, remove it
						if (std::filesystem::exists(temp_chunk_file)) {
							spdlog::warn("Temporary chunk file exists, removing: {}", temp_chunk_file.string());
							std::filesystem::remove(temp_chunk_file);
						}

						// Write the chunk to the temporary FASTA file
						std::ofstream ofs(temp_chunk_file);
						if (!ofs) {
							spdlog::error("Failed to write temporary chunk file: {}", temp_chunk_file.string());
							// Continue with next chunk even if one fails
							start = (end == seq_length) ? seq_length : start + (chunk_length - overlap_length);
							++chunk_index;
							continue;
						}
						ofs << ">" << chr << "_chunk_" << chunk_index
							<< "_start_" << start << "_end_" << end << "\n";
						ofs << sequence.substr(start, end - start) << "\n";
						ofs.close();

						// Rename temporary file to target chunk file
						std::filesystem::rename(temp_chunk_file, chunk_file);
					}
					// Else, if the target file already exists, skip writing

					// Record the chunk information
					ChunkInfo info;
					info.file_path = chunk_file;
					info.species = species;
					info.chunk_index = chunk_index;
					info.chr = chr;
					info.start = start;
					info.end = end;
					info.length = length;
					chunks.push_back(info);

					// Compute next start: move forward by (chunk_length - overlap_length)
					if (end == seq_length)
						break;
					start = start + (chunk_length - overlap_length);
					++chunk_index;
				}

				// Save generated chunks into the shared species_chunk_info_map
				species_chunk_info_map[species][chr] = chunks;
				spdlog::info("Species {}: {} chunks generated for chromosome {}", species, chunks.size(), chr);
				});
		}
	}

	// Wait for all tasks to complete
	pool.waitAllTasksDone();

	spdlog::info("Saving species_chunk_info_map to file: {}", chunk_map_filename.string());
	if (!saveSpeciesChunkInfoMap(species_chunk_info_map, chunk_map_filename)) {
		spdlog::error("Failed to save species_chunk_info_map to file: {}", chunk_map_filename.string());
		return false;
	}
	spdlog::info("Successfully saved species_chunk_info_map.");
	return true;
}



#endif