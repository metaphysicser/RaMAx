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

#endif