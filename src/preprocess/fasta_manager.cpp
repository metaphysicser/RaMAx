#include "data_process.h"

FastaManager::FastaManager(const FilePath& fasta_path, const FilePath& fai_path)
    : fasta_path_(fasta_path)
    , fai_path_(fai_path)
{
    // If fai_path is not empty, load the FAI records into an unordered_map
    if (!fai_path.empty()) {
        loadFaiRecords(fai_path);
    }
    open(); // Open FASTA file
}

FastaManager::~FastaManager() {
    if (seq_) {
        kseq_destroy(seq_);
        seq_ = nullptr;
    }
    if (fp_) {
        gzclose(fp_);
        fp_ = nullptr;
    }
}

void FastaManager::loadFaiRecords(const FilePath& fai_path)
{
    // Open the .fai file
    std::ifstream in(fai_path);
    if (!in) {
        spdlog::error("Failed to open fai file: {}", fai_path.string());
        throw std::runtime_error("Cannot open " + fai_path.string());
    }

    std::string line;
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
        fai_records[rec.seq_name] = rec;
    }

    in.close();
}

bool FastaManager::nextRecord(std::string& header, std::string& sequence) {
    int ret = kseq_read(seq_);
    if (ret < 0) {
        return false; // No more records
    }
    header.assign(seq_->name.s);
    sequence.assign(seq_->seq.s, seq_->seq.l);

    if (clean_data_) {
        cleanSequence(sequence);
    }
    return true;
}

void FastaManager::reset() {
    if (seq_) {
        kseq_destroy(seq_);
        seq_ = nullptr;
    }
    if (fp_) {
        gzclose(fp_);
        fp_ = nullptr;
    }
    open();
}

std::string FastaManager::concatRecords(char separator,
    char terminator,
    size_t limit)
{
    reset();
    std::string result;
    std::string hdr, seq;
    size_t count = 0;

    while (nextRecord(hdr, seq) && count < limit) {
        result += seq;
        // result.push_back(separator);
        ++count;
    }

    if (!result.empty()) {
        result.back() = 0x00; // Replace last separator
    }
    return result;
}

FastaManager::Stats FastaManager::getStats() {
    reset();
    Stats s;
    std::string hdr, seq;
    while (nextRecord(hdr, seq)) {
        size_t len = seq.size();
        s.record_count++;
        s.total_bases += len;
        if (len < s.min_len) {
            s.min_len = len;
        }
        if (len > s.max_len) {
            s.max_len = len;
        }
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
        if (uc != 'A' && uc != 'C' && uc != 'G' && uc != 'T' && uc != 'N') {
            uc = 'N';
        }
        c = static_cast<char>(uc);
    }
}

void FastaManager::open() {
    fp_ = gzopen(fasta_path_.string().c_str(), "r");
    if (!fp_) {
        throw std::runtime_error("Failed to open FASTA file: " + fasta_path_.string());
    }
    seq_ = kseq_init(fp_);
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

FilePath FastaManager::writeCleanedFasta(const FilePath& output_file, uint64_t line_width)
{
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

std::string FastaManager::getSubSequence(const std::string& seq_name,
    size_t start,
    size_t end)
{
    // Use the unordered_map to quickly find the FAI record by seq_name
    auto it = fai_records.find(seq_name);
    if (it == fai_records.end()) {
        throw std::runtime_error("Sequence " + seq_name + " not found in FAI index.");
    }

    const FaiRecord& rec = it->second;

    // Check for valid range
    if (start > end || end >= rec.length) {
        throw std::runtime_error("Invalid range [" + std::to_string(start) +
            "," + std::to_string(end) + "] for sequence " + seq_name);
    }

    std::ifstream in(fasta_path_, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open " + fasta_path_.string() + " for reading sub-sequence.");
    }

    size_t req_len = end - start + 1;
    std::string result;
    result.reserve(req_len);

    size_t line_bases = rec.line_bases;
    size_t line_bytes = rec.line_bytes;
    size_t seq_offset = rec.offset;

    size_t rowIndex = start / line_bases;
    size_t colIndex = start % line_bases;
    size_t filePos = seq_offset + rowIndex * line_bytes + colIndex;

    in.seekg(filePos, std::ios::beg);
    if (!in.good()) {
        throw std::runtime_error("Failed to seek to position " + std::to_string(filePos)
            + " in " + fasta_path_.string());
    }

    size_t to_read = req_len;
    while (to_read > 0) {
        char c;
        if (!in.get(c)) {
            throw std::runtime_error("Reached EOF unexpectedly while reading sub-sequence for " + seq_name);
        }
        if (c == '\n' || c == '\r') {
            continue;
        }
        result.push_back(c);
        to_read--;
    }

    in.close();
    return result;
}
