#include "SeqPro.h"
#include "data_process.h"
#include "hal/export.h"
#include "halAlignmentInstance.h"
#include "hdf5Alignment.h"
#include "hdf5ExternalArray.h"
#include "halBottomSegmentIterator.h"
#include "halColumnIterator.h"
#include "halGenome.h"
#include "halSequence.h"
#include "halTopSegmentIterator.h"
#include "ramaxHdf5BulkDna.h"
#include "ramesh.h"
#include "softmask_index.h"

#include <H5Cpp.h>
#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <omp.h>

#include <signal.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

using Managers = std::map<SpeciesName, SeqPro::SharedManagerVariant>;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        fail(message);
    }
}

std::set<std::string> readNamesWithStandardHdf5(
    const std::filesystem::path& path, const std::string& genome) {
    try {
        H5::H5File file(path.string(), H5F_ACC_RDONLY);
        auto dataset = file.openDataSet(genome + "/SEQNAME_ARRAY");
        const size_t width = dataset.getDataType().getSize();
        const size_t count =
            static_cast<size_t>(dataset.getSpace().getSimpleExtentNpoints());
        std::vector<char> bytes(width * count);
        dataset.read(bytes.data(), H5::StrType(H5::PredType::C_S1, width));
        std::set<std::string> names;
        for (size_t i = 0; i < count; ++i) {
            const char* begin = bytes.data() + i * width;
            const char* end =
                static_cast<const char*>(std::memchr(begin, '\0', width));
            require(end != nullptr, "HDF5 sequence name is not terminated");
            names.emplace(begin, static_cast<size_t>(end - begin));
        }
        return names;
    } catch (const H5::Exception& error) {
        fail("Standard HDF5 string read failed: " + error.getDetailMsg());
    }
}

class TemporaryDirectory {
public:
    TemporaryDirectory() {
        const std::string pattern =
            (std::filesystem::temp_directory_path() /
             "ramax-hal-memory-XXXXXX")
                .string();
        std::vector<char> writable(pattern.begin(), pattern.end());
        writable.push_back('\0');
        char* created = ::mkdtemp(writable.data());
        if (created == nullptr) {
            throw std::system_error(errno, std::generic_category(),
                                    "mkdtemp failed");
        }
        path_ = created;
    }

    ~TemporaryDirectory() {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }

    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

    const std::filesystem::path& path() const noexcept { return path_; }

private:
    std::filesystem::path path_;
};

class FailingChildNamesAlignment : public hal::Hdf5Alignment {
public:
    using hal::Hdf5Alignment::Hdf5Alignment;

    std::vector<std::string> getChildNames(
        const std::string& name) const override {
        if (fail_child_names) {
            throw std::bad_alloc();
        }
        return hal::Hdf5Alignment::getChildNames(name);
    }

    bool fail_child_names = false;
};

void testFailedGenomeClosePreservesCallerHandle(
    const std::filesystem::path& temp) {
    FailingChildNamesAlignment alignment(
        (temp / "close-failure.hal").string(),
        hal::CREATE_ACCESS,
        hal::hdf5DefaultFileCreatPropList(),
        hal::hdf5DefaultFileAccPropList(),
        hal::hdf5DefaultDSetCreatPropList());
    hal::Genome* root = alignment.addRootGenome("closeRoot", 0.0);
    hal::Genome* leaf =
        alignment.addLeafGenome("closeLeaf", "closeRoot", 0.1);
    alignment.fail_child_names = true;
    bool rejected = false;
    try {
        alignment.closeGenome(root);
    } catch (const std::bad_alloc&) {
        rejected = true;
    }
    alignment.fail_child_names = false;
    require(rejected, "genome close did not propagate the injected failure");
    require(root->getName() == "closeRoot",
            "failed genome close invalidated the caller's handle");
    alignment.closeGenome(root);
    alignment.closeGenome(leaf);
    alignment.close();
}

class ScopedFileSizeLimit {
public:
    explicit ScopedFileSizeLimit(rlim_t max_bytes) {
        if (::getrlimit(RLIMIT_FSIZE, &old_limit_) != 0) {
            throw std::system_error(errno, std::generic_category(),
                                    "getrlimit(RLIMIT_FSIZE) failed");
        }
        if (old_limit_.rlim_cur != RLIM_INFINITY &&
            old_limit_.rlim_cur <= max_bytes) {
            throw std::runtime_error(
                "fixture requires an initial file-size limit above its test cap");
        }

        struct sigaction ignored {};
        ignored.sa_handler = SIG_IGN;
        if (::sigemptyset(&ignored.sa_mask) != 0) {
            throw std::system_error(errno, std::generic_category(),
                                    "sigemptyset failed");
        }
        if (::sigaction(SIGXFSZ, &ignored, &old_action_) != 0) {
            throw std::system_error(errno, std::generic_category(),
                                    "sigaction(SIGXFSZ) failed");
        }
        signal_changed_ = true;

        struct rlimit limited = old_limit_;
        limited.rlim_cur = max_bytes;
        if (::setrlimit(RLIMIT_FSIZE, &limited) != 0) {
            const int saved_errno = errno;
            ::sigaction(SIGXFSZ, &old_action_, nullptr);
            signal_changed_ = false;
            throw std::system_error(saved_errno, std::generic_category(),
                                    "setrlimit(RLIMIT_FSIZE) failed");
        }
        limit_changed_ = true;
    }

    ~ScopedFileSizeLimit() { restoreNoThrow(); }

    ScopedFileSizeLimit(const ScopedFileSizeLimit&) = delete;
    ScopedFileSizeLimit& operator=(const ScopedFileSizeLimit&) = delete;

    void restore() {
        if (limit_changed_) {
            if (::setrlimit(RLIMIT_FSIZE, &old_limit_) != 0) {
                throw std::system_error(
                    errno, std::generic_category(),
                    "restoring RLIMIT_FSIZE failed");
            }
            limit_changed_ = false;
        }
        if (signal_changed_) {
            if (::sigaction(SIGXFSZ, &old_action_, nullptr) != 0) {
                throw std::system_error(
                    errno, std::generic_category(),
                    "restoring SIGXFSZ disposition failed");
            }
            signal_changed_ = false;
        }
    }

private:
    void restoreNoThrow() noexcept {
        if (limit_changed_ &&
            ::setrlimit(RLIMIT_FSIZE, &old_limit_) == 0) {
            limit_changed_ = false;
        }
        if (!limit_changed_ && signal_changed_ &&
            ::sigaction(SIGXFSZ, &old_action_, nullptr) == 0) {
            signal_changed_ = false;
        }
    }

    struct rlimit old_limit_ {};
    struct sigaction old_action_ {};
    bool limit_changed_ = false;
    bool signal_changed_ = false;
};

void writeFasta(
    const std::filesystem::path& path,
    const std::vector<std::pair<std::string, std::string>>& sequences) {
    std::ofstream output(path);
    for (const auto& [name, dna] : sequences) {
        output << '>' << name << '\n' << dna << '\n';
    }
    if (!output) {
        fail("failed to write fixture FASTA: " + path.string());
    }
}

void writeFileBytes(
    const std::filesystem::path& path,
    const std::string& contents) {
    std::ofstream output(path, std::ios::binary);
    output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    if (!output) {
        fail("failed to write fixture bytes: " + path.string());
    }
}

std::string readFileBytes(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        fail("failed to open fixture bytes: " + path.string());
    }
    std::string contents{
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()};
    if (input.bad()) {
        fail("failed to read fixture bytes: " + path.string());
    }
    return contents;
}

class ScopedOpenMpConfiguration {
public:
    explicit ScopedOpenMpConfiguration(int threads)
        : old_dynamic_(omp_get_dynamic()),
          old_max_threads_(omp_get_max_threads()) {
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
    }

    ~ScopedOpenMpConfiguration() {
        omp_set_num_threads(old_max_threads_);
        omp_set_dynamic(old_dynamic_);
    }

    ScopedOpenMpConfiguration(const ScopedOpenMpConfiguration&) = delete;
    ScopedOpenMpConfiguration& operator=(
        const ScopedOpenMpConfiguration&) = delete;

private:
    int old_dynamic_;
    int old_max_threads_;
};

std::vector<uint32_t> makeHdf5Payload(size_t size, uint32_t seed) {
    std::vector<uint32_t> payload(size);
    uint32_t state = seed;
    for (uint32_t& value : payload) {
        state = state * 1664525U + 1013904223U;
        value = state ^ (state >> 11U) ^ (state << 7U);
    }
    return payload;
}

std::vector<uint32_t> readHdf5Payload(
    H5::H5File& file, const std::string& dataset_name) {
    auto dataset = file.openDataSet(dataset_name);
    auto space = dataset.getSpace();
    require(space.getSimpleExtentNdims() == 1,
            "fixture HDF5 payload is not one-dimensional");
    hsize_t size = 0;
    space.getSimpleExtentDims(&size);
    std::vector<uint32_t> payload(static_cast<size_t>(size));
    dataset.read(payload.data(), H5::PredType::NATIVE_UINT32);
    return payload;
}

std::vector<uint32_t> readHdf5Payload(
    const std::filesystem::path& path, const std::string& dataset_name) {
    try {
        H5::H5File file(path.string(), H5F_ACC_RDONLY);
        return readHdf5Payload(file, dataset_name);
    } catch (const H5::Exception& error) {
        fail("Standard HDF5 payload read failed: " + error.getDetailMsg());
    }
}

void createFullBufferCompressedArray(
    const std::filesystem::path& path,
    const std::vector<uint32_t>& payload) {
    constexpr hsize_t chunk_elements = 64 * 1024;
    H5::H5File file(path.string(), H5F_ACC_TRUNC);
    H5::DSetCreatPropList properties;
    properties.setChunk(1, &chunk_elements);
    properties.setDeflate(6);

    hal::Hdf5ExternalArray array;
    array.create(&file, "payload", H5::PredType::NATIVE_UINT32,
                 static_cast<hsize_t>(payload.size()), &properties, 0);
    std::memcpy(array.getBuf(), payload.data(),
                payload.size() * sizeof(payload.front()));
    array.setDirty();
    array.write();
    require(readHdf5Payload(file, "payload") == payload,
            "same-file HDF5 read observed stale or malformed raw chunks");
}

void testParallelCompressedChunksAreReadableAndCoherent(
    const std::filesystem::path& temp) {
    constexpr size_t chunk_elements = 64 * 1024;
    auto expected =
        makeHdf5Payload(chunk_elements * 129 + 113, 0x42c0ffeeU);
    const auto single_path = temp / "compressed-single-thread.h5";
    const auto parallel_path = temp / "compressed-parallel.h5";

    try {
        {
            ScopedOpenMpConfiguration threads(1);
            createFullBufferCompressedArray(single_path, expected);
        }
        {
            ScopedOpenMpConfiguration threads(4);
            createFullBufferCompressedArray(parallel_path, expected);
        }

        const auto single = readHdf5Payload(single_path, "payload");
        const auto parallel = readHdf5Payload(parallel_path, "payload");
        require(single == expected && parallel == expected,
                "one-thread and multi-thread compressed arrays decode "
                "differently");

        const std::vector<size_t> rewrite_indexes{
            0, chunk_elements - 1, chunk_elements, chunk_elements + 1,
            expected.size() - 2, expected.size() - 1};
        for (size_t index : rewrite_indexes) {
            expected[index] ^= 0xa5a55a5aU;
        }
        {
            ScopedOpenMpConfiguration threads(4);
            H5::H5File file(parallel_path.string(), H5F_ACC_RDWR);
            hal::Hdf5ExternalArray array;
            array.load(&file, "payload", 0);
            for (size_t index : rewrite_indexes) {
                std::memcpy(array.getUpdate(static_cast<hsize_t>(index)),
                            &expected[index], sizeof(expected[index]));
            }
            array.write();
            require(readHdf5Payload(file, "payload") == expected,
                    "same-file read retained stale chunks after rewrite");
        }
        require(readHdf5Payload(parallel_path, "payload") == expected,
                "reopened HDF5 array lost a boundary rewrite");
    } catch (const H5::Exception& error) {
        fail("Parallel compressed-chunk fixture failed: " +
             error.getDetailMsg());
    }
}

void testCompressedArrayFallbacksPreservePayloads(
    const std::filesystem::path& temp) {
    constexpr hsize_t chunk_elements = 4096;
    const auto filtered =
        makeHdf5Payload(chunk_elements * 3 + 17, 0x12345678U);
    const auto paged =
        makeHdf5Payload(chunk_elements * 4 + 29, 0x87654321U);
    const auto filtered_path = temp / "compressed-shuffled.h5";
    const auto paged_path = temp / "compressed-paged.h5";

    try {
        {
            H5::H5File file(filtered_path.string(), H5F_ACC_TRUNC);
            H5::DSetCreatPropList properties;
            properties.setChunk(1, &chunk_elements);
            properties.setShuffle();
            properties.setDeflate(6);
            hal::Hdf5ExternalArray array;
            array.create(&file, "payload", H5::PredType::NATIVE_UINT32,
                         static_cast<hsize_t>(filtered.size()), &properties, 0);
            std::memcpy(array.getBuf(), filtered.data(),
                        filtered.size() * sizeof(filtered.front()));
            array.setDirty();
            array.write();
        }
        require(readHdf5Payload(filtered_path, "payload") == filtered,
                "multi-filter HDF5 fallback changed the payload");

        {
            H5::H5File file(paged_path.string(), H5F_ACC_TRUNC);
            H5::DSetCreatPropList properties;
            properties.setChunk(1, &chunk_elements);
            properties.setDeflate(6);
            hal::Hdf5ExternalArray array;
            array.create(&file, "payload", H5::PredType::NATIVE_UINT32,
                         static_cast<hsize_t>(paged.size()), &properties, 1);
            size_t offset = 0;
            while (offset < paged.size()) {
                const size_t count = std::min<size_t>(
                    static_cast<size_t>(chunk_elements),
                    paged.size() - offset);
                std::memcpy(
                    array.getUpdate(static_cast<hsize_t>(offset)),
                    paged.data() + offset, count * sizeof(paged.front()));
                offset += count;
            }
            array.write();
        }
        require(readHdf5Payload(paged_path, "payload") == paged,
                "partial-page compressed HDF5 fallback changed the payload");
    } catch (const H5::Exception& error) {
        fail("Compressed-array fallback fixture failed: " +
             error.getDetailMsg());
    }
}

struct FixtureInputs {
    Managers managers;
    SoftMask::PathMap softmask_paths;
};

void addLeaf(
    FixtureInputs& inputs,
    const std::filesystem::path& directory,
    const std::string& species,
    const std::vector<std::pair<std::string, std::string>>& sequences) {
    const auto source = directory / (species + ".input.fa");
    const auto uppercase = directory / (species + ".uppercase.fa");
    const auto index = directory / (species + ".softmask.bin");
    const auto marker = directory / (species + ".softmask.complete.json");
    writeFasta(source, sequences);
    SoftMask::ensureUppercaseFastaAndIndex(source, uppercase, index, marker);
    inputs.softmask_paths.emplace(species, index);
    SeqPro::ManagerVariant manager =
        std::make_unique<SeqPro::SequenceManager>(uppercase);
    inputs.managers.emplace(
        species,
        std::make_shared<SeqPro::ManagerVariant>(std::move(manager)));
}

void validateExternally(
    const std::string& hal_validate,
    const std::filesystem::path& hal_path) {
    const pid_t pid = ::fork();
    if (pid < 0) {
        throw std::system_error(errno, std::generic_category(), "fork failed");
    }
    if (pid == 0) {
        ::execlp(hal_validate.c_str(), hal_validate.c_str(),
                 hal_path.c_str(), static_cast<char*>(nullptr));
        _exit(errno == ENOENT ? 127 : 126);
    }
    int status = 0;
    while (::waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR) {
            throw std::system_error(errno, std::generic_category(),
                                    "waitpid failed");
        }
    }
    require(WIFEXITED(status) && WEXITSTATUS(status) == 0,
            "halValidate rejected " + hal_path.string());
}

RaMesh::BlockPtr makeBlock(
    const SpeciesName& reference_species,
    const ChrName& reference_chr,
    const std::vector<std::tuple<SpeciesName, ChrName, uint64_t, uint32_t,
                                 Strand, Cigar_t, RaMesh::AlignRole>>& rows) {
    auto block = RaMesh::Block::create(rows.size());
    block->ref_species = reference_species;
    block->ref_chr = reference_chr;
    for (const auto& [species, chr, start, length, strand, cigar, role] : rows) {
        block->anchors.emplace(
            RaMesh::SpeciesChrPair{species, chr},
            RaMesh::Segment::create(start, length, strand, cigar, role,
                                    RaMesh::SegmentRole::SEGMENT, block));
    }
    return block;
}

std::shared_ptr<const RaMesh::hal_export::PreparedExportInput>
prepareFixtureInput(
    const std::vector<std::weak_ptr<RaMesh::Block>>& blocks,
    const FixtureInputs& inputs,
    const std::filesystem::path& scratch_directory) {
    return RaMesh::hal_export::prepareExportInput(
        blocks, inputs.managers, scratch_directory);
}

struct MafRow {
    std::string source;
    uint64_t start = 0;
    uint64_t size = 0;
    char strand = '+';
    uint64_t source_size = 0;
    std::string dna;
};

std::vector<MafRow> readMafRows(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        fail("failed to open fixture MAF: " + path.string());
    }
    std::vector<MafRow> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.starts_with("s ")) {
            continue;
        }
        std::istringstream fields(line);
        char record_type = '\0';
        MafRow row;
        fields >> record_type >> row.source >> row.start >> row.size >>
            row.strand >> row.source_size >> row.dna;
        require(fields && record_type == 's',
                "malformed sequence row in fixture MAF");
        rows.push_back(std::move(row));
    }
    require(!input.bad(), "failed while reading fixture MAF");
    return rows;
}


void requireLeafContents(
    const hal::Genome* genome,
    const std::vector<std::pair<std::string, std::string>>& expected) {
    require(genome != nullptr, "HAL is missing an expected leaf genome");
    for (const auto& [name, dna] : expected) {
        const hal::Sequence* sequence = genome->getSequence(name);
        require(sequence != nullptr,
                "HAL leaf is missing original contig name " + name);
        require(sequence->getSequenceLength() == dna.size(),
                "HAL changed the length of contig " + name);
        std::string observed;
        sequence->getString(observed);
        require(observed == dna,
                "HAL changed sequence DNA or softmask for contig " + name);
    }
    require(genome->getNumSequences() == expected.size(),
            "HAL added or dropped a leaf contig");
}

void requireAllTopsUnmapped(const hal::Genome* genome) {
    require(genome != nullptr, "HAL is missing an expected leaf genome");
    require(genome->getNumTopSegments() != 0,
            "leaf genome has no top segments");
    for (hal_index_t index = 0;
         index < static_cast<hal_index_t>(genome->getNumTopSegments());
         ++index) {
        require(
            !genome->getTopSegmentIterator(index)->getTopSegment()->hasParent(),
            "leaf without homology acquired a parent mapping");
    }
}

uint64_t pairCoverage(
    const hal::Genome* reference,
    const hal::Genome* target) {
    std::set<const hal::Genome*> targets{target};
    auto columns = reference->getColumnIterator(
        &targets, 0, 0, hal::NULL_INDEX, false, true, false, true, false);
    uint64_t covered = 0;
    while (true) {
        const auto* column = columns->getColumnMap();
        bool present = false;
        for (const auto& [segment, copies] : *column) {
            if (segment->getGenome() == target && copies != nullptr &&
                !copies->empty()) {
                present = true;
                break;
            }
        }
        covered += static_cast<uint64_t>(present);
        if (columns->lastColumn()) {
            break;
        }
        columns->toRight();
    }
    return covered;
}

size_t requireReciprocalParentMappings(
    const hal::Genome* child,
    const hal::Genome* parent) {
    size_t mapped = 0;
    for (hal_index_t index = 0;
         index < static_cast<hal_index_t>(child->getNumTopSegments());
         ++index) {
        auto top = child->getTopSegmentIterator(index);
        if (!top->getTopSegment()->hasParent()) {
            continue;
        }
        ++mapped;
        const auto parent_index = top->getTopSegment()->getParentIndex();
        require(parent_index >= 0 &&
                    parent_index < static_cast<hal_index_t>(parent->getNumBottomSegments()),
                "top segment points outside its parent bottom array");
        auto bottom = parent->getBottomSegmentIterator(parent_index);
        require(bottom->getBottomSegment()->hasChildG(child),
                "mapped top segment has no reciprocal parent down edge");
        require(bottom->getBottomSegment()->getChildIndexG(child) == index,
                "parent down edge does not return to the child top segment");
        require(bottom->getLength() == top->getLength(),
                "parent and child homologous segments differ in length");
    }
    require(mapped != 0, "expected at least one mapped child top segment");
    return mapped;
}

void testSparseOccurrenceIdsAndRejectedTerminalEvidence() {
    using namespace RaMesh::hal_export;
    const uint64_t largest = std::numeric_limits<uint64_t>::max();
    const std::vector<uint64_t> ids{largest, 17, 0};
    const std::unordered_map<uint64_t, RunOrderKey> no_keys;
    require(buildMaximumCardinalityWeightPathCover(ids, {}, no_keys) ==
                std::vector<std::vector<uint64_t>>{{0}, {17}, {largest}},
            "sparse occurrence IDs lost numeric fallback ordering");
    require(buildMaximumCardinalityWeightPathCover({}, {}, no_keys).empty(),
            "empty occurrence domain produced a path");

    bool duplicate_rejected = false;
    try {
        buildMaximumCardinalityWeightPathCover({17, 0, 17}, {}, no_keys);
    } catch (const std::runtime_error&) {
        duplicate_rejected = true;
    }
    require(duplicate_rejected, "duplicate occurrence IDs were accepted");

    TerminalEndSupport unknown;
    unknown.end = {42, OccurrenceEndSide::RIGHT};
    unknown.supporting_lineages = {1, 2};
    bool terminal_rejected = false;
    try {
        buildAncestralSequenceAssembly(ids, {}, no_keys, {unknown}, 1);
    } catch (const std::invalid_argument&) {
        terminal_rejected = true;
    }
    require(terminal_rejected,
            "terminal evidence outside the sparse occurrence domain was accepted");
}

void testReorientedMaskedAmbiguityKeepsHistoricalConsensus(
    const std::filesystem::path& temp,
    const std::string& hal_validate) {
    FixtureInputs inputs;
    addLeaf(inputs, temp, "caseAnchor", {{"one", "N"}});
    addLeaf(inputs, temp, "caseWitness", {{"one", "N"}});
    addLeaf(inputs, temp, "caseCarrier", {{"one", "n"}});
    addLeaf(inputs, temp, "caseGuide", {{"one", "N"}});
    const Cigar_t match1{cigarToInt('M', 1)};
    auto first = makeBlock(
        "caseAnchor", "one",
        {{"caseAnchor", "one", 0, 1, Strand::FORWARD, match1,
          RaMesh::AlignRole::PRIMARY},
         {"caseWitness", "one", 0, 1, Strand::FORWARD, match1,
          RaMesh::AlignRole::PRIMARY}});
    auto overlapping = makeBlock(
        "caseGuide", "one",
        {{"caseAnchor", "one", 0, 1, Strand::REVERSE, match1,
          RaMesh::AlignRole::PRIMARY},
         {"caseCarrier", "one", 0, 1, Strand::REVERSE, match1,
          RaMesh::AlignRole::PRIMARY},
         {"caseGuide", "one", 0, 1, Strand::FORWARD, match1,
          RaMesh::AlignRole::PRIMARY}});
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks{first, overlapping};
    const auto prepared = prepareFixtureInput(
        blocks, inputs, temp / "prepared-reoriented-mask");
    RaMesh::hal_export::exportToMaf(
        *prepared, temp / "reoriented-mask.maf", inputs.managers, false);
    const auto output = temp / "reoriented-mask.hal";
    RaMesh::hal_export::exportToHal(
        *prepared, output, inputs.managers,
        NewickParser("((caseAnchor:1,caseWitness:1)caseOther:1,"
                     "(caseCarrier:0.001,caseGuide:1)casePair:1)caseRoot;"),
        "caseRoot", SoftMask::loadIndexes(inputs.softmask_paths));
    validateExternally(hal_validate, output);
    hal::AlignmentPtr alignment =
        hal::openHalAlignment(output.string(), nullptr, hal::READ_ACCESS);
    const hal::Genome* parent = alignment->openGenome("casePair");
    const hal::Genome* carrier = alignment->openGenome("caseCarrier");
    require(parent && carrier, "reoriented masked fixture lost HAL genomes");
    std::string ancestor;
    parent->getString(ancestor);
    // The existing alignment complement maps n to N. Moving masking past
    // normalization would change this consensus and can change duplicate checks.
    require(ancestor == "N",
            "shared preparation changed masked ambiguity normalization");
    requireLeafContents(carrier, {{"one", "n"}});
    requireReciprocalParentMappings(carrier, parent);
    alignment->closeGenome(carrier);
    alignment->closeGenome(parent);
}

void testMixedCaseReverseIndelAndMissingChild(
    const std::filesystem::path& temp,
    const std::string& hal_validate) {
    FixtureInputs inputs;
    const std::vector<std::pair<std::string, std::string>> leaf_a{
        {"alpha.long", "AaCCggTTACG"}, {"A_short", "tTaC"}};
    const std::vector<std::pair<std::string, std::string>> leaf_b{
        {"bUnequal", "AAttGGttCC"}};
    const std::vector<std::pair<std::string, std::string>> leaf_c{
        {"gamma_name", "aacCGGTTa"}};
    const std::vector<std::pair<std::string, std::string>> leaf_d{
        {"unrepresented", "TtAaCcG"}};
    addLeaf(inputs, temp, "leafA", leaf_a);
    addLeaf(inputs, temp, "leafB", leaf_b);
    addLeaf(inputs, temp, "leafC", leaf_c);
    addLeaf(inputs, temp, "leafD", leaf_d);

    const Cigar_t match8{cigarToInt('M', 8)};
    const Cigar_t gapped{
        cigarToInt('M', 2), cigarToInt('D', 2),
        cigarToInt('I', 2), cigarToInt('M', 4)};
    auto block = makeBlock(
        "leafA", "alpha.long",
        {{"leafA", "alpha.long", 0, 8, Strand::FORWARD, match8,
          RaMesh::AlignRole::PRIMARY},
         {"leafB", "bUnequal", 0, 8, Strand::REVERSE, gapped,
          RaMesh::AlignRole::PRIMARY},
         {"leafC", "gamma_name", 0, 8, Strand::FORWARD, match8,
          RaMesh::AlignRole::PRIMARY}});

    const auto hal_path = temp / "mixed-reverse-missing.hal";
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks{block};
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-mixed");
    RaMesh::hal_export::exportToHal(
        *prepared, hal_path, inputs.managers,
        NewickParser("((leafA:0.1,leafB:0.1)ancAB:0.1,"
                     "(leafC:0.1,leafD:0.1)ancCD:0.1)root;"),
        "root", SoftMask::loadIndexes(inputs.softmask_paths));
    validateExternally(hal_validate, hal_path);
    require(readNamesWithStandardHdf5(hal_path, "leafA") ==
                std::set<std::string>{"alpha.long", "A_short"},
            "standard HDF5 readers must recover complete sequence names");

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(hal_path.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment), "cannot reopen mixed HAL");
    const hal::Genome* a = alignment->openGenome("leafA");
    const hal::Genome* b = alignment->openGenome("leafB");
    const hal::Genome* c = alignment->openGenome("leafC");
    const hal::Genome* d = alignment->openGenome("leafD");
    const hal::Genome* anc_ab = alignment->openGenome("ancAB");
    const hal::Genome* anc_cd = alignment->openGenome("ancCD");
    const hal::Genome* root = alignment->openGenome("root");
    require(a && b && c && d && anc_ab && anc_cd && root,
            "mixed HAL is missing a tree genome");

    requireLeafContents(a, leaf_a);
    requireLeafContents(b, leaf_b);
    requireLeafContents(c, leaf_c);
    requireLeafContents(d, leaf_d);
    require(pairCoverage(a, b) == 6 && pairCoverage(b, a) == 6,
            "reverse indel changed pairwise homology coverage");
    require(pairCoverage(a, c) == 8,
            "aligned sibling subtree lost elementary runs");
    require(pairCoverage(d, a) == 0,
            "unrepresented child incorrectly acquired homology");

    size_t reverse_mapped = 0;
    for (hal_index_t index = 0;
         index < static_cast<hal_index_t>(b->getNumTopSegments());
         ++index) {
        auto top = b->getTopSegmentIterator(index);
        if (top->getTopSegment()->hasParent() &&
            top->getTopSegment()->getParentReversed()) {
            ++reverse_mapped;
        }
    }
    require(reverse_mapped != 0,
            "reverse occurrence lost its HAL parent orientation");
    requireReciprocalParentMappings(a, anc_ab);
    requireReciprocalParentMappings(b, anc_ab);
    requireReciprocalParentMappings(c, anc_cd);
    requireReciprocalParentMappings(anc_ab, root);
    requireReciprocalParentMappings(anc_cd, root);

    for (hal_index_t index = 0;
         index < static_cast<hal_index_t>(d->getNumTopSegments());
         ++index) {
        require(!d->getTopSegmentIterator(index)->getTopSegment()->hasParent(),
                "missing child run was emitted as aligned");
    }
    bool observed_sparse_child = false;
    for (hal_index_t index = 0;
         index < static_cast<hal_index_t>(anc_cd->getNumBottomSegments());
         ++index) {
        auto bottom = anc_cd->getBottomSegmentIterator(index);
        if (bottom->getBottomSegment()->hasChildG(c) &&
            !bottom->getBottomSegment()->hasChildG(d)) {
            observed_sparse_child = true;
        }
    }
    require(observed_sparse_child,
            "ancestor did not retain a run absent from one child");

    alignment->closeGenome(d);
    alignment->closeGenome(c);
    alignment->closeGenome(b);
    alignment->closeGenome(a);
    alignment->closeGenome(anc_cd);
    alignment->closeGenome(anc_ab);
    alignment->closeGenome(root);
    alignment.reset();
}

void testPhysicalCopiesUseParalogyEdges(
    const std::filesystem::path& temp,
    const std::string& hal_validate) {
    FixtureInputs inputs;
    addLeaf(inputs, temp, "copyA", {{"one", "AaCcGGTT"}});
    addLeaf(inputs, temp, "copyB", {{"copies", "AaCcAaCc"}});
    const Cigar_t match4{cigarToInt('M', 4)};
    auto block = makeBlock(
        "copyA", "one",
        {{"copyA", "one", 0, 4, Strand::FORWARD, match4,
          RaMesh::AlignRole::PRIMARY},
         {"copyB", "copies", 0, 4, Strand::FORWARD, match4,
          RaMesh::AlignRole::PRIMARY},
         {"copyB", "copies", 4, 4, Strand::FORWARD, match4,
          RaMesh::AlignRole::PRIMARY}});
    const auto hal_path = temp / "physical-copies.hal";
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks{block};
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-copies");
    RaMesh::hal_export::exportToHal(
        *prepared, hal_path, inputs.managers,
        NewickParser("(copyA:0.1,copyB:0.1)copyRoot;"), "copyRoot",
        SoftMask::loadIndexes(inputs.softmask_paths));
    validateExternally(hal_validate, hal_path);

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(hal_path.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment), "cannot reopen paralogy HAL");
    const hal::Genome* root = alignment->openGenome("copyRoot");
    const hal::Genome* copy_b = alignment->openGenome("copyB");
    require(root && copy_b, "paralogy HAL is missing a genome");
    require(root->getNumBottomSegments() == 1,
            "child-only copy created or removed an ancestral occurrence");
    std::vector<hal_index_t> aligned;
    for (hal_index_t index = 0;
         index < static_cast<hal_index_t>(copy_b->getNumTopSegments());
         ++index) {
        auto top = copy_b->getTopSegmentIterator(index);
        if (top->getTopSegment()->hasParent()) {
            aligned.push_back(index);
        }
    }
    require(aligned.size() == 2,
            "same-run physical copies were collapsed or left unmapped");
    auto first = copy_b->getTopSegmentIterator(aligned[0]);
    auto second = copy_b->getTopSegmentIterator(aligned[1]);
    require(first->getTopSegment()->getParentIndex() ==
                second->getTopSegment()->getParentIndex(),
            "paralogous child copies do not share an ancestral segment");
    require(first->getTopSegment()->hasNextParalogy() &&
                second->getTopSegment()->hasNextParalogy() &&
                first->getTopSegment()->getNextParalogyIndex() == aligned[1] &&
                second->getTopSegment()->getNextParalogyIndex() == aligned[0],
            "two physical copies do not form a complete paralogy cycle");
    require(first->getTopSegment()->isCanonicalParalog() !=
                second->getTopSegment()->isCanonicalParalog(),
            "paralogy cycle does not have exactly one canonical copy");
    auto parent = root->getBottomSegmentIterator(
        first->getTopSegment()->getParentIndex());
    const hal_index_t canonical = first->getTopSegment()->isCanonicalParalog() ? aligned[0] : aligned[1];
    require(parent->getBottomSegment()->getChildIndexG(copy_b) == canonical,
            "ancestral down edge does not select the canonical paralog");
    alignment->closeGenome(copy_b);
    alignment->closeGenome(root);
    alignment.reset();
}

void testParentProjectionPreservesChildContainers(
    const std::filesystem::path& temp,
    const std::string& hal_validate) {
    FixtureInputs inputs;
    for (const std::string species : {"projA", "projB"}) {
        addLeaf(inputs, temp, species,
                {{"left_piece", "AaAa"}, {"right_piece", "cCcC"}});
    }
    for (const std::string species : {"projC", "projD"}) {
        addLeaf(inputs, temp, species, {{"joined", "cCcCAaAa"}});
    }
    const Cigar_t match4{cigarToInt('M', 4)};
    auto make_projection_block = [&](const std::string& ab_chr,
                                     uint64_t cd_start) {
        return makeBlock(
            "projA", ab_chr,
            {{"projA", ab_chr, 0, 4, Strand::FORWARD, match4,
              RaMesh::AlignRole::PRIMARY},
             {"projB", ab_chr, 0, 4, Strand::FORWARD, match4,
              RaMesh::AlignRole::PRIMARY},
             {"projC", "joined", cd_start, 4, Strand::FORWARD, match4,
              RaMesh::AlignRole::PRIMARY},
             {"projD", "joined", cd_start, 4, Strand::FORWARD, match4,
              RaMesh::AlignRole::PRIMARY}});
    };
    auto first = make_projection_block("left_piece", 4);
    auto second = make_projection_block("right_piece", 0);
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks{first, second};
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-projection");
    const auto hal_path = temp / "parent-projection.hal";
    RaMesh::hal_export::exportToHal(
        *prepared, hal_path, inputs.managers,
        NewickParser("((projA:0.1,projB:0.1)projAB:0.1,"
                     "(projC:0.1,projD:0.1)projCD:0.1)projRoot;"),
        "projRoot", SoftMask::loadIndexes(inputs.softmask_paths));
    validateExternally(hal_validate, hal_path);

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(hal_path.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment), "cannot reopen projection HAL");
    const hal::Genome* root = alignment->openGenome("projRoot");
    const hal::Genome* ab = alignment->openGenome("projAB");
    const hal::Genome* cd = alignment->openGenome("projCD");
    require(root && ab && cd, "projection HAL is missing an ancestor");
    require(root->getNumSequences() == 1 && ab->getNumSequences() == 1 &&
                cd->getNumSequences() == 1,
            "parent projection retained incompatible temporary child containers");
    require(root->getNumBottomSegments() == 2 &&
                ab->getNumTopSegments() == 3 &&
                ab->getNumBottomSegments() == 3 &&
                cd->getNumTopSegments() == 2 &&
                cd->getNumBottomSegments() == 2,
            "top-down projection no longer tiles parent and child segments");
    requireReciprocalParentMappings(ab, root);
    requireReciprocalParentMappings(cd, root);
    alignment->closeGenome(cd);
    alignment->closeGenome(ab);
    alignment->closeGenome(root);
    alignment.reset();
}

std::string buildBalancedTree(
    const std::vector<std::string>& leaves,
    size_t begin,
    size_t end,
    size_t& internal_index,
    bool root) {
    if (end - begin == 1) {
        return leaves[begin] + ":0.1";
    }
    const size_t middle = begin + (end - begin) / 2;
    const std::string left =
        buildBalancedTree(leaves, begin, middle, internal_index, false);
    const std::string right =
        buildBalancedTree(leaves, middle, end, internal_index, false);
    const std::string name = root
        ? "wideRoot"
        : "wideAnc" + std::to_string(internal_index++);
    return "(" + left + "," + right + ")" + name +
           (root ? ";" : ":0.1");
}

void testMoreThanSixtyFourTreeNodes(
    const std::filesystem::path& temp,
    const std::string& hal_validate) {
    constexpr size_t leaf_count = 65;
    FixtureInputs inputs;
    std::vector<std::string> leaves;
    leaves.reserve(leaf_count);
    std::vector<std::tuple<SpeciesName, ChrName, uint64_t, uint32_t, Strand,
                           Cigar_t, RaMesh::AlignRole>> rows;
    rows.reserve(leaf_count);
    const Cigar_t match1{cigarToInt('M', 1)};
    for (size_t index = 0; index < leaf_count; ++index) {
        std::string leaf = "wide";
        if (index < 10) leaf += '0';
        leaf += std::to_string(index);
        leaves.push_back(leaf);
        addLeaf(inputs, temp, leaf, {{"tiny", index % 2 == 0 ? "a" : "A"}});
        rows.emplace_back(leaf, "tiny", 0, 1, Strand::FORWARD, match1,
                          RaMesh::AlignRole::PRIMARY);
    }
    size_t internal_index = 0;
    const std::string tree =
        buildBalancedTree(leaves, 0, leaves.size(), internal_index, true);
    auto block = makeBlock(leaves.front(), "tiny", rows);
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks{block};
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-wide");
    const auto hal_path = temp / "wide-tree.hal";
    RaMesh::hal_export::exportToHal(
        *prepared, hal_path, inputs.managers, NewickParser(tree), "wideRoot",
        SoftMask::loadIndexes(inputs.softmask_paths));
    validateExternally(hal_validate, hal_path);

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(hal_path.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment), "cannot reopen wide-tree HAL");
    require(alignment->getNumGenomes() == leaf_count * 2 - 1,
            ">64-node tree was truncated during HAL export");
    const hal::Genome* first = alignment->openGenome(leaves.front());
    const hal::Genome* last = alignment->openGenome(leaves.back());
    require(first && last && pairCoverage(first, last) == 1 &&
                pairCoverage(last, first) == 1,
            "leaf occurrence above bit 63 disappeared from the alignment");
    requireLeafContents(first, {{"tiny", "a"}});
    requireLeafContents(last, {{"tiny", "a"}});
    alignment->closeGenome(last);
    alignment->closeGenome(first);
    alignment.reset();
}

void testNoHomologyPreservesLeafDNAAndUnmappedTops(
    const std::filesystem::path& temp) {
    FixtureInputs inputs;
    const std::vector<std::pair<std::string, std::string>> orphan_a{
        {"left.odd", "aCgTtAA"}, {"left.second", "GGttC"}};
    const std::vector<std::pair<std::string, std::string>> orphan_b{
        {"right_only", "TtaACCGg"}};
    addLeaf(inputs, temp, "orphanA", orphan_a);
    addLeaf(inputs, temp, "orphanB", orphan_b);

    const auto hal_path = temp / "no-homology.hal";
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks;
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-no-homology");
    RaMesh::hal_export::exportToHal(
        *prepared, hal_path, inputs.managers,
        NewickParser("(orphanA:0.1,orphanB:0.1)orphanRoot;"),
        "orphanRoot", SoftMask::loadIndexes(inputs.softmask_paths));

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(hal_path.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment), "cannot reopen no-homology HAL");
    const hal::Genome* root = alignment->openGenome("orphanRoot");
    const hal::Genome* a = alignment->openGenome("orphanA");
    const hal::Genome* b = alignment->openGenome("orphanB");
    require(root && a && b, "no-homology HAL is missing a tree genome");
    require(root->getNumSequences() == 0 &&
                root->getNumTopSegments() == 0 &&
                root->getNumBottomSegments() == 0,
            "zero-sequence ancestor unexpectedly contains HAL segments");
    requireLeafContents(a, orphan_a);
    requireLeafContents(b, orphan_b);
    requireAllTopsUnmapped(a);
    requireAllTopsUnmapped(b);
    alignment->closeGenome(b);
    alignment->closeGenome(a);
    alignment->closeGenome(root);
    alignment.reset();
}

void testBulkDnaEncodingAndPackedNibbleBoundaries(
    const std::filesystem::path& temp) {
    std::string original(64, 'A');
    constexpr char cycle[] = "TgCa";
    for (size_t index = 0; index < original.size(); ++index) {
        original[index] = cycle[index % 4];
    }

    FixtureInputs inputs;
    addLeaf(inputs, temp, "packedLeaf", {{"packed", original}});
    addLeaf(inputs, temp, "packedPeer", {{"peer", "AcGt"}});
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks;
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-packed-dna");
    const auto hal_path = temp / "packed-dna.hal";
    RaMesh::hal_export::exportToHal(
        *prepared, hal_path, inputs.managers,
        NewickParser("(packedLeaf:0.1,packedPeer:0.1)packedRoot;"),
        "packedRoot", SoftMask::loadIndexes(inputs.softmask_paths));

    const std::string ambiguous_lower = "kmryuswbdhv";
    const std::string ambiguous_upper = "KMRYUSWBDHV";
    const std::string update =
        "AaCcGgTtNn" + ambiguous_lower + ambiguous_upper;
    const std::string normalized =
        "AaCcGgTtNn" + std::string(ambiguous_lower.size(), 'n') +
        std::string(ambiguous_upper.size(), 'N');
    std::string expected = original;
    expected.replace(1, update.size(), normalized);

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(
            hal_path.string(), nullptr, hal::WRITE_ACCESS);
    require(static_cast<bool>(alignment),
            "cannot reopen packed-DNA fixture for update");
    hal::Genome* leaf = alignment->openGenome("packedLeaf");
    require(leaf != nullptr, "packed-DNA fixture lost its writable leaf");
    hal::ramaxWriteBulkDna(
        *leaf, 1, update.data(), static_cast<hal_size_t>(update.size()));
    requireLeafContents(leaf, {{"packed", expected}});

    const char invalid = static_cast<char>(0xffU);
    bool rejected = false;
    try {
        hal::ramaxWriteBulkDna(*leaf, 1, &invalid, 1);
    } catch (const std::exception&) {
        rejected = true;
    }
    require(rejected, "packed-DNA lookup accepted an invalid byte");
    requireLeafContents(leaf, {{"packed", expected}});
    alignment->closeGenome(leaf);
    alignment.reset();

    alignment =
        hal::openHalAlignment(
            hal_path.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment),
            "cannot reopen packed-DNA fixture after update");
    const hal::Genome* reopened = alignment->openGenome("packedLeaf");
    requireLeafContents(reopened, {{"packed", expected}});
    alignment->closeGenome(reopened);
    alignment.reset();
}

void testLeafDNAStreamsAcrossOneMiBBoundary(
    const std::filesystem::path& temp) {
    constexpr size_t boundary = 1024 * 1024;
    constexpr size_t sequence_length = boundary + 257;
    std::string long_dna(sequence_length, 'A');
    constexpr char bases[] = "ACGT";
    for (size_t index = 0; index < long_dna.size(); ++index) {
        long_dna[index] = bases[index % 4];
    }
    const auto lowercase = [&](size_t begin, size_t end) {
        for (size_t index = begin; index < end; ++index) {
            long_dna[index] = static_cast<char>(long_dna[index] - 'A' + 'a');
        }
    };
    lowercase(boundary - 7, boundary + 12);
    lowercase(boundary + 137, boundary + 154);
    lowercase(boundary + 201, boundary + 204);
    long_dna[boundary - 1] = 'N';
    long_dna[boundary] = 'n';
    long_dna.back() = 't';

    FixtureInputs inputs;
    const std::vector<std::pair<std::string, std::string>> streamed{
        {"a.odd-prefix", "nAc"},
        {"boundary.contig", long_dna},
        {"z.odd-suffix", "Tgn"}};
    const std::vector<std::pair<std::string, std::string>> peer{
        {"peer.contig", "AcGtA"}};
    addLeaf(inputs, temp, "streamLeaf", streamed);
    addLeaf(inputs, temp, "streamPeer", peer);

    const auto hal_path = temp / "stream-boundary.hal";
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks;
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-stream");
    RaMesh::hal_export::exportToHal(
        *prepared, hal_path, inputs.managers,
        NewickParser("(streamLeaf:0.1,streamPeer:0.1)streamRoot;"),
        "streamRoot", SoftMask::loadIndexes(inputs.softmask_paths));

    hal::AlignmentPtr alignment =
        hal::openHalAlignment(hal_path.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment), "cannot reopen streaming HAL");
    const hal::Genome* leaf = alignment->openGenome("streamLeaf");
    requireLeafContents(leaf, streamed);
    alignment->closeGenome(leaf);
    const hal::Genome* peer_genome = alignment->openGenome("streamPeer");
    requireLeafContents(peer_genome, peer);
    alignment->closeGenome(peer_genome);
    alignment.reset();
}

void testPreparedInputSurvivesGraphReleaseAndReplays(
    const std::filesystem::path& temp,
    const std::string& hal_validate) {
    FixtureInputs inputs;
    const std::vector<std::pair<std::string, std::string>> lifetime_a{
        {"chr.one", "AaCCggTT"}};
    const std::vector<std::pair<std::string, std::string>> lifetime_b{
        {"chr.two", "aAccGGtt"}};
    addLeaf(inputs, temp, "lifetimeA", lifetime_a);
    addLeaf(inputs, temp, "lifetimeB", lifetime_b);

    RaMesh::RaMeshMultiGenomeGraph graph(inputs.managers);
    Anchor anchor(
        0, 0, 8, 0, 0, 8, Strand::FORWARD, 8, 8,
        Cigar_t{cigarToInt('M', 8)});
    graph.insertAnchorIntoGraph(
        *inputs.managers.at("lifetimeA"),
        *inputs.managers.at("lifetimeB"),
        "lifetimeA", "lifetimeB", anchor);
    require(graph.blocks.size() == 1,
            "lifetime fixture did not build exactly one graph block");

    std::weak_ptr<RaMesh::Block> source_block = graph.blocks.front();
    std::weak_ptr<RaMesh::Segment> source_segment;
    {
        const auto block = source_block.lock();
        require(block && !block->anchors.empty(),
                "lifetime fixture graph block has no segment");
        source_segment = block->anchors.begin()->second;
    }
    const auto prepared = prepareFixtureInput(
        graph.blocks, inputs, temp / "prepared-lifetime");
    const auto maf_first_prepared = prepareFixtureInput(
        graph.blocks, inputs, temp / "prepared-lifetime-maf-first");
    graph.clearAllGraphs();
    require(source_block.expired() && source_segment.expired(),
            "cleared graph retained source alignment objects");
    const auto& original_manager =
        *std::get<std::unique_ptr<SeqPro::SequenceManager>>(
            *inputs.managers.at("lifetimeA"));
    const auto sequence_id = original_manager.getSequenceId("chr.one");
    std::string_view borrowed_dna;
    require(original_manager.tryGetContiguousSubSequence(
                sequence_id, 0, 8, borrowed_dna),
            "lifetime fixture DNA is not physically contiguous");
    original_manager.releaseMappedPages();
    require(borrowed_dna == "AACCGGTT",
            "releasing mapped pages invalidated a borrowed DNA view");
    require(original_manager.getSubSequence(sequence_id, 1, 6) == "ACCGGT",
            "source DNA cannot be read again after releasing mapped pages");

    const auto indexes = SoftMask::loadIndexes(inputs.softmask_paths);
    const std::string tree =
        "(lifetimeA:0.1,lifetimeB:0.1)lifetimeRoot;";
    const auto first_hal = temp / "prepared-lifetime-first.hal";
    RaMesh::hal_export::exportToHal(
        *prepared, first_hal, inputs.managers, NewickParser(tree),
        "lifetimeRoot", indexes);
    validateExternally(hal_validate, first_hal);

    const auto first_maf = temp / "prepared-lifetime-first.maf";
    const auto second_maf = temp / "prepared-lifetime-second.maf";
    RaMesh::hal_export::exportToMaf(
        *prepared, first_maf, inputs.managers, false);
    RaMesh::hal_export::exportToMaf(
        *prepared, second_maf, inputs.managers, false);
    require(readFileBytes(first_maf) == readFileBytes(second_maf),
            "replaying prepared input changed consumer-visible MAF output");

    const auto maf_rows = readMafRows(first_maf);
    require(maf_rows.size() == 2,
            "prepared-input MAF did not retain both aligned leaves");
    std::map<std::string, MafRow> maf_by_source;
    for (const auto& row : maf_rows) {
        maf_by_source.emplace(row.source, row);
    }
    for (const std::string source :
         {"lifetimeA.chr.one", "lifetimeB.chr.two"}) {
        const auto row = maf_by_source.find(source);
        require(row != maf_by_source.end() &&
                    row->second.start == 0 &&
                    row->second.size == 8 &&
                    row->second.strand == '+' &&
                    row->second.source_size == 8 &&
                    row->second.dna == "AACCGGTT",
                "prepared-input MAF changed aligned DNA or coordinates for " +
                    source);
    }

    const auto second_hal = temp / "prepared-lifetime-second.hal";
    RaMesh::hal_export::exportToHal(
        *prepared, second_hal, inputs.managers, NewickParser(tree),
        "lifetimeRoot", indexes);
    hal::AlignmentPtr alignment =
        hal::openHalAlignment(second_hal.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment),
            "cannot reopen HAL from replayed prepared input");
    const hal::Genome* root = alignment->openGenome("lifetimeRoot");
    const hal::Genome* a = alignment->openGenome("lifetimeA");
    const hal::Genome* b = alignment->openGenome("lifetimeB");
    require(root && a && b,
            "replayed prepared-input HAL is missing a tree genome");
    requireLeafContents(a, lifetime_a);
    requireLeafContents(b, lifetime_b);
    require(pairCoverage(a, b) == 8 && pairCoverage(b, a) == 8,
            "replayed prepared-input HAL changed leaf relationships");
    requireReciprocalParentMappings(a, root);
    requireReciprocalParentMappings(b, root);
    alignment->closeGenome(b);
    alignment->closeGenome(a);
    alignment->closeGenome(root);
    alignment.reset();

    const auto maf_first_path = temp / "prepared-maf-first.maf";
    RaMesh::hal_export::exportToMaf(
        *maf_first_prepared, maf_first_path, inputs.managers, false);
    require(readFileBytes(maf_first_path) == readFileBytes(first_maf),
            "HAL-first preparation contaminated MAF case or coordinates");
    const auto hal_after_maf = temp / "prepared-hal-after-maf.hal";
    RaMesh::hal_export::exportToHal(
        *maf_first_prepared, hal_after_maf, inputs.managers,
        NewickParser(tree), "lifetimeRoot", indexes);
    alignment =
        hal::openHalAlignment(hal_after_maf.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment),
            "cannot reopen HAL after MAF-first preparation");
    root = alignment->openGenome("lifetimeRoot");
    a = alignment->openGenome("lifetimeA");
    b = alignment->openGenome("lifetimeB");
    require(root && a && b, "MAF-first preparation lost HAL genomes");
    requireLeafContents(a, lifetime_a);
    requireLeafContents(b, lifetime_b);
    require(pairCoverage(a, b) == 8 && pairCoverage(b, a) == 8,
            "MAF-first preparation changed HAL leaf relationships");
    requireReciprocalParentMappings(a, root);
    requireReciprocalParentMappings(b, root);
    alignment->closeGenome(b);
    alignment->closeGenome(a);
    alignment->closeGenome(root);
    alignment.reset();
}

void testMafWriteFailureIsTransactionalAndRecoverable(
    const std::filesystem::path& temp) {
    constexpr size_t sequence_length = 128 * 1024 + 3;
    std::string dna(sequence_length, 'A');
    constexpr char bases[] = "ACGT";
    for (size_t index = 0; index < dna.size(); ++index) {
        dna[index] = bases[index % 4];
    }

    FixtureInputs inputs;
    addLeaf(inputs, temp, "mafLimit", {{"long.contig", dna}});
    addLeaf(inputs, temp, "mafPeer", {{"long.contig", dna}});
    auto block = makeBlock(
        "mafLimit", "long.contig",
        {{"mafLimit", "long.contig", 0,
          static_cast<uint32_t>(sequence_length), Strand::FORWARD,
          Cigar_t{cigarToInt('M', static_cast<uint32_t>(sequence_length))},
          RaMesh::AlignRole::PRIMARY},
         {"mafPeer", "long.contig", 0,
          static_cast<uint32_t>(sequence_length), Strand::FORWARD,
          Cigar_t{cigarToInt('M', static_cast<uint32_t>(sequence_length))},
          RaMesh::AlignRole::PRIMARY}});
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks{block};
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-maf-failure");
    const auto failure_prepared = prepareFixtureInput(
        blocks, inputs, temp / "prepared-maf-cold-failure");

    const auto baseline_path = temp / "maf-baseline.maf";
    RaMesh::hal_export::exportToMaf(
        *prepared, baseline_path, inputs.managers, false);
    require(std::filesystem::file_size(baseline_path) > 4096,
            "unlimited baseline MAF does not exceed the file-size cap");

    const auto destination = temp / "maf-transaction.maf";
    constexpr char original_data[] =
        "pre-existing MAF destination\0with binary bytes";
    const std::string original_bytes(
        original_data, sizeof(original_data) - 1);
    writeFileBytes(destination, original_bytes);
    std::filesystem::path temporary = destination;
    temporary += ".tmp";

    bool caught_std_exception = false;
    bool caught_nonstd_exception = false;
    {
        ScopedFileSizeLimit file_size_limit(4096);
        try {
            RaMesh::hal_export::exportToMaf(
                *failure_prepared, destination, inputs.managers, false);
        } catch (const std::exception&) {
            caught_std_exception = true;
        } catch (...) {
            caught_nonstd_exception = true;
        }
        file_size_limit.restore();
    }
    require(!caught_nonstd_exception,
            "limited MAF preparation escaped outside std::exception");
    require(caught_std_exception,
            "limited MAF preparation or write unexpectedly completed");
    require(readFileBytes(destination) == original_bytes,
            "failed MAF export changed the pre-existing destination");
    require(!std::filesystem::exists(temporary),
            "failed MAF export retained its temporary output");

    RaMesh::hal_export::exportToMaf(
        *failure_prepared, destination, inputs.managers, false);
    const auto rows = readMafRows(destination);
    require(rows.size() == 2, "recovered MAF omitted an aligned leaf");
    for (const std::string source : {"mafLimit.long.contig", "mafPeer.long.contig"}) {
        const auto row = std::find_if(rows.begin(), rows.end(),
            [&source](const MafRow& value) { return value.source == source; });
        require(row != rows.end() &&
                    row->start == 0 &&
                    row->size == sequence_length &&
                    row->strand == '+' &&
                    row->source_size == sequence_length &&
                    row->dna == dna,
                "MAF export did not recover after the isolated write failure");
    }

    const auto recovered_hal = temp / "maf-cache-recovered.hal";
    RaMesh::hal_export::exportToHal(
        *failure_prepared, recovered_hal, inputs.managers,
        NewickParser("(mafLimit:0.1,mafPeer:0.1)mafRecoveryRoot;"),
        "mafRecoveryRoot", SoftMask::loadIndexes(inputs.softmask_paths));
    hal::AlignmentPtr alignment =
        hal::openHalAlignment(recovered_hal.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment),
            "cannot reopen HAL after failed shared preparation");
    const hal::Genome* a = alignment->openGenome("mafLimit");
    const hal::Genome* b = alignment->openGenome("mafPeer");
    requireLeafContents(a, {{"long.contig", dna}});
    requireLeafContents(b, {{"long.contig", dna}});
    require(pairCoverage(a, b) == sequence_length &&
                pairCoverage(b, a) == sequence_length,
            "failed shared preparation poisoned later HAL relationships");
    alignment->closeGenome(b);
    alignment->closeGenome(a);
    alignment.reset();
}

void testDirectoryOutputsPreserveExistingContents(
    const std::filesystem::path& temp) {
    FixtureInputs inputs;
    addLeaf(inputs, temp, "directoryA", {{"chr1", "ACGTACGT"}});
    addLeaf(inputs, temp, "directoryB", {{"chr1", "ACGTTCGT"}});
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks;
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "directory-prepared");
    const auto indexes = SoftMask::loadIndexes(inputs.softmask_paths);
    for (const bool hal_output : {false, true}) {
        const auto destination =
            temp / (hal_output ? "directory.hal" : "directory.maf");
        std::filesystem::create_directory(destination);
        const auto sentinel = destination / "keep.txt";
        const std::string contents = "pre-existing directory contents\n";
        writeFileBytes(sentinel, contents);
        bool rejected = false;
        try {
            if (hal_output) {
                RaMesh::hal_export::exportToHal(
                    *prepared, destination, inputs.managers,
                    NewickParser("(directoryA:0.1,directoryB:0.1)directoryRoot;"),
                    "directoryRoot", indexes);
            } else {
                RaMesh::hal_export::exportToMaf(
                    *prepared, destination, inputs.managers, false);
            }
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "export did not reject a directory destination");
        require(std::filesystem::is_directory(destination) &&
                    readFileBytes(sentinel) == contents,
                "rejected export changed or relocated the directory contents");
        require(!std::filesystem::exists(destination.string() + ".replace-backup") &&
                    !std::filesystem::exists(destination.string() + ".tmp"),
                "rejected directory export left transaction artifacts");
    }
}

void testHdf5WriteFailureIsTransactionalAndRecoverable(
    const std::filesystem::path& temp) {
    constexpr size_t sequence_length = 128 * 1024 + 3;
    std::string failure_dna(sequence_length, 'A');
    uint32_t state = 0x6d2b79f5U;
    constexpr char bases[] = "ACGT";
    for (char& base : failure_dna) {
        state = state * 1664525U + 1013904223U;
        base = bases[(state >> 30U) & 3U];
    }

    FixtureInputs inputs;
    const std::vector<std::pair<std::string, std::string>> limited_leaf{
        {"limit.contig", failure_dna}};
    const std::vector<std::pair<std::string, std::string>> limited_peer{
        {"limit.peer", "TgCaCC"}};
    addLeaf(inputs, temp, "limitLeaf", limited_leaf);
    addLeaf(inputs, temp, "limitPeer", limited_peer);
    const auto indexes = SoftMask::loadIndexes(inputs.softmask_paths);
    const std::vector<std::weak_ptr<RaMesh::Block>> blocks;
    const auto prepared =
        prepareFixtureInput(blocks, inputs, temp / "prepared-failure");
    const auto tree =
        "(limitLeaf:0.1,limitPeer:0.1)limitRoot;";
    RaMesh::hal_export::ExportConfig config;
    config.parallel_threads = 4;
    ScopedOpenMpConfiguration caller_threads(2);
    const int caller_thread_budget = omp_get_max_threads();

    const auto baseline_path = temp / "hdf5-baseline.hal";
    RaMesh::hal_export::exportToHal(
        *prepared, baseline_path, inputs.managers, NewickParser(tree),
        "limitRoot", indexes, config);
    {
        hal::AlignmentPtr alignment =
            hal::openHalAlignment(
                baseline_path.string(), nullptr, hal::READ_ACCESS);
        require(static_cast<bool>(alignment),
                "cannot reopen unlimited baseline HAL");
        const hal::Genome* leaf = alignment->openGenome("limitLeaf");
        requireLeafContents(leaf, limited_leaf);
        alignment->closeGenome(leaf);
    }
    require(std::filesystem::file_size(baseline_path) > 4096,
            "unlimited baseline HAL does not exceed the file-size cap");

    const auto destination = temp / "hdf5-transaction.hal";
    constexpr char original_data[] =
        "pre-existing HAL destination\0with binary bytes";
    const std::string original_bytes(
        original_data, sizeof(original_data) - 1);
    writeFileBytes(destination, original_bytes);
    std::filesystem::path temporary = destination;
    temporary += ".tmp";

    bool caught_std_exception = false;
    bool caught_nonstd_exception = false;
    {
        ScopedFileSizeLimit file_size_limit(4096);
        try {
            RaMesh::hal_export::exportToHal(
                *prepared, destination, inputs.managers, NewickParser(tree),
                "limitRoot", indexes, config);
        } catch (const std::exception&) {
            caught_std_exception = true;
        } catch (...) {
            caught_nonstd_exception = true;
        }
        file_size_limit.restore();
    }
    require(!caught_nonstd_exception,
            "limited HDF5 write escaped outside std::exception");
    require(caught_std_exception,
            "limited HDF5 write unexpectedly completed");
    require(omp_get_max_threads() == caller_thread_budget,
            "failed parallel HAL write changed the caller thread budget");
    require(readFileBytes(destination) == original_bytes,
            "failed HAL export changed the pre-existing destination");
    require(!std::filesystem::exists(temporary),
            "failed HAL export retained its temporary output");

    RaMesh::hal_export::exportToHal(
        *prepared, destination, inputs.managers, NewickParser(tree),
        "limitRoot", indexes, config);
    require(omp_get_max_threads() == caller_thread_budget,
            "recovered parallel HAL write changed the caller thread budget");
    hal::AlignmentPtr alignment =
        hal::openHalAlignment(
            destination.string(), nullptr, hal::READ_ACCESS);
    require(static_cast<bool>(alignment),
            "cannot reopen HAL after restoring file-size limit");
    const hal::Genome* leaf = alignment->openGenome("limitLeaf");
    requireLeafContents(leaf, limited_leaf);
    alignment->closeGenome(leaf);
    alignment.reset();
}

}  // namespace

int main(int argc, char** argv) {
    const std::string hal_validate = argc > 1 ? argv[1] : "halValidate";
    try {
        TemporaryDirectory temp;
        testFailedGenomeClosePreservesCallerHandle(temp.path());
        testDirectoryOutputsPreserveExistingContents(temp.path());
        testParallelCompressedChunksAreReadableAndCoherent(temp.path());
        testCompressedArrayFallbacksPreservePayloads(temp.path());
        testSparseOccurrenceIdsAndRejectedTerminalEvidence();
        testReorientedMaskedAmbiguityKeepsHistoricalConsensus(
            temp.path(), hal_validate);
        testMixedCaseReverseIndelAndMissingChild(temp.path(), hal_validate);
        testPhysicalCopiesUseParalogyEdges(temp.path(), hal_validate);
        testParentProjectionPreservesChildContainers(temp.path(), hal_validate);
        testMoreThanSixtyFourTreeNodes(temp.path(), hal_validate);
        testNoHomologyPreservesLeafDNAAndUnmappedTops(temp.path());
        testBulkDnaEncodingAndPackedNibbleBoundaries(temp.path());
        testLeafDNAStreamsAcrossOneMiBBoundary(temp.path());
        testPreparedInputSurvivesGraphReleaseAndReplays(
            temp.path(), hal_validate);
        testMafWriteFailureIsTransactionalAndRecoverable(temp.path());
        testHdf5WriteFailureIsTransactionalAndRecoverable(temp.path());
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "hal_memory_test: " << error.what() << '\n';
        return 1;
    }
}
