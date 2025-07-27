#include "sequence_utils.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <limits>

namespace SequenceUtils {

template<typename ManagerType>
void recordReferenceSequenceStats(const std::string& species_name, 
                                 const std::unique_ptr<ManagerType>& manager,
                                 SeqPro::Length& reference_min_seq_length) {
    auto seq_count = manager->getSequenceCount();
    
    if (species_name == "reference") {
        // Only compute sequence length statistics for reference
        auto seq_names = manager->getSequenceNames();
        SeqPro::Length min_seq_length = std::numeric_limits<SeqPro::Length>::max();
        SeqPro::Length max_seq_length = 0;
        
        for (const auto& seq_name : seq_names) {
            auto seq_length = manager->getSequenceLength(seq_name);
            min_seq_length = std::min(min_seq_length, seq_length);
            max_seq_length = std::max(max_seq_length, seq_length);
        }
        
        spdlog::info("[{}] Loaded {} sequences, min length: {}, max length: {}", 
                     species_name, seq_count, min_seq_length, max_seq_length);
        reference_min_seq_length = min_seq_length;
    } else {
        spdlog::info("[{}] Loaded {} sequences", species_name, seq_count);
    }
}

double buildRefGlobalCache(const SeqPro::ManagerVariant& manager_variant,
                          SeqPro::Length sampling_interval,
                          sdsl::int_vector<0>& ref_global_cache) {
    // 目前没有清空，原理上应该是不用清空的
    auto t_start_cache = std::chrono::steady_clock::now();
    
    // Get total length and calculate cache size
    auto total_length = std::visit([](auto&& manager_ptr) {
        using T = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getTotalLength();
        }
        else if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getTotalLengthWithSeparators();
        }
    }, manager_variant);
   
    auto cache_size = (total_length / sampling_interval) + 1;
    
    spdlog::info("Building ref_global_cache, sampling_interval={}, cache_size={}", 
                 sampling_interval, cache_size);
    
    ref_global_cache.resize(cache_size);
    
    std::visit([&](auto&& manager_ptr) {
        using T = std::decay_t<decltype(manager_ptr)>;
        
        // Get all sequence information, sorted by global_start_pos
        auto seq_names = manager_ptr->getSequenceNames();
        std::vector<const SeqPro::SequenceInfo*> seq_infos;
        seq_infos.reserve(seq_names.size());
        
        for (const auto& name : seq_names) {
            if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::SequenceManager>>) {
                const auto* info = manager_ptr->getIndex().getSequenceInfo(name);
                if (info) seq_infos.push_back(info);
            } else if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                const auto* info = manager_ptr->getOriginalManager().getIndex().getSequenceInfo(name);
                if (info) seq_infos.push_back(info);
            }
        }
        
        // Sort by global start position
        std::sort(seq_infos.begin(), seq_infos.end(), 
                  [](const SeqPro::SequenceInfo* a, const SeqPro::SequenceInfo* b) {
                      return a->global_start_pos < b->global_start_pos;
                  });
        
        // Sequential filling: optimized to avoid repeated binary searches
        size_t current_seq_idx = 0;
        for (SeqPro::Position i = 0; i < cache_size; ++i) {
            SeqPro::Position sample_global_pos = i * sampling_interval;
            
            if (sample_global_pos >= total_length) {
                ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
                continue;
            }
            
            // Find the sequence containing the current position
            while (current_seq_idx < seq_infos.size()) {
                const auto* current_seq = seq_infos[current_seq_idx];
                SeqPro::Position seq_end = current_seq->global_start_pos + current_seq->length + 1;
                
                // Add separator count only for MaskedSequenceManager
                if constexpr (std::is_same_v<T, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                    seq_end += manager_ptr->getSeparatorCount(current_seq->id);
                }
                
                if (sample_global_pos >= current_seq->global_start_pos && sample_global_pos < seq_end) {
                    // Found the sequence containing this position
                    ref_global_cache[i] = current_seq->id;
                    break;
                } else if (sample_global_pos >= seq_end) {
                    // Current sequence is past, move to next sequence
                    current_seq_idx++;
                } else {
                    // sample_global_pos < current_seq->global_start_pos, shouldn't happen
                    spdlog::warn("Unexpected coordinate order: sample_pos={}, seq_start={}", 
                                sample_global_pos, current_seq->global_start_pos);
                    ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
                    break;
                }
            }
            
            // If we've gone through all sequences without finding one, mark as invalid
            if (current_seq_idx >= seq_infos.size()) {
                ref_global_cache[i] = SeqPro::SequenceIndex::INVALID_ID;
            }
        }
    }, manager_variant);
    
    // Compress the cache to save memory
    sdsl::util::bit_compress(ref_global_cache);
    
    auto t_end_cache = std::chrono::steady_clock::now();
    std::chrono::duration<double> cache_time = t_end_cache - t_start_cache;
    spdlog::info("ref_global_cache building completed in {:.3f} seconds", cache_time.count());
    
    return cache_time.count();
}

void recordReferenceSequenceStats(const std::string& species_name, 
                                 const SeqPro::SharedManagerVariant& shared_manager,
                                 SeqPro::Length& reference_min_seq_length) {
    std::visit([&](auto&& manager_ptr) {
        recordReferenceSequenceStats(species_name, manager_ptr, reference_min_seq_length);
    }, *shared_manager);
}

double buildRefGlobalCache(const SeqPro::SharedManagerVariant& shared_manager_variant,
                          SeqPro::Length sampling_interval,
                          sdsl::int_vector<0>& ref_global_cache) {
    // Delegate to the ManagerVariant version by dereferencing the shared_ptr
    return buildRefGlobalCache(*shared_manager_variant, sampling_interval, ref_global_cache);
}

// Explicit template instantiations
template void recordReferenceSequenceStats<SeqPro::SequenceManager>(
    const std::string&, const std::unique_ptr<SeqPro::SequenceManager>&, SeqPro::Length&);
template void recordReferenceSequenceStats<SeqPro::MaskedSequenceManager>(
    const std::string&, const std::unique_ptr<SeqPro::MaskedSequenceManager>&, SeqPro::Length&);

} // namespace SequenceUtils 