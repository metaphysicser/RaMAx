#include "anchor.h"
#include "SeqPro.h"
#include <algorithm>
#include <spdlog/spdlog.h>

// ------------------------------------------------------------------
// 智能序列分块函数实现
// ------------------------------------------------------------------

RegionVec preAllocateChunks(const SeqPro::ManagerVariant& seq_manager,
                           uint_t chunk_size,
                           uint_t overlap_size,
                           size_t max_sequences_threshold,
                           uint_t min_chunk_size) {
    
    size_t seq_count = std::visit([](auto&& manager_ptr) -> size_t {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getSequenceCount();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getSequenceCount();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);

    uint64_t total_length = std::visit([](auto&& manager_ptr) -> uint64_t {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getTotalLength();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getTotalLength();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);
    
    // 智能策略选择
    bool use_sequence_based_strategy = false;
    
    // 策略1: 如果序列数量很多，优先按序列划分
    if (seq_count >= max_sequences_threshold) {
        use_sequence_based_strategy = true;
        spdlog::info("Using sequence-based chunking strategy: {} sequences (>= threshold {})", 
                    seq_count, max_sequences_threshold);
    }
    // 策略2: 如果序列数量较少但每个序列都很短，也按序列划分
    else if (seq_count > 0 && total_length / seq_count < chunk_size) {
        use_sequence_based_strategy = true;
        spdlog::info("Using sequence-based chunking strategy: average sequence length ({}) < chunk_size ({})", 
                    total_length / seq_count, chunk_size);
    }
    // 策略3: 默认按大小划分
    else {
        spdlog::info("Using size-based chunking strategy: {} sequences, total {} bases", 
                    seq_count, total_length);
    }
    
    if (use_sequence_based_strategy) {
        return preAllocateChunksBySequence(seq_manager);
    } else {
        return preAllocateChunksBySize(seq_manager, chunk_size, overlap_size, min_chunk_size);
    }
}

RegionVec preAllocateChunksBySequence(const SeqPro::ManagerVariant& seq_manager) {
    RegionVec chunks;
    auto seq_names = std::visit([](auto&& manager_ptr) -> std::vector<std::string> {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getSequenceNames();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getSequenceNames();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);

    chunks.reserve(seq_names.size());
    
    for (const auto& seq_name : seq_names) {
        uint64_t seq_length = std::visit([&seq_name](auto&& manager_ptr) -> uint64_t {
            using PtrType = std::decay_t<decltype(manager_ptr)>;
            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                return manager_ptr->getSequenceLength(seq_name);
            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                return manager_ptr->getSequenceLength(seq_name);
            } else {
                throw std::runtime_error("Unhandled manager type in variant.");
            }
        }, seq_manager);

        chunks.push_back({seq_name, 0, static_cast<uint_t>(seq_length)});
    }
    
    spdlog::info("Generated {} chunks by sequence (one chunk per sequence)", chunks.size());
    return chunks;
}

RegionVec preAllocateChunksBySize(const SeqPro::ManagerVariant& seq_manager,
                                 uint_t chunk_size,
                                 uint_t overlap_size,
                                 uint_t min_chunk_size) {
    RegionVec chunks;
    
    // 预估总chunk数以减少realloc
    uint64_t total_length = std::visit([](auto&& manager_ptr) -> uint64_t {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getTotalLength();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getTotalLength();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);
    size_t est_chunks = (total_length + chunk_size - 1) / chunk_size;
    chunks.reserve(est_chunks);
    
    auto seq_names = std::visit([](auto&& manager_ptr) -> std::vector<std::string> {
        using PtrType = std::decay_t<decltype(manager_ptr)>;
        if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
            return manager_ptr->getSequenceNames();
        } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
            return manager_ptr->getSequenceNames();
        } else {
            throw std::runtime_error("Unhandled manager type in variant.");
        }
    }, seq_manager);
    
    // 对每条序列分别切分
    for (const auto& seq_name : seq_names) {
        uint64_t seq_length = std::visit([&seq_name](auto&& manager_ptr) -> uint64_t {
            using PtrType = std::decay_t<decltype(manager_ptr)>;
            if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::SequenceManager>>) {
                return manager_ptr->getSequenceLength(seq_name);
            } else if constexpr (std::is_same_v<PtrType, std::unique_ptr<SeqPro::MaskedSequenceManager>>) {
                return manager_ptr->getSequenceLength(seq_name);
            } else {
                throw std::runtime_error("Unhandled manager type in variant.");
            }
        }, seq_manager);

        // 如果序列长度小于等于chunk_size或min_chunk_size，则只生成一个不重叠chunk
        if (seq_length <= chunk_size || seq_length <= min_chunk_size) {
            chunks.push_back({seq_name, 0, static_cast<uint_t>(seq_length)});
            continue;
        }
        
        // 多个chunk，需要在它们之间保留overlap_size
        uint64_t start = 0;
        while (start < seq_length) {
            // 本chunk的实际长度（最后一块可能不足chunk_size）
            uint64_t this_len = std::min<uint64_t>(chunk_size, seq_length - start);
            chunks.push_back({seq_name, static_cast<uint_t>(start), static_cast<uint_t>(this_len)});
            
            // 计算下一个chunk的起始位置：前进chunk_size，然后回退overlap_size
            if (start + this_len >= seq_length) {
                break;
            }
            start += chunk_size;
            // 保证不回退超过已经走过的距离
            start = (start >= overlap_size ? start - overlap_size : 0);
        }
    }
    
    spdlog::info("Generated {} chunks by size (chunk_size: {}, overlap_size: {})", 
                chunks.size(), chunk_size, overlap_size);
    return chunks;
} 