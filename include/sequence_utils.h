#ifndef SEQUENCE_UTILS_H
#define SEQUENCE_UTILS_H

#include "SeqPro.h"
#include "config.hpp"
#include <sdsl/int_vector.hpp>
#include <memory>
#include <chrono>

namespace SequenceUtils {

/**
 * @brief Record reference sequence statistics and compute minimum sequence length
 * 
 * @tparam ManagerType Type of sequence manager (SequenceManager or MaskedSequenceManager)
 * @param species_name Name of the species (used for logging)
 * @param manager Pointer to the sequence manager
 * @param reference_min_seq_length Output parameter for minimum sequence length
 */
template<typename ManagerType>
void recordReferenceSequenceStats(const std::string& species_name, 
                                 const std::unique_ptr<ManagerType>& manager,
                                 SeqPro::Length& reference_min_seq_length);

/**
 * @brief Record reference sequence statistics for SharedManagerVariant
 * 
 * @param species_name Name of the species (used for logging)
 * @param shared_manager Shared pointer to the manager variant
 * @param reference_min_seq_length Output parameter for minimum sequence length
 */
void recordReferenceSequenceStats(const std::string& species_name, 
                                 const SeqPro::SharedManagerVariant& shared_manager,
                                 SeqPro::Length& reference_min_seq_length);

/**
 * @brief Build ref_global_cache for optimized sequence ID lookup
 * 
 * This function pre-computes sequence IDs for sampling points to avoid repeated 
 * binary searches in globalToLocal() calls during alignment.
 * 
 * @param manager_variant Variant containing either SequenceManager or MaskedSequenceManager
 * @param sampling_interval Interval between sampling points  
 * @param ref_global_cache Output cache vector to be filled
 * @return Time taken to build the cache in seconds
 */
double buildRefGlobalCache(const SeqPro::ManagerVariant& manager_variant,
                          SeqPro::Length sampling_interval,
                          sdsl::int_vector<0>& ref_global_cache);

/**
 * @brief Build ref_global_cache for optimized sequence ID lookup (SharedManagerVariant overload)
 * 
 * This function pre-computes sequence IDs for sampling points to avoid repeated 
 * binary searches in globalToLocal() calls during alignment.
 * 
 * @param shared_manager_variant Shared pointer to variant containing either SequenceManager or MaskedSequenceManager
 * @param sampling_interval Interval between sampling points  
 * @param ref_global_cache Output cache vector to be filled
 * @return Time taken to build the cache in seconds
 */
double buildRefGlobalCache(const SeqPro::SharedManagerVariant& shared_manager_variant,
                          SeqPro::Length sampling_interval,
                          sdsl::int_vector<0>& ref_global_cache);

} // namespace SequenceUtils

#endif // SEQUENCE_UTILS_H 