/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_region_profile.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_REGION_PROFILE_H_
#define APPROXIMATE_SEARCH_REGION_PROFILE_H_

#include "approximate_search.h"
#include "gpu_buffer_fmi_bsearch.h"

/*
 * Region profile generation
 */
typedef enum {
  // Adaptive Profile reducing the total number of candidates and using cut-off conditions
  region_profile_adaptive_lightweight,
  // Adaptive Profile on not-covered and excessive long regions (to boost sensitivity)
  region_profile_adaptive_boost,
  // Adaptive Profile limiting the max. length of a region (max-region-length)
  region_profile_adaptive_limited,
  // Adaptive Profile extracting all possible unique regions (least candidates possible)
  region_profile_adaptive_delimit,
  // Adaptive Profile with less strict thresholds as to recover as much as feasible from the read
  region_profile_adaptive_recovery,
} approximate_search_region_profile_strategy_t;

/*
 * Region Partition Fixed
 */
void approximate_search_region_partition_fixed(approximate_search_t* const search);

/*
 * Region Profile Adaptive
 */
void approximate_search_region_profile_adaptive(
    approximate_search_t* const search,
    const approximate_search_region_profile_strategy_t region_profile_strategy,
    mm_stack_t* const mm_stack);

/*
 * Buffered Copy/Retrieve
 */
void approximate_search_region_profile_buffered_copy(
    approximate_search_t* const search,gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);
void approximate_search_region_profile_buffered_retrieve(
    approximate_search_t* const search,gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);

#endif /* APPROXIMATE_SEARCH_REGION_PROFILE_H_ */
