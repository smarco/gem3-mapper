/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_region_profile.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_REGION_PROFILE_H_
#define APPROXIMATE_SEARCH_REGION_PROFILE_H_

#include "approximate_search/approximate_search.h"
#include "gpu/gpu_buffer_fmi_ssearch.h"
#include "gpu/gpu_buffer_fmi_asearch.h"

/*
 * Region profile generation
 */
typedef enum {
  // Adaptive Profile reducing the number of candidates (to a few) and using cut-off conditions
  region_profile_adaptive_lightweight,
  // Adaptive Profile reducing the number of candidates (to a some) and using cut-off conditions
  region_profile_adaptive_heavyweight,
  // Adaptive Profile limiting the max. length of a region (max-region-length)
  region_profile_adaptive_limited,
  // Adaptive Profile extracting all possible unique regions (least candidates possible)
  region_profile_adaptive_delimit,
} region_profile_strategy_t;

/*
 * Region Profile Adaptive
 */
void approximate_search_region_profile_adaptive(
    approximate_search_t* const search,
    const region_profile_strategy_t strategy,
    mm_stack_t* const mm_stack);

/*
 * Region Partition Fixed
 */
void approximate_search_region_profile_static_partition(approximate_search_t* const search);
void approximate_search_region_profile_static_compute(
    approximate_search_t* const search);

/*
 * Static Buffered Copy/Retrieve
 */
void approximate_search_region_profile_static_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void approximate_search_region_profile_static_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void approximate_search_region_profile_static_buffered_recompute(
    approximate_search_t* const search);

/*
 * Adaptive Buffered Copy/Retrieve
 */
void approximate_search_region_profile_adaptive_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
void approximate_search_region_profile_adaptive_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

#endif /* APPROXIMATE_SEARCH_REGION_PROFILE_H_ */
