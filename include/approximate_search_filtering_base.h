/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_base.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_BASE_H_
#define APPROXIMATE_SEARCH_FILTERING_BASE_H_

#include "approximate_search.h"

/*
 * Region profile generation
 */
typedef enum {
  region_profile_adaptive_lightweight, // Adaptive Profile reducing the total number of candidates and using cut-off conditions
  region_profile_adaptive_boost,       // Adaptive Profile on not-covered and excessive long regions (to boost sensitivity)
  region_profile_adaptive_limited,     // Adaptive Profile limiting the max. length of a region (max-region-length)
  region_profile_adaptive_delimit,     // Adaptive Profile extracting all possible unique regions (least candidates possible)
  region_profile_adaptive_recovery,    // Adaptive Profile with less strict thresholds as to recover as much as feasible from the read
} region_profiling_strategy_t;

/*
 * Region profile generation
 */
GEM_INLINE void approximate_search_generate_region_profile(
    approximate_search_t* const search,const region_profiling_strategy_t rp_strategy,
    mm_stack_t* const mm_stack);

/*
 * Generate Candidates
 */
GEM_INLINE void approximate_search_generate_exact_candidates(
    approximate_search_t* const search,matches_t* const matches);
GEM_INLINE void approximate_search_generate_inexact_candidates(
    approximate_search_t* const search,const bool dynamic_scheduling,
    const bool verify_ahead,matches_t* const matches);

/*
 * Verify Candidates
 */
GEM_INLINE void approximate_search_verify_candidates(approximate_search_t* const search,matches_t* const matches);

/*
 * Neighborhood Generation (Inexact Search)
 */
GEM_INLINE void approximate_search_neighborhood_exact_search(approximate_search_t* const search,matches_t* const matches);
GEM_INLINE void approximate_search_neighborhood_inexact_search(approximate_search_t* const search,matches_t* const matches);


#endif /* APPROXIMATE_SEARCH_FILTERING_BASE_H_ */
