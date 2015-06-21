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
  region_profile_adaptive_lightweight, // Adaptive Profile stopping whenever any cut-off cond. is reached
  region_profile_adaptive_limited,     // Adaptive Profile limiting the max. length of a region (max-region-length)
  region_profile_adaptive_extensive,   // Adaptive Profile extensive extracting all possible regions
} region_profile_generation_strategy_t;

/*
 * Stats
 */
GEM_INLINE void approximate_search_generate_region_profile_minimal_stats(approximate_search_t* const search);
GEM_INLINE void approximate_search_generate_region_profile_boost_stats(approximate_search_t* const search);
GEM_INLINE void approximate_search_generate_region_profile_delimit_stats(approximate_search_t* const search);

/*
 * Region profile generation
 */
GEM_INLINE void approximate_search_generate_region_profile(
    approximate_search_t* const search,const region_profile_model_t* const profile_model,
    const region_profile_generation_strategy_t region_profile_generation_strategy,
    const uint64_t min_regions,const bool allow_zero_regions);

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
