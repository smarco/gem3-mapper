/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_STAGES_H_
#define APPROXIMATE_SEARCH_FILTERING_STAGES_H_

#include "essentials.h"
#include "approximate_search.h"
#include "approximate_search_filtering_base.h"

/*
 * Basic Cases
 */
void approximate_search_adaptive_mapping_basic_cases(approximate_search_t* const search);

/*
 * Exact Filtering Fixed
 */
void approximate_search_exact_filtering_fixed(approximate_search_t* const search);

/*
 * Exact Filtering Adaptive
 */
void approximate_search_exact_filtering_adaptive_lightweight(
    approximate_search_t* const search,matches_t* const matches);
void approximate_search_exact_filtering_adaptive_recovery(
    approximate_search_t* const search,matches_t* const matches);
void approximate_search_exact_filtering_adaptive_cutoff(
    approximate_search_t* const search,matches_t* const matches);

/*
 * Boost Exact Filtering Adaptive
 */
void approximate_search_exact_filtering_boost(approximate_search_t* const search,matches_t* const matches);

/*
 * Inexact Filtering Adaptive
 */
void approximate_search_inexact_filtering(approximate_search_t* const search,matches_t* const matches);

/*
 * Unbound Filtering
 */
void approximate_search_unbounded_align(approximate_search_t* const search,matches_t* const matches);

/*
 * Filtering Verification (+ realign)
 */
GEM_INLINE void approximate_search_verify(approximate_search_t* const search,matches_t* const matches);
GEM_INLINE void approximate_search_verify_using_bpm_buffer(
    approximate_search_t* const search,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end,
    matches_t* const matches);


#endif /* APPROXIMATE_SEARCH_FILTERING_STAGES_H_ */
