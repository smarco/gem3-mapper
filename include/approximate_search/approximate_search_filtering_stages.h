/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_STAGES_H_
#define APPROXIMATE_SEARCH_FILTERING_STAGES_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"

/*
 * Exact Filtering Adaptive
 */
void approximate_search_exact_filtering_adaptive_lightweight(
    approximate_search_t* const restrict search,
    matches_t* const restrict matches);
void approximate_search_exact_filtering_adaptive_heavyweight(
    approximate_search_t* const restrict search,
    matches_t* const restrict matches);
void approximate_search_exact_filtering_adaptive_cutoff(
    approximate_search_t* const restrict search,
    matches_t* const restrict matches);

/*
 * Inexact Filtering Adaptive
 */
void approximate_search_inexact_filtering(
    approximate_search_t* const restrict search,
    matches_t* const restrict matches);

/*
 * Filtering Verification (+ realign)
 */
void approximate_search_verify(
    approximate_search_t* const restrict search,
    matches_t* const restrict matches);

/*
 * Unbound Filtering (+ realign)
 */
void approximate_search_align_local(
    approximate_search_t* const restrict search,
    matches_t* const restrict matches);

/*
 * End of the search
 */
void approximate_search_end(
    approximate_search_t* const restrict search,
    matches_t* const restrict matches);

#endif /* APPROXIMATE_SEARCH_FILTERING_STAGES_H_ */
