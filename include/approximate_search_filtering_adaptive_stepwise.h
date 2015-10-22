/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_adaptive_stepwise.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_STEPWISE_H_
#define APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_STEPWISE_H_

#include "essentials.h"
#include "approximate_search.h"

/*
 * Adaptive mapping (Stepwise)
 */
void approximate_search_filtering_adaptive_stepwise_generate_regions(approximate_search_t* const search);
void approximate_search_filtering_adaptive_stepwise_generate_candidates(approximate_search_t* const search);
void approximate_search_filtering_adaptive_stepwise_finish(approximate_search_t* const search,matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_STEPWISE_H_ */
