/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_control.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_CONTROL_H_
#define APPROXIMATE_SEARCH_FILTERING_CONTROL_H_

#include "essentials.h"
#include "approximate_search.h"

void approximate_search_adjust_max_differences_using_strata(
    approximate_search_t* const search,matches_t* const matches) ;

void asearch_compute_predictors(
    approximate_search_t* const search,matches_t* const matches,
    matches_predictors_t* const predictors);

bool asearch_fulfilled(approximate_search_t* const search,matches_t* const matches);

bool asearch_filter_ahead_candidates(approximate_search_t* const search,matches_t* const matches);

/*
 * Control DFA States
 */
bool asearch_control_trigger_boost(approximate_search_t* const search,matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_FILTERING_CONTROL_H_ */
