/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_neighborhood.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_NEIGHBORHOOD_H_
#define APPROXIMATE_SEARCH_NEIGHBORHOOD_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"

/*
 * Neighborhood Generation (Inexact Search)
 */
void approximate_search_neighborhood_exact_search(
    approximate_search_t* const search,
    matches_t* const matches);
void approximate_search_neighborhood_inexact_search(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Neighborhood Search
 */
void approximate_search_neighborhood_search(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_NEIGHBORHOOD_H_ */
