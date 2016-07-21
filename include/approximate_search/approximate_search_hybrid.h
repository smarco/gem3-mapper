/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_hybrid.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_HYBRID_H_
#define APPROXIMATE_SEARCH_HYBRID_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"

/*
 * Hybrid mapping [GEM-workflow 5.0]
 */
void approximate_search_hybrid(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Approximate Complete-Search based on PURE filtering+NS-search
 */
void approximate_search_hybrid_complete_search(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_HYBRID_H_ */
