/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_complete.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_COMPLETE_H_
#define APPROXIMATE_SEARCH_FILTERING_COMPLETE_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"

/*
 * Complete Search
 */
void approximate_search_filtering_complete(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_FILTERING_COMPLETE_H_ */
