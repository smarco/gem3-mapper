/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_adaptive.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_H_
#define APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_H_

#include "essentials.h"
#include "approximate_search.h"

/*
 * [GEM-workflow 4.0] Adaptive mapping
 *
 *   Filtering-only approach indented to adjust the degree of filtering w.r.t
 *   the structure of the read. Thus, in general terms, a read with many regions
 *   will enable this approach to align the read up to more mismatches than a read
 *   with less number of regions.
 *   Fast-mapping (in all its kinds) tries to detect the proper degree of filtering
 *   to achieve a compromise between speed and depth of the search (max_mismatches)
 */
void approximate_search_filtering_adaptive(approximate_search_t* const search,matches_t* const matches);

// Reduced workflows (step-wise as to use callbacks)
void approximate_search_filtering_adaptive_generate_regions(approximate_search_t* const search);
void approximate_search_filtering_adaptive_generate_candidates(approximate_search_t* const search);

#endif /* APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_H_ */
