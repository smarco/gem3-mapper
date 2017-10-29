/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Approximate-String-Matching (ASM) module encapsulating
 *   the basic search-stages that many ASM approaches use.
 */

#ifndef APPROXIMATE_SEARCH_STAGES_H_
#define APPROXIMATE_SEARCH_STAGES_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"

/*
 * Basic Cases Selector
 */
void approximate_search_begin(
    approximate_search_t* const search);

/*
 * Exact Filtering Adaptive
 */
void approximate_search_exact_filtering_adaptive(
    approximate_search_t* const search,
    matches_t* const matches);
void approximate_search_exact_filtering_adaptive_cutoff(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Filtering Verification (+ realign)
 */
void approximate_search_verify(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Neighborhood Search
 */
void approximate_search_neighborhood(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Unbound Filtering (+ realign)
 */
void approximate_search_align_local(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * End of the search
 */
void approximate_search_end(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_STAGES_H_ */
