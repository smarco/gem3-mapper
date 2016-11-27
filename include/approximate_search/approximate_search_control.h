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
 *   Approximate-String-Matching (ASM) search control functions (regulate
 *   the depth of the search)
 */

#ifndef APPROXIMATE_SEARCH_CONTROL_H_
#define APPROXIMATE_SEARCH_CONTROL_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"
#include "matches/classify/matches_predictors.h"

/*
 * Pattern test
 */
bool asearch_control_test_pattern(
    approximate_search_t* const search);

/*
 * Search Limits
 */
void asearch_control_adjust_current_max_error(
    approximate_search_t* const search,
    matches_t* const matches);
bool asearch_control_max_matches_reached(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Accuracy test
 */
bool asearch_control_test_accuracy__adjust_depth(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Local-alignment
 */
bool asearch_control_test_local_alignment(
    approximate_search_t* const search,
    matches_t* const matches);


#endif /* APPROXIMATE_SEARCH_CONTROL_H_ */
