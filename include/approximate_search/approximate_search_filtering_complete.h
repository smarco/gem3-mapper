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
 *   Approximate-String-Matching (ASM) using modified adaptive-filtering techniques
 *   as to produce complete results (complete AF)
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_COMPLETE_H_
#define APPROXIMATE_SEARCH_FILTERING_COMPLETE_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"

/*
 * Approximate Complete-Search based on filtering
 */
void approximate_search_filtering_complete(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Approximate Complete-Search based on filtering+NS-search
 */
void approximate_search_hybrid_complete_search(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_FILTERING_COMPLETE_H_ */
