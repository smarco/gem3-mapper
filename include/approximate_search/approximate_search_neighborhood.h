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
 *   Approximate-String-Matching (ASM) using neighborhood-search (NS)
 */

#ifndef APPROXIMATE_SEARCH_NEIGHBORHOOD_H_
#define APPROXIMATE_SEARCH_NEIGHBORHOOD_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"

/*
 * Exact Search
 */
void approximate_search_neighborhood_exact_search(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Neighborhood Search Brute Force
 */
void approximate_search_neighborhood_search_brute_force(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Neighborhood Search (Using partitions)
 */
void approximate_search_neighborhood_search_partition(
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Neighborhood Search (Using partitions + region-profile preconditioned)
 */
void approximate_search_neighborhood_search_partition_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_NEIGHBORHOOD_H_ */
