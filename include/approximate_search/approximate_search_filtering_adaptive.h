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
 *   Approximate-String-Matching (ASM) using adaptive-filtering techniques (AF)
 */

#ifndef APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_H_
#define APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"

/*
 * Adaptive mapping [GEM-workflow 4.0]
 *
 *   Filtering-only approach indented to adjust the degree of filtering w.r.t
 *   the structure of the read. Thus, in general terms, a read with many regions
 *   will enable this approach to align the read up to more mismatches than a read
 *   with less number of regions.
 *   Fast-mapping (in all its kinds) tries to detect the proper degree of filtering
 *   to achieve a compromise between speed and depth of the search (max_mismatches)
 */
void approximate_search_filtering_adaptive(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_FILTERING_ADAPTIVE_H_ */
