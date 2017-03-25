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
 *   Filtering region module provides functions to configure all the parameters
 *   involved in a full-alignment process of a filtering-region against a text-region
 */

#ifndef FILTERING_REGION_ALIGN_CONFIGURE_H_
#define FILTERING_REGION_ALIGN_CONFIGURE_H_

#include "filtering/region/filtering_region.h"
#include "align/pattern/pattern.h"
#include "archive/search/archive_search_se_parameters.h"

/*
 * Configure Basic Alignment
 */
void filtering_region_align_configure_exact(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern);
void filtering_region_align_configure_hamming(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern);
void filtering_region_align_configure_levenshtein(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    const bool left_gap_alignment,
    mm_allocator_t* const mm_allocator);

/*
 * Configure SWG-based Alignment
 */
void filtering_region_align_configure_swg(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    const bool left_gap_alignment,
    const bool local_alignment,
    mm_allocator_t* const mm_allocator);

#endif /* FILTERING_REGION_ALIGN_CONFIGURE_H_ */
