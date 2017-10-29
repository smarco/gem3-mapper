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
 */

#ifndef FILTERING_CANDIDATES_ACCESSORS_H_
#define FILTERING_CANDIDATES_ACCESSORS_H_

#include "utils/essentials.h"
#include "archive/archive.h"
#include "archive/search/archive_search_se_parameters.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/region/filtering_region.h"
#include "filtering/region/filtering_region_cache.h"
#include "matches/align/match_alignment_region.h"

/*
 * Counting
 */
uint64_t filtering_candidates_count_regions_by_status(
    const filtering_candidates_t* const filtering_candidates,
    const filtering_region_status_t filtering_region_status);

/*
 * Adding Positions (Candidate Positions)
 */
void filtering_candidates_add_positions_from_interval(
    filtering_candidates_t* const filtering_candidates,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    const uint64_t interval_lo,
    const uint64_t interval_hi,
    const uint64_t region_begin_pos,
    const uint64_t region_end_pos,
    const uint64_t region_errors);

/*
 * Adding Region (filtering regions)
 */
void filtering_candidates_add_region_from_group_positions(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const uint64_t first_candidate_idx,
    const uint64_t last_candidate_idx,
    const uint64_t align_distance,
    const uint64_t align_offset,
    const bool compose_alignment_regions,
    const bool run_length_text);
void filtering_candidates_add_region_verified(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const uint64_t text_begin_offset,
    const uint64_t text_end_offset,
    const uint64_t begin_position,
    const uint64_t end_position,
    const uint64_t align_distance);

#endif /* FILTERING_CANDIDATES_ACCESSORS_H_ */
