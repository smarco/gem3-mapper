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
 */

#include "filtering/candidates/filtering_candidates.h"
#include "align/pattern/pattern.h"

/*
 * Filtering candidates classify as subdominant match
 * (Expected best score worse than all the found matches so far)
 */
bool filtering_candidates_classify_subdominant_match(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches);

/*
 * Filtering candidates classify as subdominant region
 * compared with a reference score/distance
 */
bool filtering_candidates_classify_subdominant_region_edit(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t filtering_region_rank,
    const uint64_t filtering_region_edit_bound,
    const uint64_t reference_edit_bound);
bool filtering_candidates_classify_subdominant_region_swg(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t filtering_region_rank,
    const int32_t filtering_region_score_bound,
    const int32_t reference_score_bound);

