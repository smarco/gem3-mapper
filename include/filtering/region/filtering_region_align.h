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
 *   Filtering region module provides functions to produce a full-alignment
 *   of a filtering-region against its corresponding text-region
 */

#ifndef FILTERING_REGION_ALIGN_H_
#define FILTERING_REGION_ALIGN_H_

#include "filtering/candidates/filtering_candidates.h"
#include "filtering/region/filtering_region.h"
#include "align/pattern/pattern.h"
#include "matches/matches.h"

/*
 * Region (Re)Align by clone previous
 */
void filtering_region_align_clone(
    match_trace_t* const match_trace_src,
    match_trace_t* const match_trace_dst,
    filtering_region_t* const filtering_region_dst);

/*
 * Region (Re)Align
 */
bool filtering_region_align(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches,
    match_trace_t* const match_trace);
bool filtering_region_align_local(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches,
    match_trace_t* const match_trace);

#endif /* FILTERING_REGION_ALIGN_H_ */
