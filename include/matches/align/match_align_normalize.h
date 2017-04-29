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

#ifndef MATCH_ALIGN_NORMALIZE_H_
#define MATCH_ALIGN_NORMALIZE_H_

#include "utils/essentials.h"
#include "archive/search/archive_search_se_parameters.h"
#include "matches/matches.h"

/*
 * SWG Normalize CIGAR & Adjust Position (Translate RL-CIGAR if required)
 */
void match_align_normalize(
    matches_t* const matches,
    match_trace_t* const match_trace,
    search_parameters_t* const search_parameters);

#endif /* MATCH_ALIGN_NORMALIZE_H_ */
