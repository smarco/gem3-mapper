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

#ifndef MATCH_ALIGNMENT_REGION_RL_H_
#define MATCH_ALIGNMENT_REGION_RL_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/align/match_alignment_region.h"

/*
 * RL-Space Translate CIGAR
 */
void match_alignment_region_rl_translate(
    match_alignment_region_t* const match_alignment_region,
    uint32_t* const rl_key_runs_acc,
    uint32_t* const rl_text_runs_acc,
    const bool left_gap_alignment,
    vector_t* const cigar_vector);

#endif /* MATCH_ALIGNMENT_REGION_RL_H_ */
