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
 *   Smith-Waterman-Gotoh (SWG) alignment module
 */

#ifndef ALIGN_SWG_BANDED_H_
#define ALIGN_SWG_BANDED_H_

#include "utils/essentials.h"
#include "align/align_swg_score.h"
#include "matches/align/match_alignment.h"

/*
 * SWG Banded
 */
void align_swg_banded(
    const uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    uint64_t text_length,
    const swg_penalties_t* const swg_penalties,
    const uint64_t max_bandwidth,
    const bool begin_free,
    const bool end_free,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator);
void align_swg_banded_extend(
    const uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    uint64_t text_length,
    const bool reverse_extension,
    const bool local_extension,
    const int32_t base_score,
    const swg_penalties_t* const swg_penalties,
    const uint64_t max_bandwidth,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    uint64_t* const max_key_aligned,
    uint64_t* const max_text_aligned,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator);

#endif /* ALIGN_SWG_BANDED_H_ */
