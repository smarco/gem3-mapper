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
 *   Efficient Smith-Waterman-Gotoh (SWG) alignment module using SIMD instructions
 */

#ifndef ALIGN_SWG_SIMD_H_
#define ALIGN_SWG_SIMD_H_

#include "utils/essentials.h"
#include "align/align_swg_score.h"

/*
 * Constants
 */
#define SWG_SCORE_MIN (INT16_MIN)

/*
 * SWG Query Profile
 */
typedef struct {
  // 8-bits cells
  uint8_t match_bias_uint8;  // Query Profile Bias
  uint8_t matrix_bias_uint8; // Full matrix Bias
  bool overflow_uint8;
  uint64_t segment_length_uint8;
  uint64_t key_effective_length_uint8;
  uint8_t* query_profile_uint8[DNA__N_RANGE];
  // 16-bits cells
  uint64_t segment_length_int16;
  uint64_t key_effective_length_int16;
  int16_t* query_profile_int16[DNA__N_RANGE];
} swg_query_profile_t;

/*
 * Init SWG Query Profile
 */
void align_swg_query_profile_init(
    swg_query_profile_t* const swg_query_profile,
    const swg_penalties_t* swg_penalties,
    const uint64_t max_expected_key_length,
    mm_allocator_t* const mm_allocator);
bool align_swg_query_profile_compile_uint8(
    swg_query_profile_t* const swg_query_profile,
    const swg_penalties_t* swg_penalties,
    const uint8_t* const key,
    const uint64_t key_length);
bool align_swg_query_profile_compile_int16(
    swg_query_profile_t* const swg_query_profile,
    const swg_penalties_t* swg_penalties,
    const uint8_t* const key,
    const uint64_t key_length);

#endif /* ALIGN_SWG_SIMD_H_ */
