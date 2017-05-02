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

#ifndef MATCH_ALIGNMENT_H_
#define MATCH_ALIGNMENT_H_

#include "utils/essentials.h"

/*
 * Match Alignment
 */
typedef struct {
  uint64_t match_text_offset; // Match text offset (wrt beginning of text-candidate)
  uint64_t match_position;    // Match position
  uint64_t cigar_offset;      // CIGAR offset in buffer
  uint64_t cigar_length;      // CIGAR length
  int64_t effective_length;   // Match effective length
  int32_t score;              // Score assigned by the aligner
} match_alignment_t;
/*
 * Alignment Model
 */
typedef enum {
  match_alignment_model_none,
  match_alignment_model_hamming,
  match_alignment_model_levenshtein,
  match_alignment_model_gap_affine
} match_alignment_model_t;

/*
 * Check
 */
bool match_alignment_check(
    FILE* const stream,
    match_alignment_t* const match_alignment,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    const bool verbose,
    mm_allocator_t* const mm_allocator);
/*
 *  Display
 */
void match_alignment_print_pretty(
    FILE* const stream,
    match_alignment_t* const match_alignment,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator);

#endif /* MATCH_ALIGNMENT_H_ */
