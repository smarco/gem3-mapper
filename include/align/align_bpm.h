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
 *   Alignment module using BPM-algorithm to compute levenshtein full-alignment
 *   (Myers' Fast Bit-Vector algorithm to compute levenshtein full-alignment)
 */

#ifndef ALIGN_BPM_H_
#define ALIGN_BPM_H_

#include "utils/essentials.h"
#include "matches/align/match_alignment.h"
#include "matches/align/match_align_dto.h"

/*
 * BPM matrix
 */
typedef struct {
  uint64_t* Pv;
  uint64_t* Mv;
  uint64_t min_score;
  uint64_t min_score_column;
} bpm_align_matrix_t;

/*
 * BPM. Compute BPM-DP-Matrix
 */
void align_bpm_compute_matrix(
    match_align_input_t* const align_input,
    uint64_t max_distance,
    bpm_align_matrix_t* const bpm_align_matrix,
    mm_stack_t* const mm_stack);
/*
 * BPM. Recover CIGAR from a matching string
 */
void align_bpm_backtrace_matrix(
    match_align_input_t* const align_input,
    const bool left_gap_alignment,
    bpm_align_matrix_t* const bpm_align_matrix,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector);
/*
 * BPM Align match
 */
void align_bpm_match(
    match_align_input_t* const align_input,
    const uint64_t max_distance,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack);

#endif /* ALIGN_BPM_H_ */
