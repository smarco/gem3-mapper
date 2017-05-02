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

#ifndef ALIGN_SWG_H_
#define ALIGN_SWG_H_

#include "utils/essentials.h"
#include "align/align_swg.h"
#include "align/align_swg_score.h"
#include "matches/align/match_alignment.h"

/*
 * Constants
 */
#define SWG_SCORE_MIN (INT16_MIN)

/*
 * SWG Cell
 */
typedef struct {
  int32_t M; // Alignment matching/mismatching
  int32_t I; // Alignment ends with a gap in the reference (insertion)
  int32_t D; // Alignment ends with a gap in the read (deletion)
} swg_cell_t;

/*
 * Smith-waterman-gotoh - Init
 */
swg_cell_t** align_swg_allocate_table(
    const uint64_t num_columns,
    const uint64_t num_rows,
    mm_allocator_t* const mm_allocator);
void align_swg_init_table_banded(
    swg_cell_t** const dp,
    const uint64_t num_columns,
    const uint64_t num_rows,
    const uint64_t column_start_band,
    const uint64_t band_low_offset,
    const bool begin_free,
    const int32_t single_gap,
    const int32_t gap_extension);
void align_swg_allocate__init_column_banded(
    swg_cell_t** const dp,
    const uint64_t column_idx,
    const uint64_t column_start_band,
    const uint64_t band_low_offset,
    const uint64_t band_high_offset,
    const int32_t single_gap,
    const int32_t gap_extension,
    mm_allocator_t* const mm_allocator);

/*
 * Smith-waterman-gotoh Base (i.e. no-optimizations)
 */
void align_swg_base(
    const uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    const swg_penalties_t* swg_penalties,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator);

/*
 * Smith-Waterman-Gotoh - Main procedure (Dispatcher)
 */
void align_swg(
    const uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    const swg_penalties_t* const swg_penalties,
    const uint64_t max_bandwidth,
    const bool begin_free,
    const bool end_free,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator);

/*
 * Smith-Waterman-Gotoh - BackTrace
 */
void align_swg_traceback(
    /* Input key/Text */
    const uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    const bool reverse_strings,
    /* DP Computed */
    swg_cell_t** const dp,
    const int32_t max_score,
    const uint64_t max_score_row,
    const uint64_t max_score_column,
    /* SWG Scoring */
    const int32_t single_gap,
    const swg_matching_score_t* const matching_score,
    /* Alignment Parameters */
    const bool begin_free,
    const bool left_gap_alignment,
    /* Match-Alignment Results */
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector);

/*
 * Display
 */
void align_swg_print_table(
    swg_cell_t** const dp,
    const uint64_t num_columns,
    const uint64_t num_rows);
void align_swg_print_input(
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length);

#endif /* ALIGN_SWG_H_ */
