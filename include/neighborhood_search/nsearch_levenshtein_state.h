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

#ifndef NSEARCH_LEVENSHTEIN_STATE_H_
#define NSEARCH_LEVENSHTEIN_STATE_H_

#include "utils/essentials.h"
#include "neighborhood_search/dp_matrix.h"

/*
 * Constants
 */
#define NO_ACTIVE_COLUMN UINT64_MAX

/*
 * Levenshtein Search State
 */
typedef struct {
  dp_matrix_t dp_matrix;
  bool supercondensed;
} nsearch_levenshtein_state_t;

/*
 * Setup
 */
void nsearch_levenshtein_state_init(
    nsearch_levenshtein_state_t* const nsearch_levenshtein_state,
    const uint64_t num_rows,
    const uint64_t num_columns,
    mm_allocator_t* const mm_allocator);

/*
 * Prepare
 */
void nsearch_levenshtein_state_prepare(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool supercondensed);
void nsearch_levenshtein_state_prepare_column(
    dp_column_t* const column,
    const uint64_t column_position,
    const bool supercondensed);

/*
 * Accessors
 */
uint64_t nsearch_levenshtein_state_get_global_align_distance(
    nsearch_levenshtein_state_t* const nsearch_state,
    const uint64_t key_length,
    const uint64_t text_length,
    const uint64_t max_error);
uint64_t nsearch_levenshtein_state_get_local_align_distance(
    nsearch_levenshtein_state_t* const nsearch_state,
    const uint64_t local_key_length,
    const uint64_t global_key_length,
    const uint64_t global_text_length,
    const uint64_t max_error);

/*
 * Compute DP
 */
void nsearch_levenshtein_state_compute_chararacter(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t text_position,
    const uint8_t text_char_enc,
    const uint64_t max_error,
    uint64_t* const min_val,
    uint64_t* const align_distance,
    mm_allocator_t* const mm_allocator);

/*
 * Compute DP-Banded
 */
void nsearch_levenshtein_state_compute_text_banded(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error,
    uint64_t* const min_align_distance,
    uint64_t* const min_align_distance_column,
    mm_allocator_t* const mm_allocator);
void nsearch_levenshtein_state_compute_chararacter_banded(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t text_position,
    const uint8_t text_char_enc,
    const uint64_t max_error,
    uint64_t* const min_val,
    uint64_t* const align_distance,
    mm_allocator_t* const mm_allocator);

#endif /* NSEARCH_LEVENSHTEIN_STATE_H_ */
