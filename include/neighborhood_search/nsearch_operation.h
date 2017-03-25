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

#ifndef NSEARCH_OPERATION_H_
#define NSEARCH_OPERATION_H_

#include "neighborhood_search/nsearch_levenshtein_state.h"

/*
 * Neighborhood Search Operation
 */
typedef enum { direction_forward, direction_backward } search_direction_t;
typedef struct {
  // Search direction
  search_direction_t search_direction;
  // Error
  uint64_t min_local_error;  // Minimum error to the local partition/segment
  uint64_t min_global_error; // Minimum error accumulated (including current partition/segment)
  uint64_t max_global_error; // Maximum error accumulated (including current partition/segment)
  // Local Key
  uint64_t local_key_begin;
  uint64_t local_key_end;
  // Global Key
  uint64_t global_key_begin;
  uint64_t global_key_end;
  // Text (Operation Search String)
  uint8_t* text;
  uint64_t text_position;
  // Search State
  nsearch_levenshtein_state_t nsearch_state;
} nsearch_operation_t;

/*
 * Setup
 */
void nsearch_operation_init(
    nsearch_operation_t* const nsearch_operation,
    const uint64_t max_key_length,
    const uint64_t max_text_length,
    mm_allocator_t* const mm_allocator);

/*
 * Prepare Operation Chained
 */
bool nsearch_operation_chained_prepare(
    nsearch_operation_t* const current_nsearch_operation,
    nsearch_operation_t* const next_nsearch_operation,
    uint8_t* const key,
    const uint64_t key_length,
    const bool reverse_sequence,
    mm_allocator_t* const mm_allocator);

/*
 * Utils
 */
bool nsearch_operation_state_text_eq(
    nsearch_operation_t* const nsearch_operation,
    char* const text,
    mm_allocator_t* const mm_allocator);

/*
 * Display
 */
void nsearch_operation_print(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation);
void nsearch_operation_state_print(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation,
    const uint8_t* const key);
void nsearch_operation_state_print_global_text(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation);
void nsearch_operation_state_print_local_text(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation);

#endif /* NSEARCH_OPERATION_H_ */
