/*
 * PROJECT: GEMMapper
 * FILE: nsearch_operation.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
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
  uint64_t max_local_error;
  uint64_t min_local_error;
  uint64_t max_global_error;
  uint64_t min_global_error;
  uint64_t max_text_length;
  // Local Key
  uint64_t local_key_begin;
  uint64_t local_key_end;
  // Global Key
  uint64_t global_key_begin;
  uint64_t global_key_end;
  // Local Text (Local search string)
  uint8_t* text;               // Current local text
  uint64_t text_position;      // Current text position
  uint64_t context_length;     // Context length
  // Global Text (Global search string)
  uint8_t* global_text;
  uint64_t global_text_length;
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
    mm_stack_t* const mm_stack);

/*
 * Prepare Sequence
 */
void nsearch_operation_chained_prepare_sequence(
    nsearch_operation_t* const current_operation,
    nsearch_operation_t* const next_operation);
void nsearch_operation_chained_prepare_context(
    nsearch_operation_t* const current_operation,
    nsearch_operation_t* const next_operation,
    const uint64_t first_active_column);
void nsearch_operation_chained_prepare_reverse_sequence(
    nsearch_operation_t* const current_operation,
    nsearch_operation_t* const next_operation);
void nsearch_operation_chained_prepare_reverse_context(
    nsearch_operation_t* const next_operation,
    const uint64_t first_active_column,
    const uint64_t last_active_column);

/*
 * Display
 */
void nsearch_operation_print(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation);
void nsearch_operation_state_print(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const uint8_t* const key);
void nsearch_operation_state_print_global_text(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search);
void nsearch_operation_state_print_local_text(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search);

#endif /* NSEARCH_OPERATION_H_ */
