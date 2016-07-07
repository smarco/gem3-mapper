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
    mm_stack_t* const mm_stack);

/*
 * Prepare Operation Chained
 */
bool nsearch_operation_chained_prepare(
    nsearch_operation_t* const current_nsearch_operation,
    nsearch_operation_t* const next_nsearch_operation,
    uint8_t* const key,
    const uint64_t key_length,
    const bool reverse_sequence);

/*
 * Utils
 */
bool nsearch_operation_state_text_eq(
    nsearch_operation_t* const nsearch_operation,
    char* const text,
    mm_stack_t* const mm_stack);

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
