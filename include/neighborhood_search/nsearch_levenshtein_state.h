/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein_state.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NSEARCH_LEVENSHTEIN_STATE_H_
#define NSEARCH_LEVENSHTEIN_STATE_H_

#include "utils/essentials.h"
#include "neighborhood_search/dp_matrix.h"

/*
 * Levenshtein Search State
 */
typedef struct {
  // DP-Search State
  dp_matrix_t dp_matrix;
  // Current key chunk
  uint64_t key_begin;
  uint64_t key_end;
  // Current text (search string)
  uint8_t* local_text;
  uint64_t local_text_length;
  // Global text (search string)
  uint8_t* global_text;
  uint64_t global_text_length;
} nsearch_levenshtein_state_t;

/*
 * Setup
 */
void nsearch_levenshtein_state_init(
    nsearch_levenshtein_state_t* const nsearch_state,
    const uint64_t num_columns,const uint64_t column_length,
    mm_stack_t* const mm_stack);

/*
 * Prepare DP
 */
void nsearch_levenshtein_state_prepare_full(
    nsearch_levenshtein_state_t* const nsearch_state,const uint64_t key_begin,
    const uint64_t key_end,const uint64_t max_error);
void nsearch_levenshtein_state_prepare_supercondensed(
    nsearch_levenshtein_state_t* const nsearch_state,const uint64_t max_error);

/*
 * Compute DP
 */
void nsearch_levenshtein_state_compute_chararacter(
    nsearch_levenshtein_state_t* const nsearch_state,const bool forward_search,
    const uint8_t* const key,const uint64_t key_begin,const uint64_t key_end,
    const uint64_t text_offset,const uint8_t text_char_enc,
    uint64_t* const min_val,uint64_t* const align_distance);
void nsearch_levenshtein_state_compute_sequence(
    nsearch_levenshtein_state_t* const nsearch_state,
    nsearch_levenshtein_state_t* const next_nsearch_state,
    const bool forward_search);

/*
 * Display
 */
void nsearch_levenshtein_state_print(
    FILE* const stream,nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,const uint8_t* const key);
void nsearch_levenshtein_state_print_search_text(
    FILE* const stream,nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search);
void nsearch_levenshtein_state_print_local_text(
    FILE* const stream,nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search);

#endif /* NSEARCH_LEVENSHTEIN_STATE_H_ */
