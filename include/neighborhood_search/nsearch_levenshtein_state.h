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
    mm_stack_t* const mm_stack);

/*
 * Prepare
 */
void nsearch_levenshtein_state_prepare(
    nsearch_levenshtein_state_t* const nsearch_state,
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
    uint64_t* const align_distance);

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
    uint64_t* const min_align_distance_column);
void nsearch_levenshtein_state_compute_chararacter_banded(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t text_position,
    const uint8_t text_char_enc,
    const uint64_t max_error,
    uint64_t* const min_val,
    uint64_t* const align_distance);

#endif /* NSEARCH_LEVENSHTEIN_STATE_H_ */
