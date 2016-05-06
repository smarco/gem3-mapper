/*
 * PROJECT: GEMMapper
 * FILE: align_bpm.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: BPM-alignment (BitParalellMyers, Bit-compressed Alignment)
 *   Myers' Fast Bit-Vector algorithm to compute levenshtein alignment (CIGAR)
 */

#ifndef ALIGN_BPM_H_
#define ALIGN_BPM_H_

#include "utils/essentials.h"
#include "matches/match_alignment.h"
#include "matches/match_align_dto.h"

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
