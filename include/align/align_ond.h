/*
 * PROJECT: GEMMapper
 * FILE: align_ond.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ALIGN_OND_H_
#define ALIGN_OND_H_

#include "utils/essentials.h"
#include "matches/match_align_dto.h"

typedef struct {
  int32_t** contour;
  uint64_t lcs_distance;
  uint64_t match_end_column;
} align_ond_contours_t;

/*
 * O(ND) Compute LCS
 */
void align_ond_compute_lcs_distance(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    uint64_t* const lcs_distance,
    uint64_t* const match_end_column,
    mm_stack_t* const mm_stack);

/*
 * O(ND) Align
 */
void align_ond_compute_contours(
    match_align_input_t* const align_input,
    const int32_t max_distance,
    align_ond_contours_t* const align_ond_contours,
    mm_stack_t* const mm_stack);
void align_ond_backtrace_contours(
    match_align_input_t* const align_input,
    align_ond_contours_t* const align_ond_contours,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector);
void align_ond_match(
    match_align_input_t* const align_input,
    const int32_t max_distance,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack);

/*
 * Display
 */
void align_ond_print_contour(
    FILE* const stream,
    const int32_t* const contour,
    const int32_t begin_contour,
    const int32_t end_contour,
    const int32_t distance);

#endif /* ALIGN_OND_H_ */
