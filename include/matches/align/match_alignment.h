/*
 * PROJECT: GEMMapper
 * FILE: match_alignment.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGNMENT_H_
#define MATCH_ALIGNMENT_H_

#include "utils/essentials.h"

/*
 * Match Alignment
 */
typedef struct {
  uint64_t match_text_offset; // Match text offset (wrt beginning of text-candidate)
  uint64_t match_position;    // Match position
  uint64_t cigar_offset;      // CIGAR offset in buffer
  uint64_t cigar_length;      // CIGAR length
  int64_t effective_length;   // Match effective length
  int32_t score;              // Score assigned by the aligner
} match_alignment_t;
/*
 * Alignment Model
 */
typedef enum {
  match_alignment_model_none,
  match_alignment_model_hamming,
  match_alignment_model_levenshtein,
  match_alignment_model_gap_affine
} match_alignment_model_t;

/*
 * Display
 */
void match_alignment_print_pretty(
    FILE* const stream,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    mm_stack_t* const mm_stack);

#endif /* MATCH_ALIGNMENT_H_ */
