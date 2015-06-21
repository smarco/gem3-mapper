/*
 * PROJECT: GEMMapper
 * FILE: match_align.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGN_H_
#define MATCH_ALIGN_H_

#include "essentials.h"
#include "matches.h"
#include "match_elements.h"
#include "match_scaffold.h"
#include "bpm_align.h"
#include "swg_align.h"

typedef enum {
  alignment_model_none,
  alignment_model_hamming,
  alignment_model_levenshtein,
  alignment_model_gap_affine
} alignment_model_t;

/*
 * Exact Alignment
 *   @align_input->key_length
 *   @align_input->text_position
 *   @align_input->text_trace_offset
 *   @align_input->text
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->swg_penalties
 *   @match_trace->match_alignment.score
 */
GEM_INLINE void match_align_exact(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters);

/*
 * Hamming Alignment (Only mismatches)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text_position
 *   @align_input->text_trace_offset
 *   @align_input->text
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->allowed_enc
 */
GEM_INLINE void match_align_hamming(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters);

/*
 * Levenshtein Alignment (edit distance)
 *   @align_input->key
 *   @align_input->bpm_pattern
 *   @align_input->text_trace_offset
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_input->text_length
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->max_error
 *   @align_parameters->left_gap_alignment
 */
GEM_INLINE void match_align_levenshtein(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    mm_stack_t* const mm_stack);

/*
 * Smith-Waterman-Gotoh Alignment (Gap-affine)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text_trace_offset
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->swg_penalties
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @match_scaffold->scaffold_regions
 *   @match_scaffold->num_scaffold_regions
 */
GEM_INLINE void match_align_smith_waterman_gotoh(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack);

#endif /* MATCH_ALIGN_H_ */
