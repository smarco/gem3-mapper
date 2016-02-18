/*
 * PROJECT: GEMMapper
 * FILE: align_swg.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ALIGN_SWG_H_
#define ALIGN_SWG_H_

#include "utils/essentials.h"
#include "matches/match_align_dto.h"

/*
 * Constants
 */
#define SWG_SCORE_MIN (INT16_MIN)

/*
 * SWG Cell
 */
typedef struct {
  int32_t M; // Alignment matching/mismatching
  int32_t I; // Alignment ends with a gap in the reference (insertion)
  int32_t D; // Alignment ends with a gap in the read (deletion)
} swg_cell_t;

/*
 * Smith-waterman-gotoh Base (ie. no-optimizations)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @match_alignment->match_position (Adjusted)
 */
void align_swg_base(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack);
/*
 * Smith-Waterman-Gotoh - Main procedure (Dispatcher)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @match_alignment->match_position (Adjusted)
 *   @match_alignment->cigar_length (Cumulative)
 */
void align_swg(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    const bool begin_free,const bool end_free,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack);

/*
 * Display
 */
void align_swg_print_table(
    swg_cell_t** const dp,
    const uint64_t num_columns,
    const uint64_t num_rows);
void align_swg_print_input(
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length);

#endif /* ALIGN_SWG_H_ */
