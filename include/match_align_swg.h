/*
 * PROJECT: GEMMapper
 * FILE: match_align_swg.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGN_SWG_H_
#define MATCH_ALIGN_SWG_H_

#include "essentials.h"
#include "match_align_dto.h"
#include "match_alignment.h"
#include "match_scaffold.h"
#include "matches.h"

/*
 * SWG Align Matching Region (@region_matching)
 *   @align_input->key
 *   @align_input->text
 *   @align_parameters->swg_penalties
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @match_alignment->cigar_length (Cumulative)
 */
void match_align_swg_add_region_matching(
    const region_matching_t* const region_matching,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,mm_stack_t* const mm_stack);
void match_align_swg_add_gap(
    matches_t* const matches,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_length,const uint64_t text_chunk_length,
    const bool trim);

/*
 * SWG Align region
 */
bool match_align_swg_middle_region(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,const uint64_t text_chunk_length,
    mm_stack_t* const mm_stack);
bool match_align_swg_begin_region(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,const uint64_t text_chunk_length,
    mm_stack_t* const mm_stack);
bool match_align_swg_end_region(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,const uint64_t text_chunk_length,
    mm_stack_t* const mm_stack);

/*
 * SWG Post-Alignment (Normalize cigar, compute metrics, filter bad-alignments)
 */
void match_align_swg_post_alignment(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters);

/*
 * SWG Chained Alignment
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @match_scaffold->scaffold_regions
 *   @match_scaffold->num_scaffold_regions
 *   @match_alignment->cigar_length (Cumulative)
 *   @match_alignment->effective_length (Cumulative)
 *   @match_alignment->score (Cumulative)
 */
void match_align_swg_chain_scaffold(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack);

#endif /* MATCH_ALIGN_SWG_H_ */
