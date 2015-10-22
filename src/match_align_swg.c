/*
 * PROJECT: GEMMapper
 * FILE: match_align_swg.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_align_swg.h"
#include "match_align_dto.h"
#include "match_align.h"
#include "align.h"
#include "align_swg.h"
#include "matches.h"
#include "matches_cigar.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

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
GEM_INLINE void match_align_swg_add_region_matching(
    const region_matching_t* const region_matching,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,mm_stack_t* const mm_stack) {
  // Parameters
  uint8_t* const key = align_input->key;
  const uint64_t key_matching_length = region_matching->key_end - region_matching->key_begin;
  uint8_t* const text = align_input->text;
  const uint64_t text_matching_length = region_matching->text_end - region_matching->text_begin;
  // Select matching-region type
  switch (region_matching->matching_type) {
    case region_matching_exact:
      matches_cigar_vector_append_match(cigar_vector,&match_alignment->cigar_length,key_matching_length,cigar_attr_none);
      match_alignment->score = align_swg_score_match(align_parameters->swg_penalties,key_matching_length);
      break;
    case region_matching_approximate:
      if (region_matching->cigar_length > 0) {
        // Copy the CIGAR from the matching region
        const uint64_t cigar_buffer_offset = region_matching->cigar_buffer_offset;
        uint64_t i;
        match_alignment->score = 0;
        for (i=0;i<region_matching->cigar_length;++i) {
          cigar_element_t* const scaffolding_elm = vector_get_elm(cigar_vector,cigar_buffer_offset+i,cigar_element_t);
          match_alignment->score += align_swg_score_cigar_element(align_parameters->swg_penalties,scaffolding_elm);
          matches_cigar_vector_append_cigar_element(cigar_vector,&match_alignment->cigar_length,scaffolding_elm);
        }
      } else {
        // Force (re)computing the CIGAR from the matching region (missing CIGAR or never computed)
        match_align_input_t align_chunk_input = {
            .key = key+region_matching->key_begin,
            .key_length = key_matching_length,
            .text = text+region_matching->text_begin,
            .text_length = text_matching_length,
        };
        match_align_parameters_t align_chunk_parameters = {
            .max_bandwidth = region_matching->error+1, // FIXME Include CIGAR when generating approximate regions
            .left_gap_alignment = align_parameters->left_gap_alignment,
            .allowed_enc = align_parameters->allowed_enc,
            .swg_penalties = align_parameters->swg_penalties,
        };
        align_swg(&align_chunk_input,&align_chunk_parameters,false,false,match_alignment,cigar_vector,mm_stack);
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * SWG Bride the gap
 */
GEM_INLINE void match_align_swg_add_gap(
    matches_t* const matches,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_length,const uint64_t text_chunk_length,const bool trim) {
  vector_t* const cigar_vector = matches->cigar_vector;
  if (trim) {
    // Trim the alignment
    if (key_chunk_length > 0) {
      matches_cigar_vector_append_deletion(cigar_vector,
          &match_alignment->cigar_length,key_chunk_length,cigar_attr_trim);
    }
  } else {
    // Delete the read chunk
    if (key_chunk_length > 0) {
      matches_cigar_vector_append_deletion(cigar_vector,
          &match_alignment->cigar_length,key_chunk_length,cigar_attr_none);
    }
    // Insert the text chunk
    if (text_chunk_length > 0) {
      matches_cigar_vector_append_insertion(cigar_vector,
          &match_alignment->cigar_length,text_chunk_length,cigar_attr_none);
    }
  }
}
/*
 * SWG Align region
 */
GEM_INLINE bool match_align_swg_region(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,const uint64_t text_chunk_length,
    const bool begin_free,const bool end_free,
    const bool force_swg_threshold,mm_stack_t* const mm_stack) {
  vector_t* const cigar_vector = matches->cigar_vector;
  // Align gap between regions matching
  match_align_input_t align_chunk_input;
  align_chunk_input.key = align_input->key+key_chunk_begin_offset;
  align_chunk_input.key_length = key_chunk_length;
  align_chunk_input.text = align_input->text+text_chunk_begin_offset;
  align_chunk_input.text_length = text_chunk_length;
  // Keep CIGAR state
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t cigar_vector_used = vector_get_used(cigar_vector);
  // SWG-Align
  align_swg(&align_chunk_input,align_parameters,begin_free,end_free,match_alignment,cigar_vector,mm_stack);
  // Check alignment result
  if (match_alignment->score == SWG_SCORE_MIN || (force_swg_threshold && match_alignment->score < align_parameters->swg_threshold)) {
    // Restore the CIGAR
    match_alignment->cigar_length = cigar_length;
    vector_set_used(cigar_vector,cigar_vector_used);
    // Bride the gap
    match_align_swg_add_gap(matches,match_alignment,key_chunk_length,text_chunk_length,begin_free||end_free);
    return false;
  } else {
    return true;
  }
}
GEM_INLINE bool match_align_swg_middle_region(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,const uint64_t text_chunk_length,
    const bool force_swg_threshold,mm_stack_t* const mm_stack) {
  return match_align_swg_region(matches,align_input,align_parameters,match_alignment,
      key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,
      false,false,force_swg_threshold,mm_stack);
}
GEM_INLINE bool match_align_swg_begin_region(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,const uint64_t text_chunk_length,
    const bool force_swg_threshold,mm_stack_t* const mm_stack) {
  return match_align_swg_region(matches,align_input,align_parameters,match_alignment,
      key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,
      true,false,force_swg_threshold,mm_stack);
}
GEM_INLINE bool match_align_swg_end_region(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,const uint64_t text_chunk_length,
    const bool force_swg_threshold,mm_stack_t* const mm_stack) {
  return match_align_swg_region(matches,align_input,align_parameters,match_alignment,
      key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,
      false,true,force_swg_threshold,mm_stack);
}
/*
 * SWG Post-Alignment (Curate cigar, compute metrics, filter bad-alignments)
 */
GEM_INLINE void match_align_swg_post_alignment(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    const bool local_alignment) {
  vector_t* const cigar_vector = matches->cigar_vector;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Check for bad alignments (discarded)
  if (match_alignment->score == SWG_SCORE_MIN) {
    match_trace->swg_score = SWG_SCORE_MIN;
    match_trace->distance = ALIGN_DISTANCE_INF;
    return;
  }
  // Curate alignment
  if (align_parameters->cigar_curation) {
    match_align_curate_cigar(match_trace,cigar_vector,align_parameters);
  }
  // Discard bad alignments & Compute edit distance
  const uint64_t cigar_offset = match_alignment->cigar_offset;
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t matching_bases = matches_cigar_compute_matching_bases(cigar_vector,cigar_offset,cigar_length);
  if (matching_bases < align_parameters->min_identity) {
    match_trace->swg_score = SWG_SCORE_MIN;
    match_trace->distance = ALIGN_DISTANCE_INF;
    return;
  }
  match_trace->edit_distance = matches_cigar_compute_edit_distance(cigar_vector,cigar_offset,cigar_length);
  // Compute score & discard alignments below threshold
//  match_trace->swg_score = (!local_alignment) ?
//      swg_score_cigar(align_parameters->swg_penalties,cigar_vector,cigar_offset,cigar_length) :
//      swg_score_cigar__excluding_clipping(align_parameters->swg_penalties,cigar_vector,cigar_offset,cigar_length);
  match_trace->swg_score = align_swg_score_cigar(align_parameters->swg_penalties,cigar_vector,cigar_offset,cigar_length);
  if (!local_alignment && match_trace->swg_score < align_parameters->swg_threshold) {
    match_trace->swg_score = SWG_SCORE_MIN;
    match_trace->distance = ALIGN_DISTANCE_INF;
    return;
  }
  // Compute distance + effective-length
  match_trace->distance = matches_cigar_compute_event_distance(cigar_vector,cigar_offset,cigar_length);
  match_alignment->effective_length = matches_cigar_effective_length(cigar_vector,cigar_offset,cigar_length);
}
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
GEM_INLINE void match_align_swg_chain_scaffold(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t key_length = align_input->key_length;
  const uint64_t text_length = align_input->text_length;
  const uint64_t max_bandwidth = align_parameters->max_bandwidth;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  region_matching_t* const scaffold_regions = match_scaffold->scaffold_regions;
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  vector_t* const cigar_vector = matches->cigar_vector;
  // Auxiliary Variables
  uint64_t key_chunk_begin_offset, key_chunk_length, i;
  uint64_t text_chunk_end_offset, text_chunk_begin_offset, text_chunk_length, match_begin_position;
  // Chain to first matching region
  const region_matching_t* const first_region_matching = scaffold_regions;
  key_chunk_length = first_region_matching->key_begin;
  key_chunk_begin_offset = 0;
  text_chunk_length = BOUNDED_ADDITION(key_chunk_length,max_bandwidth,first_region_matching->text_begin);
  text_chunk_begin_offset = first_region_matching->text_begin - text_chunk_length;
  match_alignment->match_position += text_chunk_begin_offset; // Offset match position
  match_begin_position = match_alignment->match_position; // Save match position
  const bool feasible_alignment = match_align_swg_begin_region(matches,
      align_input,align_parameters,match_alignment,key_chunk_begin_offset,
      key_chunk_length,text_chunk_begin_offset,text_chunk_length,false,mm_stack);
  match_begin_position = (feasible_alignment) ? match_alignment->match_position : match_begin_position + text_chunk_length;
  // Add first matching region
  match_align_swg_add_region_matching(first_region_matching,
      align_input,align_parameters,match_alignment,cigar_vector,mm_stack);
  // Chain matching regions
  for (i=1;i<num_scaffold_regions;++i) {
    // Align the gap
    const region_matching_t* const prev_region_matching = scaffold_regions + (i-1);
    const region_matching_t* const current_region_matching = scaffold_regions + i;
    const uint64_t key_chunk_begin_offset = prev_region_matching->key_end;
    const uint64_t key_chunk_length = current_region_matching->key_begin - key_chunk_begin_offset;
    const uint64_t text_chunk_begin_offset = prev_region_matching->text_end;
    const uint64_t text_chunk_length = current_region_matching->text_begin - text_chunk_begin_offset;
    match_align_swg_middle_region(matches,align_input,align_parameters,match_alignment,
        key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,false,mm_stack);
    // Add matching region
    match_align_swg_add_region_matching(current_region_matching,
        align_input,align_parameters,match_alignment,cigar_vector,mm_stack);
  }
  // Chain from last matching region
  const region_matching_t* const last_region_matching = scaffold_regions + (num_scaffold_regions-1);
  key_chunk_begin_offset = last_region_matching->key_end;
  key_chunk_length = key_length - last_region_matching->key_end;
  text_chunk_end_offset = BOUNDED_ADDITION(last_region_matching->text_end,key_chunk_length+max_bandwidth,text_length);
  text_chunk_begin_offset = last_region_matching->text_end;
  text_chunk_length = text_chunk_end_offset-last_region_matching->text_end;
  match_align_swg_end_region(matches,align_input,align_parameters,match_alignment,
      key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,false,mm_stack);
  // Post-processing
  match_alignment->match_position = match_begin_position; // Restore match position
}
