/*
 * PROJECT: GEMMapper
 * FILE: match_align_local.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_align_local.h"
#include "match_align_swg.h"
#include "align.h"
#include "align_swg.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Chained local-SWG Alignment
 */
GEM_INLINE void match_align_local_swg_chain_local_max(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t key_length = align_input->key_length;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  region_matching_t* const scaffold_regions = match_scaffold->scaffold_regions;
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  vector_t* const cigar_vector = matches->cigar_vector;
  // Auxiliary Variables
  const uint64_t match_text_base_position = match_alignment->match_position;
  uint64_t i, key_chunk_begin_offset, key_chunk_length;
  uint64_t text_chunk_begin_offset, text_chunk_length, match_begin_position;
  bool begin_region = true;
  // Keep CIGAR state
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t cigar_vector_used = vector_get_used(cigar_vector);
  // Chain matching regions
  for (i=0;i<num_scaffold_regions;++i) {
    const region_matching_t* const current_region_matching = scaffold_regions + i;
    key_chunk_begin_offset = current_region_matching->key_begin;
    key_chunk_length = current_region_matching->key_end - current_region_matching->key_begin;
    text_chunk_begin_offset = current_region_matching->text_begin;
    text_chunk_length = current_region_matching->text_end - text_chunk_begin_offset;
    if (begin_region) {
      // Trim the alignment (to first matching region)
      match_alignment->match_position = match_text_base_position + current_region_matching->text_begin;
      match_align_swg_add_gap(matches,match_alignment,current_region_matching->key_begin,0,true);
      // Align region
      begin_region = !match_align_swg_begin_region(matches,align_input,align_parameters,match_alignment,
          key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,true,mm_stack);
      if (begin_region) {
        // Restore the CIGAR
        match_alignment->cigar_length = cigar_length;
        vector_set_used(cigar_vector,cigar_vector_used);
      } else {
        match_begin_position = match_alignment->match_position;
      }
    } else {
      // Align region
      if (i+1 < num_scaffold_regions) {
        match_align_swg_middle_region(matches,align_input,align_parameters,match_alignment,
            key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,true,mm_stack);
      } else {
        match_align_swg_end_region(matches,align_input,align_parameters,match_alignment,
            key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,true,mm_stack);
      }
    }
    // Bride the gap
    if (!begin_region) {
      if (i+1 < num_scaffold_regions) {
        // Bride the gap
        const region_matching_t* const next_region_matching = current_region_matching + 1;
        key_chunk_begin_offset = current_region_matching->key_end;
        key_chunk_length = next_region_matching->key_begin - key_chunk_begin_offset;
        text_chunk_begin_offset = current_region_matching->text_end;
        text_chunk_length = next_region_matching->text_begin - text_chunk_begin_offset;
        match_align_swg_add_gap(matches,match_alignment,key_chunk_length,text_chunk_length,false);
      } else {
        // Trim the alignment (from last matching region)
        const region_matching_t* const last_region_matching = scaffold_regions + (num_scaffold_regions-1);
        key_chunk_length = key_length - last_region_matching->key_end;
        match_align_swg_add_gap(matches,match_alignment,key_chunk_length,0,true);
      }
    }
  }
  // Post-processing
  if (begin_region) {
    match_alignment->score = SWG_SCORE_MIN;
  } else {
    match_alignment->match_position = match_begin_position; // Restore match position
  }
}
/*
 * Local Smith-Waterman-Gotoh Alignment (Gap-affine)
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
 *   @align_parameters->swg_threshold
 *   @align_parameters->min_identity
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @align_parameters->scaffolding
 *   @align_parameters->scaffolding_min_coverage
 *   @align_parameters->scaffolding_matching_min_length
 *   @align_parameters->scaffolding_homopolymer_min_context
 *   @align_parameters->cigar_curation
 *   @align_parameters->cigar_curation_min_end_context
 *   @match_scaffold->scaffold_regions
 *   @match_scaffold->num_scaffold_regions
 */
GEM_INLINE void match_align_local_smith_waterman_gotoh(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_MATCHES_ALIGN_LOCAL_SWG,PROFILE_LEVEL);
  // Configure match-trace
  match_trace->trace_offset = align_input->text_trace_offset;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  // Scaffold the alignment
  if (!match_scaffold_smith_waterman_gotoh(matches,align_input,align_parameters,match_scaffold,mm_stack)) {
    match_trace->distance=ALIGN_DISTANCE_INF;
    PROFILE_STOP(GP_MATCHES_ALIGN_LOCAL_SWG,PROFILE_LEVEL);
    return;
  }
#ifdef GEM_DEBUG
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  match_trace->match_scaffold = (match_scaffold!=NULL && num_scaffold_regions > 0) ? match_scaffold : NULL;
#endif
  // Configure match-alignment
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  const uint64_t base_match_position = align_input->text_position;
  match_alignment->match_position = base_match_position;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  // Chain SWG-matching regions
  match_align_local_swg_chain_local_max(matches,match_trace,align_input,align_parameters,match_scaffold,mm_stack);
  // Post alignment checks & setup
  match_align_swg_post_alignment(matches,match_trace,align_input,align_parameters,true);
//  // Check correctness
//  align_check_match(stderr,align_input->key,align_input->key_length,
//      align_input->text+(match_alignment->match_position - base_match_position),match_alignment->effective_length,
//      matches->cigar_vector,match_trace->match_alignment.cigar_offset,match_trace->match_alignment.cigar_length,true);
  // Adjust text
  match_trace->text = align_input->text + (match_alignment->match_position - base_match_position);
  match_trace->text_length = match_alignment->effective_length;
  PROFILE_STOP(GP_MATCHES_ALIGN_LOCAL_SWG,PROFILE_LEVEL);
}
