/*
 * PROJECT: GEMMapper
 * FILE: match_align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_align.h"
#include "match_align_swg.h"
#include "match_align_dto.h"
#include "matches_cigar.h"
#include "align.h"
#include "align_swg.h"
#include "align_bpm.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Curate Alignment
 */
bool match_align_curate_cigar_trim(
    match_align_parameters_t* const align_parameters,const cigar_element_t* const cigar_element,
    uint64_t* const trim_length,uint64_t* const match_position) {
  switch (cigar_element->type) {
    case cigar_match:
      // Small Match
      if (cigar_element->length < align_parameters->cigar_curation_min_end_context) {
        *trim_length += cigar_element->length;
        if (match_position!=NULL) *match_position += cigar_element->length;
        return true;
      }
    break;
    case cigar_mismatch:
      ++(*trim_length);
      if (match_position!=NULL) ++(*match_position);
      return true;
      break;
    case cigar_del:
      *trim_length += cigar_element->length;
      return true;
      break;
    case cigar_ins:
      if (match_position!=NULL) *match_position += cigar_element->length;
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}
void match_align_curate_cigar(
    match_trace_t* const match_trace,vector_t* const cigar_vector,
    match_align_parameters_t* const align_parameters) {
  // Parameters
  const bool left_gap_alignment = align_parameters->left_gap_alignment;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  const uint64_t cigar_length = match_alignment->cigar_length;
  cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,match_alignment->cigar_offset,cigar_element_t);
  uint64_t curated_cigar_length = 0, indel_length = 0, i = 0, j = 0;
  // Trim the beginning of the read
  while (j < cigar_length) {
    cigar_element_t* const cigar_element = cigar_buffer + j;
    if (!match_align_curate_cigar_trim(align_parameters,cigar_element,&indel_length,&match_alignment->match_position)) {
      break;
    }
    // Trim
    cigar_element->type = cigar_del;
    cigar_element->attributes = cigar_attr_trim;
    cigar_element->length = indel_length;
    i = j++;
  }
  // Traverse all CIGAR elements
  while (i < cigar_length) {
    cigar_element_t* cigar_element = cigar_buffer + i;
    switch (cigar_element->type) {
      case cigar_mismatch:
        cigar_buffer[curated_cigar_length++] = *cigar_element;
        ++i;
        break;
      case cigar_match:
        if (i+1<cigar_length && cigar_buffer[i+1].type == cigar_match) {
          cigar_buffer[i+1].length += cigar_element->length;
        } else {
          cigar_buffer[curated_cigar_length++] = *cigar_element;
        }
        ++i;
        break;
      case cigar_del:
      case cigar_ins: {
        cigar_element_t accum_del = { .type = cigar_del, .attributes = cigar_attr_none, .length = 0 };
        cigar_element_t accum_ins = { .type = cigar_ins, .attributes = cigar_attr_none, .length = 0 };
        if (cigar_element->type==cigar_del && cigar_element->attributes==cigar_attr_trim) {
          accum_del.attributes = cigar_attr_trim;
        }
        // Compact all deletions/insertions in a row
        while (i<cigar_length) {
          cigar_element = cigar_buffer + i;
          if (cigar_element->type == cigar_del) {
            accum_del.length += cigar_element->length;
            ++i;
          } else if (cigar_element->type == cigar_ins) {
            accum_ins.length += cigar_element->length;
            ++i;
          } else {
            break;
          }
        }
        // Copy curated Ins/Del
        if (left_gap_alignment) {
          if (accum_del.length > 0) cigar_buffer[curated_cigar_length++] = accum_del;
          if (accum_ins.length > 0) cigar_buffer[curated_cigar_length++] = accum_ins;
        } else {
          if (accum_ins.length > 0) cigar_buffer[curated_cigar_length++] = accum_ins;
          if (accum_del.length > 0) cigar_buffer[curated_cigar_length++] = accum_del;
        }
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Handle last CIGAR element(s)
  if (curated_cigar_length > 0) {
    // Check match or deletion at the end
    cigar_element_t* last_cigar_element = cigar_buffer + (curated_cigar_length-1);
    uint64_t indel_length = 0;
    if (match_align_curate_cigar_trim(align_parameters,last_cigar_element,&indel_length,NULL)) {
      // Chain all the mismatches & deletions at the end
      while (last_cigar_element > cigar_buffer) {
        --last_cigar_element;
        if (!match_align_curate_cigar_trim(align_parameters,last_cigar_element,&indel_length,NULL)) break;
        --curated_cigar_length;
      }
      // Merge all of them
      last_cigar_element = cigar_buffer + (curated_cigar_length-1);
      last_cigar_element->type = cigar_del;
      last_cigar_element->attributes = cigar_attr_trim;
      last_cigar_element->length = indel_length;
      if (indel_length==0) --curated_cigar_length;
    }
  }
  // Set curated-CIGAR length
  match_alignment->cigar_length = curated_cigar_length;
}
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
void match_align_exact(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters) {
  PROFILE_START(GP_MATCHES_ALIGN_EXACT,PROFILE_LEVEL);
  // Parameters
  const uint64_t key_length = align_input->key_length;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Configure match-trace
  match_trace->text_trace_offset = align_input->text_trace_offset;
  match_trace->text = (match_alignment->score==0) ?
      align_input->text + align_input->text_offset_begin : align_input->text;
  match_trace->text_length = align_input->key_length;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  match_trace->distance = match_alignment->score;
  match_trace->edit_distance = match_alignment->score;
  match_trace->swg_score = align_swg_score_match(align_parameters->swg_penalties,(int32_t)key_length);
  // Insert exact-match CIGAR
  match_alignment->match_offset = (match_alignment->score==0) ? align_input->text_offset_begin : 0;
  match_alignment->match_position = align_input->text_position + match_alignment->match_offset; // Adjust position
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  match_alignment->effective_length = key_length;
  matches_cigar_vector_append_match(matches->cigar_vector,&match_alignment->cigar_length,key_length,cigar_attr_none);
  PROFILE_STOP(GP_MATCHES_ALIGN_EXACT,PROFILE_LEVEL);
}
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
void match_align_hamming(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters) {
  PROFILE_START(GP_MATCHES_ALIGN_HAMMING,PROFILE_LEVEL);
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  const uint64_t text_offset_begin = align_input->text_offset_begin;
  const uint64_t text_offset_end = align_input->text_offset_end;
  const bool* const allowed_enc = align_parameters->allowed_enc;
  // Configure match-trace
  match_trace->text_trace_offset = align_input->text_trace_offset;
  match_trace->text = text + text_offset_begin;
  match_trace->text_length = key_length;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment->match_position = align_input->text_position;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  match_alignment->match_offset = 0;
  // Hamming Alignment
  uint64_t i, mismatches;
  for (i=text_offset_begin,mismatches=0;i<text_offset_end;++i) {
    const uint8_t candidate_enc = text[i];
    // Check Mismatch
    if (!allowed_enc[candidate_enc] || candidate_enc != key[i]) {
      ++mismatches;
      matches_cigar_vector_append_mismatch(matches->cigar_vector,&match_alignment->cigar_length,candidate_enc,cigar_attr_none);
    } else {
      matches_cigar_vector_append_match(matches->cigar_vector,&match_alignment->cigar_length,1,cigar_attr_none);
    }
  }
  match_alignment->effective_length = key_length;
  match_trace->distance = mismatches;
  match_trace->edit_distance = mismatches;
  match_trace->swg_score = align_swg_score_cigar(align_parameters->swg_penalties,
      matches->cigar_vector,match_alignment->cigar_offset,match_alignment->cigar_length);
  PROFILE_STOP(GP_MATCHES_ALIGN_HAMMING,PROFILE_LEVEL);
}
/*
 * Levenshtein
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text_trace_offset
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->max_error
 *   @align_parameters->left_gap_alignment
 */
void match_align_levenshtein(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    mm_stack_t* const mm_stack) {
  PROFILE_START(GP_MATCHES_ALIGN_LEVENSHTEIN,PROFILE_LEVEL);
  // Configure match-trace
  match_trace->text_trace_offset = align_input->text_trace_offset;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  // Levenshtein Align
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  align_input->text += align_input->text_offset_begin;
  align_input->text_length = align_input->text_offset_end - align_input->text_offset_begin;
  const uint64_t match_position = align_input->text_position + align_input->text_offset_begin;
  match_alignment->match_position = match_position;
  align_bpm_match(align_input,align_parameters->max_error,
      align_parameters->left_gap_alignment,match_alignment,matches->cigar_vector,mm_stack);
  match_trace->distance = match_alignment->score;
  match_trace->edit_distance = match_alignment->score;
  match_trace->swg_score = align_swg_score_cigar(align_parameters->swg_penalties,
      matches->cigar_vector,match_alignment->cigar_offset,match_alignment->cigar_length);
  // Store matching text
  match_alignment->match_offset = match_alignment->match_position - match_position;
  match_trace->text = align_input->text + match_alignment->match_offset;
  match_trace->text_length = match_alignment->effective_length;
  PROFILE_STOP(GP_MATCHES_ALIGN_LEVENSHTEIN,PROFILE_LEVEL);
}
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
 *   @align_parameters->swg_threshold
 *   @align_parameters->min_identity
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @align_parameters->cigar_curation
 *   @align_parameters->cigar_curation_min_end_context
 *   @match_scaffold->scaffold_regions
 *   @match_scaffold->num_scaffold_regions
 */
void match_align_smith_waterman_gotoh(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
  // Configure match-trace
  match_trace->text_trace_offset = align_input->text_trace_offset;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  // Scaffold the alignment
  if (align_parameters->scaffolding) {
    match_scaffold_alignment(matches,align_input,align_parameters,match_scaffold,mm_stack);
//    // Check minimum identity
//    if (align_parameters->min_identity + match_scaffold->match_alignment.score > align_input->key_length) {
//      match_scaffold->num_scaffold_regions = 0;
//      match_scaffold->scaffolding_coverage = 0;
//    }
  }
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
#ifdef GEM_DEBUG
  match_trace->match_scaffold = (match_scaffold!=NULL && num_scaffold_regions > 0) ? match_scaffold : NULL;
#endif
  // Configure match-alignment
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  const uint64_t base_match_position = align_input->text_position;
  match_alignment->match_position = base_match_position;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  // Check the number of matching regions
  if (match_scaffold!=NULL && num_scaffold_regions > 0) {
    // Chain matching regions and align gaps (SWG)
    match_align_swg_chain_scaffold(matches,match_trace,align_input,align_parameters,match_scaffold,mm_stack);
  } else {
    // Full SWG
    match_align_input_t align_chunk_input;
    align_chunk_input.key = align_input->key;
    align_chunk_input.key_length = align_input->key_length;
    align_chunk_input.text = align_input->text + align_input->text_offset_begin;
    align_chunk_input.text_length = align_input->text_offset_end - align_input->text_offset_begin;
    match_alignment->match_position += align_input->text_offset_begin;
    align_swg(&align_chunk_input,align_parameters,true,true,
        &match_trace->match_alignment,matches->cigar_vector,mm_stack);
  }
  // Post alignment checks & setup
  match_align_swg_post_alignment(matches,match_trace,align_input,align_parameters,false);
  // Adjust text
  match_alignment->match_offset = match_alignment->match_position - base_match_position;
  match_trace->text = align_input->text + match_alignment->match_offset;
  match_trace->text_length = match_alignment->effective_length;
  PROFILE_STOP(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
}

