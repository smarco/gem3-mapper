/*
 * PROJECT: GEMMapper
 * FILE: match_align_swg_local.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/match_align_swg_local.h"

/*
 * SWG-Local MAX
 */
typedef struct {
  uint64_t local_begin_key_pos;
  uint64_t local_end_key_pos;
  uint64_t local_begin_text_pos;
  uint64_t local_end_text_pos;
  uint64_t local_begin_cigar;
  uint64_t local_end_cigar;
  int64_t local_score;
  uint64_t local_identity;
} match_align_swg_local_max_t;

void match_align_swg_local_alignment_add_gap(
    cigar_element_t** const local_cigar_buffer,
    const uint64_t begin_key_pos,
    const uint64_t begin_text_pos,
    const uint64_t end_key_pos,
    const uint64_t end_text_pos) {
  // Compute sentinels differences
  const uint64_t key_pos_diff = end_key_pos - begin_key_pos;
  const uint64_t text_pos_diff = end_text_pos - begin_text_pos;
  // Delete the read chunk
  if (key_pos_diff > 0) {
    (*local_cigar_buffer)->type = cigar_del;
    (*local_cigar_buffer)->length = key_pos_diff;
    (*local_cigar_buffer)->attributes = cigar_attr_none;
    ++(*local_cigar_buffer);
  }
  // Insert the text chunk
  if (text_pos_diff > 0) {
    (*local_cigar_buffer)->type = cigar_ins;
    (*local_cigar_buffer)->length = text_pos_diff;
    (*local_cigar_buffer)->attributes = cigar_attr_none;
    ++(*local_cigar_buffer);
  }
}
void match_align_swg_local_alignment_add_local_chunk(
    cigar_element_t* const global_cigar_buffer,
    cigar_element_t** const local_cigar_buffer,
    match_align_swg_local_max_t* const local_max) {
  // Add local-CIGAR
  uint64_t i;
  for (i=local_max->local_begin_cigar;i<=local_max->local_end_cigar;++i) {
     **local_cigar_buffer = global_cigar_buffer[i];
     ++(*local_cigar_buffer);
  }
}
/*
 * SWG-Local Alignment
 */
void match_align_swg_local_alignment(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters) {
  // Parameters
  const swg_penalties_t* const swg_penalties = align_parameters->swg_penalties;
  const int64_t local_min_swg_threshold = align_parameters->local_min_swg_threshold;
  const uint64_t local_min_identity = align_parameters->local_min_identity;
  vector_t* const cigar_vector = matches->cigar_vector;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Prepare CIGARs
  const uint64_t global_cigar_length = match_alignment->cigar_length;
  vector_reserve_additional(cigar_vector,global_cigar_length);
  cigar_element_t* const global_cigar_buffer = vector_get_elm(cigar_vector,match_alignment->cigar_offset,cigar_element_t);
  cigar_element_t* local_cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t);
  const cigar_element_t* const local_cigar_buffer_base = local_cigar_buffer;
  // Local Max
  match_align_swg_local_max_t local_max = {
    .local_begin_key_pos = 0,
    .local_end_key_pos = 0,
    .local_begin_text_pos = 0,
    .local_end_text_pos = 0,
    .local_begin_cigar = 0,
    .local_end_cigar = 0,
    .local_score = 0,
    .local_identity = 0,
  };
  // Traverse all CIGAR elements
  uint64_t key_pos = 0, text_pos = 0;
  uint64_t global_key_pos = 0, global_text_pos = 0;
  uint64_t local_identity = 0, i;
  int64_t local_score = 0, global_score = 0;
  for (i=0;i<global_cigar_length;++i) {
    // Add CIGAR operation
    switch (global_cigar_buffer[i].type) {
      case cigar_match: {
        const int32_t match_length = global_cigar_buffer[i].length;
        local_score += align_swg_score_match(swg_penalties,match_length);
        local_identity += match_length;
        text_pos += match_length;
        key_pos += match_length;
        break;
      }
      case cigar_mismatch:
        local_score += align_swg_score_mismatch(swg_penalties);
        ++text_pos;
        ++key_pos;
        break;
      case cigar_ins: {
        const int32_t indel_length = global_cigar_buffer[i].length;
        local_score += align_swg_score_insertion(swg_penalties,indel_length);
        text_pos += indel_length;
        break;
      }
      case cigar_del: {
        const int32_t indel_length = global_cigar_buffer[i].length;
        local_score += align_swg_score_deletion(swg_penalties,indel_length);
        key_pos += indel_length;
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
    ++i;
    // Check local score
    if (local_score > local_max.local_score) {
      // Max local score
      local_max.local_score = local_score;
      local_max.local_identity = local_identity;
      local_max.local_end_key_pos = key_pos;
      local_max.local_end_text_pos = text_pos;
      local_max.local_end_cigar = i;
    } else if (local_score < 0) {
      // Store local-max alignment chunk
      if (local_max.local_score >= local_min_swg_threshold &&
          local_max.local_identity >= local_min_identity) {
        match_align_swg_local_alignment_add_gap(&local_cigar_buffer,
            global_key_pos,global_text_pos,local_max.local_begin_key_pos,local_max.local_begin_text_pos);
        match_align_swg_local_alignment_add_local_chunk(global_cigar_buffer,&local_cigar_buffer,&local_max);
        global_key_pos = local_max.local_end_key_pos;
        global_text_pos = local_max.local_end_text_pos;
        global_score += local_max.local_score;
      }
      // Reset
      local_score = 0;
      local_identity = 0;
      local_max.local_score = 0;
      local_max.local_identity = 0;
      local_max.local_begin_key_pos = key_pos;
      local_max.local_begin_text_pos = text_pos;
      local_max.local_begin_cigar = i;
    }
  }
  // Check local score
  if (local_max.local_score >= local_min_swg_threshold &&
      local_max.local_identity >= local_min_identity) {
    match_align_swg_local_alignment_add_gap(&local_cigar_buffer,
        global_key_pos,global_text_pos,local_max.local_begin_key_pos,local_max.local_begin_text_pos);
    match_align_swg_local_alignment_add_local_chunk(global_cigar_buffer,&local_cigar_buffer,&local_max);
    global_key_pos = local_max.local_end_key_pos;
    global_text_pos = local_max.local_end_text_pos;
    global_score += local_max.local_score;
  }
  // Add final trim (if needed)
  if (global_key_pos != align_input->key_length) {
    match_align_swg_local_alignment_add_gap(
        &local_cigar_buffer,global_key_pos,0,align_input->key_length,0);
  }
  // Set local-CIGAR buffer used
  const uint64_t local_cigar_used = local_cigar_buffer-local_cigar_buffer_base;
  vector_add_used(cigar_vector,local_cigar_used);
  // Set SWG-Score
  match_trace->swg_score = global_score;
}
