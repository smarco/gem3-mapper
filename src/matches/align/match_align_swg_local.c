/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "matches/align/match_align_swg_local.h"

/*
 * SWG-Local
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
} match_align_swg_local_t;

/*
 * SWG-Local Compose Alignment
 */
void match_align_swg_local_alignment_add_match_compose_cigar(
    matches_t* const matches,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_align_swg_local_t* const local_alignment,
    cigar_element_t* const global_cigar_buffer,
    match_alignment_t* const match_alignment) {
  // Parameters
  vector_t* const cigar_vector = matches->cigar_vector;
  // Reserve CIGAR
  match_alignment->cigar_offset = vector_get_used(cigar_vector);
  const uint64_t max_local_cigar_length = local_alignment->local_end_cigar-local_alignment->local_begin_cigar+2;
  vector_reserve_additional(cigar_vector,max_local_cigar_length);
  cigar_element_t* local_cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t);
  uint64_t cigar_length = 0;
  // Begin trim
  const uint64_t begin_trim = pattern->clip_left + local_alignment->local_begin_key_pos;
  if (begin_trim > 0) {
    local_cigar_buffer->type = cigar_del;
    local_cigar_buffer->length = begin_trim;
    local_cigar_buffer->attributes = cigar_attr_trim;
    ++local_cigar_buffer;
    ++cigar_length;
  }
  // Compose local CIGAR
  uint64_t i;
  for (i=local_alignment->local_begin_cigar;i<local_alignment->local_end_cigar;++i) {
    *local_cigar_buffer = global_cigar_buffer[i];
    ++local_cigar_buffer;
    ++cigar_length;
  }
  // End trim
  const uint64_t end_trim = pattern->clip_right + (pattern->key_length-local_alignment->local_end_key_pos);
  if (end_trim > 0) {
    local_cigar_buffer->type = cigar_del;
    local_cigar_buffer->length = end_trim;
    local_cigar_buffer->attributes = cigar_attr_trim;
    ++local_cigar_buffer;
    ++cigar_length;
  }
  vector_add_used(cigar_vector,cigar_length);
  match_alignment->cigar_length = cigar_length;
}
void match_align_swg_local_alignment_add_match(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    const uint64_t match_position,
    match_align_swg_local_t* const local_alignment,
    cigar_element_t* const global_cigar_buffer) {
  // Parameters
  vector_t* const cigar_vector = matches->cigar_vector;
  match_trace_t match_trace;
  match_alignment_t* const match_alignment = &match_trace.match_alignment;
  // Configure match-trace
  match_trace.type = match_type_local;
  match_trace.text_trace = text_trace;
  match_trace.sequence_name = NULL;
  match_trace.text_position = UINT64_MAX;
  match_trace.match_scaffold = match_scaffold;
  // Position
  match_alignment->match_position = match_position + local_alignment->local_begin_text_pos;
  match_alignment->match_text_offset = match_alignment->match_position - text_trace->position;
  // Compose CIGAR
  match_align_swg_local_alignment_add_match_compose_cigar(
      matches,pattern,text_trace,local_alignment,
      global_cigar_buffer,match_alignment);
  // Compute distances
  const uint64_t cigar_offset = match_alignment->cigar_offset;
  const uint64_t cigar_length = match_alignment->cigar_length;
  match_trace.event_distance = matches_cigar_compute_event_distance(cigar_vector,cigar_offset,cigar_length);
  match_trace.edit_distance = matches_cigar_compute_edit_distance(cigar_vector,cigar_offset,cigar_length);
  match_trace.swg_score = align_swg_score_cigar_excluding_clipping(
      &search_parameters->swg_penalties,cigar_vector,cigar_offset,cigar_length);
  match_trace.error_quality =
      matches_cigar_compute_error_quality(
          cigar_vector,cigar_offset,cigar_length,
          pattern->quality_mask,pattern->key_length);
  match_alignment->effective_length = matches_cigar_compute_effective_length(cigar_vector,cigar_offset,cigar_length);
  match_trace.text = text_trace->text + match_alignment->match_text_offset;
  match_trace.text_length = match_alignment->effective_length;
  // Add to local pending
  matches_add_match_trace_local_pending(matches,&match_trace);
}
/*
 * SWG-Local Alignment
 */
void match_align_swg_local_alignment(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    match_trace_t* const match_trace) {
  // Parameters
  const swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const int32_t local_min_swg_threshold = search_parameters->alignment_local_min_swg_threshold_nominal;
  const uint64_t local_min_identity = search_parameters->alignment_local_min_identity_nominal;
  vector_t* const cigar_vector = matches->cigar_vector;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  const uint64_t match_position = match_alignment->match_position;
  // Prepare CIGARs
  const uint64_t cigar_length = match_alignment->cigar_length;
  cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,match_alignment->cigar_offset,cigar_element_t);
  // Local Max
  match_align_swg_local_t local_max = {
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
  uint64_t key_end = 0, text_end = 0;
  uint64_t i, local_identity = 0;
  int64_t local_score = 0, max_local_score = SWG_SCORE_MIN;
  cigar_t previous_type = cigar_null;
  for (i=0;i<cigar_length;) {
    // Add CIGAR operation
    switch (cigar_buffer[i].type) {
      case cigar_match: {
        const int32_t match_length = cigar_buffer[i].length;
        local_score += align_swg_score_match(swg_penalties,match_length);
        local_identity += match_length;
        text_end += match_length;
        key_end += match_length;
        break;
      }
      case cigar_mismatch:
        local_score += align_swg_score_mismatch(swg_penalties);
        ++text_end;
        ++key_end;
        break;
      case cigar_ins: {
        const int32_t indel_length = cigar_buffer[i].length;
        local_score += align_swg_score_insertion(swg_penalties,indel_length);
        text_end += indel_length;
        if (previous_type == cigar_del) {
          local_score = 0;
          local_identity = 0;
          local_max.local_score = 0;
          local_max.local_identity = 0;
          local_max.local_begin_key_pos = key_end;
          local_max.local_begin_text_pos = text_end;
          previous_type = cigar_buffer[i].type;
          local_max.local_begin_cigar = ++i;
          continue; // Reset
        }
        break;
      }
      case cigar_del: {
        const int32_t indel_length = cigar_buffer[i].length;
        local_score += align_swg_score_deletion(swg_penalties,indel_length);
        key_end += indel_length;
        if (previous_type == cigar_ins) {
          local_score = 0;
          local_identity = 0;
          local_max.local_score = 0;
          local_max.local_identity = 0;
          local_max.local_begin_key_pos = key_end;
          local_max.local_begin_text_pos = text_end;
          previous_type = cigar_buffer[i].type;
          local_max.local_begin_cigar = ++i;
          continue; // Reset
        }
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
    // Next
    previous_type = cigar_buffer[i].type;
    ++i;
    // Check local score
    if (local_score > local_max.local_score) {
      // Max local score
      local_max.local_score = local_score;
      local_max.local_identity = local_identity;
      local_max.local_end_key_pos = key_end;
      local_max.local_end_text_pos = text_end;
      local_max.local_end_cigar = i;
    } else if (local_score < 0) {
      // Check local swg-score and identity
      if (local_max.local_score >= local_min_swg_threshold &&
          local_max.local_identity >= local_min_identity) {
        match_align_swg_local_alignment_add_match(
            matches,search_parameters,pattern,text_trace,
            match_scaffold,match_position,&local_max,cigar_buffer);
        // Update positions & score
        max_local_score = MAX(max_local_score,local_max.local_score);
      }
      // Reset
      local_score = 0;
      local_identity = 0;
      local_max.local_score = 0;
      local_max.local_identity = 0;
      local_max.local_begin_key_pos = key_end;
      local_max.local_begin_text_pos = text_end;
      local_max.local_begin_cigar = i;
    }
  }
  // Check local swg-score and identity
  if (local_max.local_score >= local_min_swg_threshold &&
      local_max.local_identity >= local_min_identity) {
    match_align_swg_local_alignment_add_match(
        matches,search_parameters,pattern,text_trace,
        match_scaffold,match_position,&local_max,cigar_buffer);
    max_local_score = MAX(max_local_score,local_max.local_score);
  }
  // Finish
  match_trace->swg_score = max_local_score;
}
