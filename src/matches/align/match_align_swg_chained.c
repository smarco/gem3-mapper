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

#include "matches/align/match_align_swg_chained.h"
#include "align/align_swg.h"
#include "align/align_swg_banded.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * SWG Chained Add Region
 */
void match_align_swg_chained_add_region(
    const match_alignment_region_t* const match_alignment_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const uint64_t region_key_begin = match_alignment_region_get_key_begin(match_alignment_region);
  const uint64_t region_key_end = match_alignment_region_get_key_end(match_alignment_region);
  const uint64_t region_text_begin = match_alignment_region_get_text_begin(match_alignment_region);
  const uint64_t region_text_end = match_alignment_region_get_text_end(match_alignment_region);
  const uint64_t region_error = match_alignment_region_get_error(match_alignment_region);
  const uint64_t region_cigar_buffer_offset = match_alignment_region_get_cigar_buffer_offset(match_alignment_region);
  const uint64_t region_cigar_length = match_alignment_region_get_cigar_length(match_alignment_region);
  const uint64_t key_matching_length = region_key_end - region_key_begin;
  const uint64_t text_matching_length = region_text_end - region_text_begin;
  // Select matching-region type
  switch (match_alignment_region_get_type(match_alignment_region)) {
    case match_alignment_region_exact:
      matches_cigar_vector_append_match(cigar_vector,&match_alignment->cigar_length,key_matching_length,cigar_attr_none);
      match_alignment->score = align_swg_score_match(swg_penalties,key_matching_length);
      break;
    case match_alignment_region_approximate:
      if (region_cigar_length > 0) {
        // Copy the CIGAR from the alignment-region
        uint64_t i;
        match_alignment->score = 0;
        for (i=0;i<region_cigar_length;++i) {
          cigar_element_t* const cigar_element = vector_get_elm(cigar_vector,region_cigar_buffer_offset+i,cigar_element_t);
          match_alignment->score += align_swg_score_cigar_element(swg_penalties,cigar_element);
          matches_cigar_vector_append_cigar_element(cigar_vector,&match_alignment->cigar_length,cigar_element);
        }
      } else {
        // Parameters
        const uint8_t* const key = pattern->key + region_key_begin;
        const uint64_t key_length = key_matching_length;
        uint8_t* const text = text_trace->text + region_text_begin;
        const uint64_t text_length = text_matching_length;
        const uint64_t max_bandwidth = region_error+1; // TODO Include CIGAR when generating approximate regions
        const bool left_gap_alignment = (text_trace->strand == Forward);
        // Force (re)computing the CIGAR from the alignment-region (missing CIGAR or never computed)
        align_swg(
            key,key_length,text,text_length,
            swg_penalties,max_bandwidth,
            false,false,left_gap_alignment,
            match_alignment,cigar_vector,mm_allocator);
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * SWG Align Gaps/Endings
 */
void match_align_swg_chained_gap_full(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    const bool local_alignment,
    match_alignment_t* const match_alignment,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  vector_t* const cigar_vector = matches->cigar_vector;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
  // Parameters Gap
  const uint8_t* const key = pattern->key + key_chunk_begin_offset;
  const uint64_t key_length = key_chunk_length;
  uint8_t* const text = text_trace->text + text_chunk_begin_offset;
  const uint64_t text_length = text_chunk_length;
  const bool left_gap_alignment = (text_trace->strand == Forward);
  // Keep CIGAR state
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t cigar_vector_used = vector_get_used(cigar_vector);
  // Align gap between alignment-regions
  align_swg(
      key,key_length,text,text_length,
      swg_penalties,max_bandwidth,
      false,false,left_gap_alignment,
      match_alignment,cigar_vector,mm_allocator);
  // Check alignment result
  if (match_alignment->score == SWG_SCORE_MIN) {
    // Restore the CIGAR
    match_alignment->cigar_length = cigar_length;
    vector_set_used(cigar_vector,cigar_vector_used);
    // Bride the gap
    if (key_chunk_length > 0) { // Delete the read chunk
      matches_cigar_vector_append_deletion(cigar_vector,
          &match_alignment->cigar_length,key_chunk_length,cigar_attr_none);
    }
    if (text_chunk_length > 0) { // Insert the text chunk
      matches_cigar_vector_append_insertion(cigar_vector,
          &match_alignment->cigar_length,text_chunk_length,cigar_attr_none);
    }
  }
}
bool match_align_swg_chained_gap_extending(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    const bool local_alignment,
    const int32_t left_base_score,
    const int32_t right_base_score,
    match_alignment_t* const match_alignment,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  vector_t* const cigar_vector = matches->cigar_vector;
  const bool left_gap_alignment = (text_trace->strand == Forward);
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
  // Parameters Extend Forward
  const uint8_t* const key_forward = pattern->key + key_chunk_begin_offset;
  const uint64_t key_length_forward = key_chunk_length;
  uint8_t* const text_forward = text_trace->text + text_chunk_begin_offset;
  const uint64_t text_length_forward = text_chunk_length;
  // SWG-Align Extend Forward
  uint64_t total_key_aligned_forward, total_text_aligned_forward;
  align_swg_banded_extend(
      key_forward,key_length_forward,
      text_forward,text_length_forward,
      false,local_alignment,left_base_score,
      swg_penalties,max_bandwidth,
      left_gap_alignment,match_alignment,
      &total_key_aligned_forward,
      &total_text_aligned_forward,
      cigar_vector,mm_allocator);
  // Check bases aligned
  if (total_key_aligned_forward == key_length_forward &&
      total_text_aligned_forward == text_length_forward) {
    return true; // Gap fully bridged
  }
  // Parameters Extend Reverse
  const uint8_t* const key_reverse = key_forward + total_key_aligned_forward;
  const uint64_t key_length_reverse = key_length_forward - total_key_aligned_forward;
  uint8_t* const text_reverse = text_forward + total_text_aligned_forward;
  const uint64_t text_length_reverse = text_length_forward - total_text_aligned_forward;
  const uint64_t cigar_gap_position = match_alignment->cigar_offset + match_alignment->cigar_length;
  // SWG-Align Extend Reverse
  int32_t swg_score_forward = match_alignment->score;
  int32_t swg_score_reverse = 0;
  uint64_t total_key_aligned_reverse = 0, total_text_aligned_reverse = 0;
  if (key_length_reverse > 0 && text_length_reverse > 0) {
    align_swg_banded_extend(
        key_reverse,key_length_reverse,
        text_reverse,text_length_reverse,
        true,local_alignment,right_base_score,
        swg_penalties,max_bandwidth,
        left_gap_alignment,match_alignment,
        &total_key_aligned_reverse,
        &total_text_aligned_reverse,
        cigar_vector,mm_allocator);
    swg_score_reverse = match_alignment->score;
  }
  // Add key gap
  const uint64_t key_gap = key_chunk_length - total_key_aligned_forward - total_key_aligned_reverse;
  int32_t swg_score_key_gap = 0;
  if (key_gap > 0) {
    cigar_element_t cigar_element;
    cigar_element.type = cigar_del;
    cigar_element.length = key_gap;
    cigar_element.attributes = cigar_attr_none;
    matches_cigar_vector_insert_cigar_element(cigar_vector,cigar_gap_position,&cigar_element);
    ++(match_alignment->cigar_length);
    swg_score_key_gap = align_swg_score_deletion(swg_penalties,key_gap);
  }
  // Add text gap
  const uint64_t text_gap = text_chunk_length - total_text_aligned_forward - total_text_aligned_reverse;
  int32_t swg_score_text_gap = 0;
  if (text_gap > 0) {
    cigar_element_t cigar_element;
    cigar_element.type = cigar_ins;
    cigar_element.length = text_gap;
    cigar_element.attributes = cigar_attr_none;
    matches_cigar_vector_insert_cigar_element(cigar_vector,cigar_gap_position,&cigar_element);
    ++(match_alignment->cigar_length);
    swg_score_text_gap = align_swg_score_insertion(swg_penalties,text_gap);
  }
  // Return
  const int32_t total_score = swg_score_forward + swg_score_reverse + swg_score_key_gap + swg_score_text_gap;
  if ((total_score <= 0) || (key_gap > 0 && text_gap > 0)) {
    match_alignment->score = swg_score_reverse;
    return false; // Gap not bridged
  } else {
    match_alignment->score = total_score;
    return true; // Gap fully bridged
  }
}
uint64_t match_align_swg_end_extending(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    const bool reverse_extension,
    const bool local_alignment,
    const int32_t extension_base_score,
    match_alignment_t* const match_alignment,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  vector_t* const cigar_vector = matches->cigar_vector;
  // Parameters
  const uint8_t* const key = pattern->key + key_chunk_begin_offset;
  const uint64_t key_length = key_chunk_length;
  uint8_t* const text = text_trace->text + text_chunk_begin_offset;
  const uint64_t text_length = text_chunk_length;
  const bool left_gap_alignment = (text_trace->strand == Forward);
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
  // SWG-Align Extend
  uint64_t total_key_aligned, total_text_aligned;
  align_swg_banded_extend(
      key,key_length,text,text_length,reverse_extension,
      local_alignment,extension_base_score,swg_penalties,
      max_bandwidth,left_gap_alignment,match_alignment,
      &total_key_aligned,&total_text_aligned,
      cigar_vector,mm_allocator);
  // Check key-aligned (key-trim)
  if (total_key_aligned < key_length) {
    // Add final key-trim (Local maximum)
    const uint64_t trim_pos = reverse_extension ?
        match_alignment->cigar_offset :
        match_alignment->cigar_offset + match_alignment->cigar_length;
    cigar_element_t cigar_element;
    cigar_element.type = cigar_del;
    cigar_element.length = key_length-total_key_aligned;
    cigar_element.attributes = cigar_attr_none;
    matches_cigar_vector_insert_cigar_element(cigar_vector,trim_pos,&cigar_element);
    ++(match_alignment->cigar_length);
  }
  // Return match offset (text offset)
  return text_length - total_text_aligned;
}
/*
 * SWG Chained
 */
uint64_t match_align_swg_chained_to_first(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    const bool local_alignment,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment_region_t* const alignment_regions = match_scaffold->alignment_regions;
  const match_alignment_region_t* const first_alignment_region = alignment_regions;
  const uint64_t first_key_begin = match_alignment_region_get_key_begin(first_alignment_region);
  const uint64_t first_key_end = match_alignment_region_get_key_end(first_alignment_region);
  const uint64_t first_text_begin = match_alignment_region_get_text_begin(first_alignment_region);
  // Offsets
  const uint64_t key_chunk_length = first_key_begin;
  const uint64_t key_chunk_begin_offset = 0;
  const uint64_t text_chunk_length = BOUNDED_ADDITION(key_chunk_length,max_bandwidth,first_text_begin);
  const uint64_t text_chunk_begin_offset = first_text_begin - text_chunk_length;
  // Chain to first alignment-region
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const int32_t generic_match_score = swg_penalties->generic_match_score;
  const int32_t extension_base_score = (first_key_end-first_key_begin)*generic_match_score;
  uint64_t match_begin_position = match_alignment->match_position; // Save match position
  const uint64_t match_position_offset =
      match_align_swg_end_extending(
          matches,search_parameters,pattern,text_trace,
          key_chunk_begin_offset,key_chunk_length,
          text_chunk_begin_offset,text_chunk_length,
          true,local_alignment,extension_base_score,
          match_alignment,mm_allocator);
  // Set position
  return match_begin_position + text_chunk_begin_offset + match_position_offset;
}
void match_align_swg_chained_from_last(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    const bool local_alignment,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  const uint64_t text_length = text_trace->text_length;
  const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment_region_t* const alignment_regions = match_scaffold->alignment_regions;
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  const match_alignment_region_t* const last_alignment_region = alignment_regions + (num_alignment_regions-1);
  const uint64_t last_key_begin = match_alignment_region_get_key_begin(last_alignment_region);
  const uint64_t last_key_end = match_alignment_region_get_key_end(last_alignment_region);
  const uint64_t last_text_end = match_alignment_region_get_text_end(last_alignment_region);
  // Offsets
  const uint64_t key_chunk_begin_offset = last_key_end;
  const uint64_t key_chunk_length = key_length - last_key_end;
  const uint64_t text_chunk_end_offset = BOUNDED_ADDITION(last_text_end,key_chunk_length+max_bandwidth,text_length);
  const uint64_t text_chunk_begin_offset = last_text_end;
  const uint64_t text_chunk_length = text_chunk_end_offset-last_text_end;
  // Chain from last alignment-region
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const int32_t generic_match_score = swg_penalties->generic_match_score;
  const int32_t extension_base_score = (last_key_end-last_key_begin)*generic_match_score;
  match_align_swg_end_extending(
      matches,search_parameters,pattern,text_trace,
      key_chunk_begin_offset,key_chunk_length,
      text_chunk_begin_offset,text_chunk_length,
      false,local_alignment,extension_base_score,
      match_alignment,mm_allocator);
}
bool match_align_swg_chained_gap(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    const uint64_t match_alignment_region_idx,
    const bool local_alignment,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment_region_t* const alignment_regions = match_scaffold->alignment_regions;
  const match_alignment_region_t* const prev_alignment_region = alignment_regions + (match_alignment_region_idx-1);
  const match_alignment_region_t* const current_alignment_region = alignment_regions + match_alignment_region_idx;
  const uint64_t prev_key_begin = match_alignment_region_get_key_begin(prev_alignment_region);
  const uint64_t prev_key_end = match_alignment_region_get_key_end(prev_alignment_region);
  const uint64_t prev_text_end = match_alignment_region_get_text_end(prev_alignment_region);
  const uint64_t current_key_begin = match_alignment_region_get_key_begin(current_alignment_region);
  const uint64_t current_key_end = match_alignment_region_get_key_end(current_alignment_region);
  const uint64_t current_text_begin = match_alignment_region_get_text_begin(current_alignment_region);
  // Offsets
  const uint64_t key_chunk_begin_offset = prev_key_end;
  const uint64_t key_chunk_length = current_key_begin - key_chunk_begin_offset;
  const uint64_t text_chunk_begin_offset = prev_text_end;
  const uint64_t text_chunk_length = current_text_begin - text_chunk_begin_offset;
  // Align the gap
  if (local_alignment) {
    swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
    const int32_t generic_match_score = swg_penalties->generic_match_score;
    const int32_t right_base_score = (current_key_end-current_key_begin)*generic_match_score;
    const int32_t left_base_score = (prev_key_end-prev_key_begin)*generic_match_score;
    const bool gap_filled = match_align_swg_chained_gap_extending(
          matches,search_parameters,pattern,text_trace,
          key_chunk_begin_offset,key_chunk_length,
          text_chunk_begin_offset,text_chunk_length,
          local_alignment,left_base_score,right_base_score,
          match_alignment,mm_allocator);
    return gap_filled;
  } else {
    match_align_swg_chained_gap_full(
          matches,search_parameters,pattern,text_trace,
          key_chunk_begin_offset,key_chunk_length,
          text_chunk_begin_offset,text_chunk_length,
          local_alignment,match_alignment,mm_allocator);
    return true;
  }
}
void match_align_swg_chained(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    const bool local_alignment,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment_region_t* const alignment_regions = match_scaffold->alignment_regions;
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  vector_t* const cigar_vector = matches->cigar_vector;
  // Chain to first alignment-region
  const uint64_t match_position =
      match_align_swg_chained_to_first(
          matches,search_parameters,pattern,text_trace,
          match_scaffold,local_alignment,match_trace,mm_allocator);
  // Add first alignment-region
  match_align_swg_chained_add_region(
      alignment_regions,search_parameters,pattern,
      text_trace,match_alignment,cigar_vector,mm_allocator);
  // Chain alignment-regions
  uint64_t i;
  for (i=1;i<num_alignment_regions;++i) {
    // Align the gap
    match_align_swg_chained_gap(
        matches,search_parameters,pattern,text_trace,
        match_scaffold,i,local_alignment,match_trace,mm_allocator);
    if (match_alignment->score == SWG_SCORE_MIN) return;
    // Add alignment-region
    match_align_swg_chained_add_region(
        alignment_regions+i,search_parameters,pattern,
        text_trace,match_alignment,cigar_vector,mm_allocator);
  }
  // Chain from last alignment-region
  match_align_swg_chained_from_last(
      matches,search_parameters,pattern,text_trace,
      match_scaffold,local_alignment,match_trace,mm_allocator);
  // Restore match position
  match_alignment->match_position = match_position;
  match_alignment->match_text_offset = match_position - text_trace->position;
}
