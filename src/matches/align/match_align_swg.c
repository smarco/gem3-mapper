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

#include "align/alignment.h"
#include "matches/align/match_align_swg.h"
#include "matches/align/match_align.h"
#include "align/align_swg.h"
#include "archive/archive_text_rl.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * SWG Align alignment-region (@match_alignment_region)
 */
void match_align_swg_add_alignment_region(
    const match_alignment_region_t* const match_alignment_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator) {
  // Parameters
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
      break;
    case match_alignment_region_approximate:
      if (region_cigar_length > 0) {
        // Copy the CIGAR from the alignment-region
        uint64_t i;
        for (i=0;i<region_cigar_length;++i) {
          cigar_element_t* const scaffolding_elm =
              vector_get_elm(cigar_vector,region_cigar_buffer_offset+i,cigar_element_t);
          matches_cigar_vector_append_cigar_element(cigar_vector,&match_alignment->cigar_length,scaffolding_elm);
        }
      } else {
        // Parameters
        const uint8_t* const key = pattern->key + region_key_begin;
        const uint64_t key_length = key_matching_length;
        uint8_t* const text = text_trace->text + region_text_begin;
        const uint64_t text_length = text_matching_length;
        const uint64_t max_bandwidth = region_error+1; // TODO Include CIGAR when generating approximate regions
        const bool left_gap_alignment = (text_trace->strand == Forward);
        swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
        // Force (re)computing the CIGAR from the alignment-region (missing CIGAR or never computed)
        align_swg(
            key,key_length,text,text_length,swg_penalties,
            false,false,max_bandwidth,left_gap_alignment,
            match_alignment,cigar_vector,mm_allocator);
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
void match_align_swg_add_gap(
    matches_t* const matches,
    match_alignment_t* const match_alignment,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_length,
    const bool trim) {
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
bool match_align_swg_region(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    const bool begin_free,
    const bool end_free,
    match_alignment_t* const match_alignment,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint64_t max_aligned_gap_length = search_parameters->alignment_max_aligned_gap_length_nominal;
  vector_t* const cigar_vector = matches->cigar_vector;
  // Keep CIGAR state
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t cigar_vector_used = vector_get_used(cigar_vector);
  // SWG-Align
  if (key_chunk_length > max_aligned_gap_length || text_chunk_length > max_aligned_gap_length) {
    match_alignment->score = SWG_SCORE_MIN;
  } else {
    // Parameters
    const uint8_t* const key = pattern->key + key_chunk_begin_offset;
    const uint64_t key_length = key_chunk_length;
    uint8_t* const text = text_trace->text + text_chunk_begin_offset;
    const uint64_t text_length = text_chunk_length;
    const bool left_gap_alignment = (text_trace->strand == Forward);
    swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
    const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
    // Align gap between alignment-regions
    align_swg(
        key,key_length,text,text_length,swg_penalties,
        begin_free,end_free,max_bandwidth,left_gap_alignment,
        match_alignment,cigar_vector,mm_allocator);
  }
  // Check alignment result
  if (match_alignment->score == SWG_SCORE_MIN) {
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
bool match_align_swg_middle_region(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    match_alignment_t* const match_alignment,
    mm_allocator_t* const mm_allocator) {
  return match_align_swg_region(
      matches,search_parameters,pattern,text_trace,
      key_chunk_begin_offset,key_chunk_length,
      text_chunk_begin_offset,text_chunk_length,
      false,false,match_alignment,mm_allocator);
}
bool match_align_swg_begin_region(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    match_alignment_t* const match_alignment,
    mm_allocator_t* const mm_allocator) {
  return match_align_swg_region(
      matches,search_parameters,pattern,text_trace,
      key_chunk_begin_offset,key_chunk_length,
      text_chunk_begin_offset,text_chunk_length,
      true,false,match_alignment,mm_allocator);
}
bool match_align_swg_end_region(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    match_alignment_t* const match_alignment,
    mm_allocator_t* const mm_allocator) {
  return match_align_swg_region(
      matches,search_parameters,pattern,text_trace,
      key_chunk_begin_offset,key_chunk_length,
      text_chunk_begin_offset,text_chunk_length,
      false,true,match_alignment,mm_allocator);
}
/*
 * Compute alignment type (local/global) wrt identity/score thresholds
 */
void match_align_swg_compute_alignment_type(
    matches_t* const matches,
    match_trace_t* const match_trace,
    search_parameters_t* const search_parameters) {
  // Parameters
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const uint64_t global_min_identity = search_parameters->alignment_global_min_identity_nominal;
  const uint64_t global_min_swg_threshold = search_parameters->alignment_global_min_swg_threshold_nominal;
  vector_t* const cigar_vector = matches->cigar_vector;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Compute matching bases (identity) + score
  const uint64_t cigar_offset = match_alignment->cigar_offset;
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t matching_bases = matches_cigar_compute_matching_bases(cigar_vector,cigar_offset,cigar_length);
  if (matching_bases < global_min_identity) {
    PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_MATCHING_BASES);
    match_trace->type = match_type_local; // Local (by identity)
  } else {
    match_trace->swg_score = align_swg_score_cigar(swg_penalties,cigar_vector,cigar_offset,cigar_length);
    if (match_trace->swg_score < global_min_swg_threshold) {
      PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_SWG_THRESHOLD);
      match_trace->type = match_type_local; // Local (by score)
    } else {
      match_trace->type = match_type_regular; // Regular Alignment (global)
    }
  }
}
/*
 * SWG Alignment
 */
uint64_t match_align_swg_chain_scaffold_to_first(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment_region_t* const alignment_regions = match_scaffold->alignment_regions;
  const match_alignment_region_t* const first_alignment_region = alignment_regions;
  const uint64_t first_key_begin = match_alignment_region_get_key_begin(first_alignment_region);
  const uint64_t first_text_begin = match_alignment_region_get_text_begin(first_alignment_region);
  // Offsets
  const uint64_t key_chunk_length = first_key_begin;
  const uint64_t key_chunk_begin_offset = 0;
  const uint64_t text_chunk_length = BOUNDED_ADDITION(key_chunk_length,max_bandwidth,first_text_begin);
  const uint64_t text_chunk_begin_offset = first_text_begin - text_chunk_length;
  // Chain to first alignment-region
  match_alignment->match_position += text_chunk_begin_offset; // Offset match position
  uint64_t match_begin_position = match_alignment->match_position; // Save match position
  const bool feasible_alignment =
      match_align_swg_begin_region(
          matches,search_parameters,pattern,text_trace,
          key_chunk_begin_offset,key_chunk_length,
          text_chunk_begin_offset,text_chunk_length,
          match_alignment,mm_allocator);
  match_begin_position = (feasible_alignment) ?
      match_alignment->match_position :
      match_begin_position + text_chunk_length;
  return match_begin_position;
}
void match_align_swg_chain_scaffold_from_last(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
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
  const uint64_t last_key_end = match_alignment_region_get_key_end(last_alignment_region);
  const uint64_t last_text_end = match_alignment_region_get_text_end(last_alignment_region);
  // Offsets
  const uint64_t key_chunk_begin_offset = last_key_end;
  const uint64_t key_chunk_length = key_length - last_key_end;
  const uint64_t text_chunk_end_offset = BOUNDED_ADDITION(last_text_end,key_chunk_length+max_bandwidth,text_length);
  const uint64_t text_chunk_begin_offset = last_text_end;
  const uint64_t text_chunk_length = text_chunk_end_offset-last_text_end;
  // Chain from last alignment-region
  match_align_swg_end_region(
      matches,search_parameters,pattern,text_trace,
      key_chunk_begin_offset,key_chunk_length,
      text_chunk_begin_offset,text_chunk_length,
      match_alignment,mm_allocator);
}
void match_align_swg_chain_scaffold(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment_region_t* const alignment_regions = match_scaffold->alignment_regions;
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  vector_t* const cigar_vector = matches->cigar_vector;
  // Chain to first alignment-region
  const uint64_t match_begin_position =
      match_align_swg_chain_scaffold_to_first(
          matches,search_parameters,pattern,text_trace,
          match_scaffold,match_trace,mm_allocator);
  // Add first alignment-region
  match_align_swg_add_alignment_region(
      alignment_regions,search_parameters,pattern,text_trace,
      match_alignment,cigar_vector,mm_allocator);
  // Chain alignment-regions
  uint64_t i;
  for (i=1;i<num_alignment_regions;++i) {
    // Parameters
    const match_alignment_region_t* const prev_alignment_region = alignment_regions + (i-1);
    const match_alignment_region_t* const current_alignment_region = alignment_regions + i;
    const uint64_t prev_key_end = match_alignment_region_get_key_end(prev_alignment_region);
    const uint64_t prev_text_end = match_alignment_region_get_text_end(prev_alignment_region);
    const uint64_t current_key_begin = match_alignment_region_get_key_begin(current_alignment_region);
    const uint64_t current_text_begin = match_alignment_region_get_text_begin(current_alignment_region);
    // Offsets
    const uint64_t key_chunk_begin_offset = prev_key_end;
    const uint64_t key_chunk_length = current_key_begin - key_chunk_begin_offset;
    const uint64_t text_chunk_begin_offset = prev_text_end;
    const uint64_t text_chunk_length = current_text_begin - text_chunk_begin_offset;
    // Align the gap
    match_align_swg_middle_region(
        matches,search_parameters,pattern,text_trace,
        key_chunk_begin_offset,key_chunk_length,
        text_chunk_begin_offset,text_chunk_length,
        match_alignment,mm_allocator);
    // Add alignment-region
    match_align_swg_add_alignment_region(
        current_alignment_region,search_parameters,pattern,text_trace,
        match_alignment,cigar_vector,mm_allocator);
  }
  // Chain from last alignment-region
  match_align_swg_chain_scaffold_from_last(
      matches,search_parameters,pattern,text_trace,
      match_scaffold,match_trace,mm_allocator);
  // Post-processing
  match_alignment->score = 0; // Deferred calculation
  match_alignment->match_position = match_begin_position; // Restore match position
}
void match_align_swg(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  // Configure match-alignment
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment->match_position = text_trace->position;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  // Check the number of alignment-regions
  if (search_parameters->alignment_force_full_swg) {
    vector_t* const cigar_vector = matches->cigar_vector;
    // Add left trim
    if (text_trace->text_padded_left > 0) {
      matches_cigar_vector_append_deletion(
          cigar_vector,&match_alignment->cigar_length,
          text_trace->text_padded_left,cigar_attr_trim);
    }
    // Parameters
    const uint8_t* const key = pattern->key + text_trace->text_padded_left;
    const uint64_t key_length = pattern->key_length - text_trace->text_padded_left - text_trace->text_padded_right;
    uint8_t* const text = text_trace->text;
    const uint64_t text_length = text_trace->text_length;
    const bool left_gap_alignment = (text_trace->strand == Forward);
    swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
    const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
    // Full SWG
    align_swg(
        key,key_length,text,text_length,swg_penalties,
        true,true,max_bandwidth,left_gap_alignment,
        match_alignment,cigar_vector,mm_allocator);
    // Add right trim
    if (text_trace->text_padded_right > 0) {
      matches_cigar_vector_append_deletion(
          cigar_vector,&match_alignment->cigar_length,
          text_trace->text_padded_right,cigar_attr_trim);
    }
  } else {
    // Check null scaffolding
    if (match_scaffold->num_alignment_regions==0) {
      match_alignment->score = SWG_SCORE_MIN;
      return;
    }
    // Chain alignment-regions and align gaps (SWG)
    match_align_swg_chain_scaffold(
        matches,search_parameters,pattern,text_trace,
        match_scaffold,match_trace,mm_allocator);
  }
}
