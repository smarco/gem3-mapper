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
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack) {
  // Parameters
  uint8_t* const key = align_input->key;
  uint8_t* const text = align_input->text;
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
        // Force (re)computing the CIGAR from the alignment-region (missing CIGAR or never computed)
        match_align_input_t align_chunk_input = {
            .key = key + region_key_begin,
            .key_length = key_matching_length,
            .text = text + region_text_begin,
            .text_length = text_matching_length,
        };
        match_align_parameters_t align_chunk_parameters = {
            .max_bandwidth = region_error+1, // FIXME Include CIGAR when generating approximate regions
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
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    const bool begin_free,const bool end_free,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t max_aligned_gap_length = align_parameters->max_aligned_gap_length;
  vector_t* const cigar_vector = matches->cigar_vector;
  // Align gap between alignment-regions
  match_align_input_t align_chunk_input;
  align_chunk_input.key = align_input->key+key_chunk_begin_offset;
  align_chunk_input.key_length = key_chunk_length;
  align_chunk_input.text = align_input->text + text_chunk_begin_offset;
  align_chunk_input.text_length = text_chunk_length;
  // Keep CIGAR state
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t cigar_vector_used = vector_get_used(cigar_vector);
  // SWG-Align
  if (key_chunk_length > max_aligned_gap_length || text_chunk_length > max_aligned_gap_length) {
    match_alignment->score = SWG_SCORE_MIN;
  } else {
    align_swg(&align_chunk_input,align_parameters,begin_free,end_free,match_alignment,cigar_vector,mm_stack);
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
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    mm_stack_t* const mm_stack) {
  return match_align_swg_region(matches,align_input,align_parameters,match_alignment,
      key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,
      false,false,mm_stack);
}
bool match_align_swg_begin_region(
    matches_t* const matches,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    mm_stack_t* const mm_stack) {
  return match_align_swg_region(matches,align_input,align_parameters,match_alignment,
      key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,
      true,false,mm_stack);
}
bool match_align_swg_end_region(
    matches_t* const matches,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,
    const uint64_t key_chunk_begin_offset,
    const uint64_t key_chunk_length,
    const uint64_t text_chunk_begin_offset,
    const uint64_t text_chunk_length,
    mm_stack_t* const mm_stack) {
  return match_align_swg_region(matches,align_input,align_parameters,match_alignment,
      key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,
      false,true,mm_stack);
}
/*
 * Compute alignment type (local/global) wrt identity/score thresholds
 */
void match_align_swg_compute_alignment_type(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_parameters_t* const align_parameters) {
  // Parameters
  vector_t* const cigar_vector = matches->cigar_vector;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Compute matching bases (identity) + score
  const uint64_t cigar_offset = match_alignment->cigar_offset;
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t matching_bases = matches_cigar_compute_matching_bases(cigar_vector,cigar_offset,cigar_length);
  if (matching_bases < align_parameters->global_min_identity) {
    PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_MATCHING_BASES);
    match_trace->type = match_type_local; // Local (by identity)
  } else {
    match_trace->swg_score = align_swg_score_cigar(
        align_parameters->swg_penalties,cigar_vector,cigar_offset,cigar_length);
    if (match_trace->swg_score < align_parameters->global_min_swg_threshold) {
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
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t max_bandwidth = align_parameters->max_bandwidth;
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
  const bool feasible_alignment = match_align_swg_begin_region(matches,
      align_input,align_parameters,match_alignment,key_chunk_begin_offset,
      key_chunk_length,text_chunk_begin_offset,text_chunk_length,mm_stack);
  match_begin_position = (feasible_alignment) ?
      match_alignment->match_position : match_begin_position + text_chunk_length;
  return match_begin_position;
}
void match_align_swg_chain_scaffold_from_last(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t key_length = align_input->key_length;
  const uint64_t text_length = align_input->text_length;
  const uint64_t max_bandwidth = align_parameters->max_bandwidth;
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
  match_align_swg_end_region(matches,align_input,align_parameters,match_alignment,
      key_chunk_begin_offset,key_chunk_length,text_chunk_begin_offset,text_chunk_length,mm_stack);
}
void match_align_swg_chain_scaffold(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,
    mm_stack_t* const mm_stack) {
  // Parameters
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment_region_t* const alignment_regions = match_scaffold->alignment_regions;
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  vector_t* const cigar_vector = matches->cigar_vector;
  // Chain to first alignment-region
  const uint64_t match_begin_position = match_align_swg_chain_scaffold_to_first(
      matches,match_trace,align_input,align_parameters,match_scaffold,mm_stack);
  // Add first alignment-region
  match_align_swg_add_alignment_region(alignment_regions,
      align_input,align_parameters,match_alignment,cigar_vector,mm_stack);
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
        matches,align_input,align_parameters,match_alignment,
        key_chunk_begin_offset,key_chunk_length,
        text_chunk_begin_offset,text_chunk_length,mm_stack);
    // Add alignment-region
    match_align_swg_add_alignment_region(current_alignment_region,
        align_input,align_parameters,match_alignment,cigar_vector,mm_stack);
  }
  // Chain from last alignment-region
  match_align_swg_chain_scaffold_from_last(
      matches,match_trace,align_input,align_parameters,match_scaffold,mm_stack);
  // Post-processing
  match_alignment->score = 0; // Deferred calculation
  match_alignment->match_position = match_begin_position; // Restore match position
}
void match_align_swg(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,
    mm_stack_t* const mm_stack) {
  // Configure match-alignment
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment->match_position = align_input->text_position;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  // Check the number of alignment-regions
  if (align_parameters->alignment_force_full_swg) {
    // Add left trim
    if (align_input->key_trim_left > 0) {
      matches_cigar_vector_append_deletion(matches->cigar_vector,
          &match_alignment->cigar_length,align_input->key_trim_left,cigar_attr_trim);
    }
    // Full SWG
    match_align_input_t align_chunk_input;
    align_chunk_input.key = align_input->key + align_input->key_trim_left;
    align_chunk_input.key_length = align_input->key_length - align_input->key_trim_left - align_input->key_trim_right;
    align_chunk_input.text = align_input->text;
    align_chunk_input.text_length = align_input->text_length;
    align_swg(&align_chunk_input,align_parameters,true,true,
        &match_trace->match_alignment,matches->cigar_vector,mm_stack);
    // Add right trim
    if (align_input->key_trim_right > 0) {
      matches_cigar_vector_append_deletion(matches->cigar_vector,
          &match_alignment->cigar_length,align_input->key_trim_right,cigar_attr_trim);
    }
  } else {
    // Check null scaffolding
    if (match_scaffold->num_alignment_regions==0) {
      match_alignment->score = SWG_SCORE_MIN;
      return;
    }
    // Chain alignment-regions and align gaps (SWG)
    match_align_swg_chain_scaffold(matches,match_trace,
        align_input,align_parameters,match_scaffold,mm_stack);
  }
}
