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

#include "matches/align/match_align_swg.h"
#include "matches/align/match_align_swg_chained.h"
#include "matches/align/match_align.h"
#include "align/alignment.h"
#include "align/align_swg.h"
#include "align/align_swg_banded.h"
#include "archive/archive_text_rl.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Constants
 */
//#define MATCH_ALIGN_SWG_CHECK

/*
 * SWG Align Check
 */
void match_align_swg_check(
    matches_t* const matches,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_alignment_t* const match_alignment,
    mm_allocator_t* const mm_allocator) {
  match_alignment->effective_length = matches_cigar_compute_effective_length(
      matches->cigar_vector,match_alignment->cigar_offset,match_alignment->cigar_length);
  const bool check_correct = match_alignment_check(
      stderr,match_alignment,
      pattern->key,pattern->key_length,
      text_trace->text,text_trace->text_length,
      matches->cigar_vector,true,mm_allocator);
  if (!check_correct) {
    match_alignment_print_pretty(
      stderr,match_alignment,
      pattern->key,pattern->key_length,
      text_trace->text+match_alignment->match_text_offset,
      match_alignment->effective_length,
      matches->cigar_vector,mm_allocator);
  }
}
/*
 * Compute alignment type (local/global) wrt identity/score thresholds
 */
void match_align_swg_compute_alignment_type(
    matches_t* const matches,
    match_trace_t* const match_trace,
    pattern_t* const pattern,
    search_parameters_t* const search_parameters) {
  // Parameters
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  const uint64_t global_min_identity = search_parameters->alignment_global_min_identity_nominal;
  const int32_t global_min_swg_threshold = search_parameters->alignment_global_min_swg_threshold_nominal;
  vector_t* const cigar_vector = matches->cigar_vector;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Compute distance + edit distance + effective-length
  const uint64_t cigar_offset = match_alignment->cigar_offset;
  const uint64_t cigar_length = match_alignment->cigar_length;
  const uint64_t matching_bases = matches_cigar_compute_matching_bases(cigar_vector,cigar_offset,cigar_length);
  match_alignment->effective_length = matches_cigar_compute_effective_length(cigar_vector,cigar_offset,cigar_length);
  match_trace->event_distance = matches_cigar_compute_event_distance(cigar_vector,cigar_offset,cigar_length);
  match_trace->edit_distance = matches_cigar_compute_edit_distance(cigar_vector,cigar_offset,cigar_length);
  match_trace->text_length = match_alignment->effective_length;
  match_trace->swg_score = align_swg_score_cigar_excluding_clipping(
      swg_penalties,cigar_vector,cigar_offset,cigar_length);
  match_trace->error_quality =
      matches_cigar_compute_error_quality(
          cigar_vector,cigar_offset,cigar_length,
          pattern->quality_mask,pattern->key_length);
  // Classify
  if (matching_bases < global_min_identity) {
    PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_MATCHING_BASES);
    match_trace->type = match_type_local; // Local (by identity)
  } else {
    if (match_trace->swg_score < global_min_swg_threshold) {
      PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_SWG_THRESHOLD);
      match_trace->type = match_type_local; // Local (by score)
    } else {
      match_trace->type = match_type_regular; // Regular Alignment (global)
    }
  }
}
/*
 * SWG Align
 */
void match_align_swg(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    match_scaffold_t* const match_scaffold,
    const bool local_alignment,
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
        key,key_length,text,text_length,
        swg_penalties,max_bandwidth,
        true,true,left_gap_alignment,
        match_alignment,cigar_vector,mm_allocator);
    // Add right trim
    if (text_trace->text_padded_right > 0) {
      matches_cigar_vector_append_deletion(
          cigar_vector,&match_alignment->cigar_length,
          text_trace->text_padded_right,cigar_attr_trim);
    }
    // Adjust offset
    match_alignment->match_text_offset = match_alignment->match_position - text_trace->position;
  } else {
    // Check null scaffolding
    if (match_scaffold->num_alignment_regions==0) {
      match_alignment->score = SWG_SCORE_MIN;
      return;
    }
    // Chain alignment-regions and align gaps (SWG)
    match_align_swg_chained(
        matches,search_parameters,pattern,text_trace,
        match_scaffold,local_alignment,match_trace,mm_allocator);
  }
  // Check
#ifdef MATCH_ALIGN_SWG_CHECK
  if (match_alignment->score != SWG_SCORE_MIN) {
    match_align_swg_check(matches,pattern,
        text_trace,match_alignment,mm_allocator);
  }
#endif
}
