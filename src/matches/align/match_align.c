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
#include "matches/align/match_align.h"
#include "matches/align/match_align_swg.h"
#include "matches/align/match_align_swg_local.h"
#include "matches/align/match_align_normalize.h"
#include "matches/matches_cigar.h"
#include "align/align_swg.h"
#include "align/align_bpm.h"
#include "archive/archive_text_rl.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Match Clipping
 */
void match_aling_add_clipping(
    match_trace_t* const match_trace,
    vector_t* const cigar_vector,
    const uint64_t sequence_clip_left,
    const uint64_t sequence_clip_right) {
  // Left trim
  if (sequence_clip_left > 0) {
    // Parameters
    match_alignment_t* const match_alignment = &match_trace->match_alignment;
    const uint64_t cigar_offset = match_alignment->cigar_offset;
    cigar_element_t* cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
    if (cigar_buffer[0].type==cigar_del) {
      cigar_buffer[0].length += sequence_clip_left;
      cigar_buffer[0].attributes = cigar_attr_trim;
    } else {
      // Reserve additional
      vector_reserve_additional(cigar_vector,1);
      vector_inc_used(cigar_vector); // Increment used
      cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
      // Shift CIGAR one position right
      uint64_t i;
      for (i=match_alignment->cigar_length;i>0;--i) {
        cigar_buffer[i] = cigar_buffer[i-1];
      }
      // Add clip
      cigar_buffer[0].type = cigar_del;
      cigar_buffer[0].attributes = cigar_attr_trim;
      cigar_buffer[0].length = sequence_clip_left;
      ++(match_alignment->cigar_length);
    }
  }
  // Right trim
  if (sequence_clip_right > 0) {
    // Parameters
    match_alignment_t* const match_alignment = &match_trace->match_alignment;
    const uint64_t cigar_offset = match_alignment->cigar_offset;
    cigar_element_t* cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
    const uint64_t cigar_length = match_alignment->cigar_length;
    if (cigar_buffer[cigar_length-1].type==cigar_del) {
      cigar_buffer[cigar_length-1].length += sequence_clip_right;
      cigar_buffer[cigar_length-1].attributes = cigar_attr_trim;
    } else {
      // Reserve additional
      vector_reserve_additional(cigar_vector,1);
      vector_inc_used(cigar_vector); // Increment used
      cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
      // Add clip
      cigar_buffer[cigar_length].type = cigar_del;
      cigar_buffer[cigar_length].attributes = cigar_attr_trim;
      cigar_buffer[cigar_length].length = sequence_clip_right;
      ++(match_alignment->cigar_length);
    }
  }
}
/*
 * Exact Alignment
 */
void match_align_exact(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    match_trace_t* const match_trace) {
  PROFILE_START(GP_MATCHES_ALIGN_EXACT,PROFILE_LEVEL);
  // Parameters
  const uint64_t key_length = pattern->key_length;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Configure match-trace
  match_trace->type = match_type_regular;
  match_trace->text_trace = text_trace;
  match_trace->text = NULL;
  match_trace->text_length = pattern->key_length;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->event_distance = 0;
  match_trace->edit_distance = 0;
  match_trace->swg_score = align_swg_score_match(&search_parameters->swg_penalties,(int32_t)key_length);
  match_trace->error_quality = (float)SEQUENCE_QUALITIES_MAX;
  // Set position/distance
  match_alignment->match_text_offset = alignment->alignment_tiles->text_begin_offset;
  match_alignment->match_position = text_trace->position + match_alignment->match_text_offset; // Adjust position
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  match_alignment->score = 0;
  match_alignment->effective_length = key_length;
  // Insert exact-match CIGAR
  matches_cigar_vector_append_match(
      matches->cigar_vector,&match_alignment->cigar_length,key_length,cigar_attr_none);
  PROFILE_STOP(GP_MATCHES_ALIGN_EXACT,PROFILE_LEVEL);
}
/*
 * Pseudo-Alignment
 */
void match_align_pseudoalignment(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    match_trace_t* const match_trace) {
  // Parameters
  const uint64_t text_offset_begin = alignment->alignment_tiles->text_begin_offset;
  const uint64_t key_length = pattern->key_length;
  uint8_t* const text = text_trace->text + text_offset_begin;
  // Configure match-trace
  match_trace->type = match_type_regular;
  match_trace->text_trace = text_trace;
  match_trace->text = text;
  match_trace->text_length = key_length;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment->match_text_offset = text_offset_begin;
  match_alignment->match_position = text_trace->position + text_offset_begin;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  match_alignment->score = alignment->distance_min_bound;
  // Set CIGAR
  matches_cigar_vector_append_match(matches->cigar_vector,
      &match_alignment->cigar_length,key_length,cigar_attr_none);
  match_alignment->effective_length = key_length;
  // Set distances
  match_trace->edit_distance = match_alignment->score;
  match_trace->event_distance = match_alignment->score;
  match_trace->swg_score = match_alignment->score;
  match_trace->error_quality = (float)SEQUENCE_QUALITIES_MAX;
  // Store matching text
  match_trace->match_scaffold = NULL;
  match_trace->text = text_trace->text + match_alignment->match_text_offset;
  match_trace->text_length = match_alignment->effective_length;
  // Add clipping to CIGAR
  match_aling_add_clipping(
      match_trace,matches->cigar_vector,
      pattern->clip_left,pattern->clip_right);
}
/*
 * Hamming Alignment (Only mismatches)
 */
void match_align_hamming(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    match_trace_t* const match_trace) {
  PROFILE_START(GP_MATCHES_ALIGN_HAMMING,PROFILE_LEVEL);
  // Parameters
  const uint64_t text_offset_begin = alignment->alignment_tiles->text_begin_offset;
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  uint8_t* const text = text_trace->text + text_offset_begin;
  // Configure match-trace
  match_trace->type = match_type_regular;
  match_trace->text_trace = text_trace;
  match_trace->text = text;
  match_trace->text_length = key_length;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment->match_text_offset = text_offset_begin;
  match_alignment->match_position = text_trace->position + text_offset_begin;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 0;
  // Hamming Alignment
  uint64_t i, mismatches;
  for (i=0,mismatches=0;i<key_length;++i) {
    const uint8_t candidate_enc = text[i];
    // Check Mismatch
    if (candidate_enc == ENC_DNA_CHAR_N || candidate_enc != key[i]) {
      ++mismatches;
      matches_cigar_vector_append_mismatch(matches->cigar_vector,
          &match_alignment->cigar_length,candidate_enc,cigar_attr_none);
    } else {
      matches_cigar_vector_append_match(matches->cigar_vector,
          &match_alignment->cigar_length,1,cigar_attr_none);
    }
  }
  match_alignment->effective_length = key_length;
  // Set distances
  match_trace->edit_distance = mismatches;
  match_trace->event_distance = mismatches;
  match_trace->swg_score = align_swg_score_cigar_excluding_clipping(
      &search_parameters->swg_penalties,matches->cigar_vector,
      match_alignment->cigar_offset,match_alignment->cigar_length);
  match_trace->error_quality = matches_cigar_compute_error_quality(
      matches->cigar_vector,match_alignment->cigar_offset,
      match_alignment->cigar_length,pattern->quality_mask,
      pattern->key_length);
  match_trace->match_scaffold = NULL;
  // Add clipping to CIGAR
  match_aling_add_clipping(
      match_trace,matches->cigar_vector,
      pattern->clip_left,pattern->clip_right);
  PROFILE_STOP(GP_MATCHES_ALIGN_HAMMING,PROFILE_LEVEL);
}
/*
 * Levenshtein
 */
void match_align_levenshtein(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  PROFILE_START(GP_MATCHES_ALIGN_LEVENSHTEIN,PROFILE_LEVEL);
  // Parameters
  uint8_t* const key = pattern->key;
  bpm_pattern_t* const bpm_pattern = &pattern->pattern_tiled.bpm_pattern;
  uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->text_length;
  const bool left_gap_alignment = (text_trace->strand == Forward);
  vector_t* const cigar_vector = matches->cigar_vector;
  // Configure match-trace
  match_trace->type = match_type_regular;
  match_trace->text_trace = text_trace;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  // Levenshtein Align
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment->match_position = text_trace->position;
  const uint64_t edit_bound = alignment->distance_min_bound;
  align_bpm_match(
      bpm_pattern,key,text,text_length,
      edit_bound,left_gap_alignment,
      match_alignment,cigar_vector,mm_allocator);
  // Set distances
  if (match_alignment->score==-1) {
    PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_SWG);
    match_trace->edit_distance = ALIGN_DISTANCE_INF;
    match_trace->event_distance = ALIGN_DISTANCE_INF;
    match_trace->swg_score = SWG_SCORE_MIN;
  } else {
    match_trace->edit_distance = match_alignment->score;
    match_trace->event_distance = matches_cigar_compute_event_distance(
        cigar_vector,match_alignment->cigar_offset,match_alignment->cigar_length);
    match_trace->swg_score = align_swg_score_cigar_excluding_clipping(
        &search_parameters->swg_penalties,cigar_vector,
        match_alignment->cigar_offset,match_alignment->cigar_length);
    match_trace->error_quality = matches_cigar_compute_error_quality(
        matches->cigar_vector,match_alignment->cigar_offset,
        match_alignment->cigar_length,pattern->quality_mask,
        pattern->key_length);
    // Store matching text
    match_trace->match_scaffold = NULL;
    match_alignment->match_text_offset = match_alignment->match_position - text_trace->position;
    match_trace->text = text + match_alignment->match_text_offset;
    match_trace->text_length = match_alignment->effective_length;
    // Add clipping to CIGAR
    match_aling_add_clipping(
        match_trace,matches->cigar_vector,
        pattern->clip_left,pattern->clip_right);
  }
  PROFILE_STOP(GP_MATCHES_ALIGN_LEVENSHTEIN,PROFILE_LEVEL);
}
/*
 * Smith-Waterman-Gotoh Alignment (Gap-affine)
 */
void match_align_smith_waterman_gotoh(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    match_scaffold_t* const match_scaffold,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  PROFILE_START(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
  // Configure match-trace
  match_trace->type = match_type_regular;
  match_trace->text_trace = text_trace;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  // Scaffold the alignment
  if (!search_parameters->alignment_force_full_swg) {
    PROFILE_PAUSE(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
    match_scaffold_adaptive(
        match_scaffold,pattern,text_trace,alignment,
        search_parameters->alignment_global_min_identity_nominal,
        search_parameters->alignment_scaffolding_min_matching_length_nominal,
        matches,mm_allocator);
    PROFILE_CONTINUE(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
  }
  match_trace->match_scaffold = match_scaffold;
  // Align SWG
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_align_swg(
      matches,search_parameters,pattern,text_trace,
      match_scaffold,false,match_trace,mm_allocator);
  if (match_alignment->score == SWG_SCORE_MIN) {
    match_trace->swg_score = SWG_SCORE_MIN;
    PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_SWG);
    PROFILE_STOP(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
    return;
  }
  // Normalize CIGAR (Adjust Position, Translate RL, ...)
  match_align_normalize(matches,match_trace,search_parameters);
  match_align_swg_compute_alignment_type(
      matches,match_trace,pattern,search_parameters);
  if (match_trace->type == match_type_local) {
    match_trace->swg_score = SWG_SCORE_MIN;
    PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_SWG);
    PROFILE_STOP(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
    return;
  }
  // Add clipping to CIGAR
  match_aling_add_clipping(
      match_trace,matches->cigar_vector,
      pattern->clip_left,pattern->clip_right);
  PROFILE_STOP(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
}
void match_align_smith_waterman_gotoh_local(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    alignment_t* const alignment,
    match_scaffold_t* const match_scaffold,
    match_trace_t* const match_trace,
    mm_allocator_t* const mm_allocator) {
  PROFILE_START(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
  // Configure match-trace
  match_trace->type = match_type_regular;
  match_trace->text_trace = text_trace;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->match_scaffold = match_scaffold;
  // Align SWG
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_align_swg(
      matches,search_parameters,pattern,text_trace,
      match_scaffold,true,match_trace,mm_allocator);
  if (match_alignment->score == SWG_SCORE_MIN) {
    match_trace->swg_score = SWG_SCORE_MIN;
    PROF_INC_COUNTER(GP_ALIGNED_DISCARDED_SWG);
    PROFILE_STOP(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
    return;
  }
  // Normalize CIGAR (Adjust Position, Translate RL, ...)
  match_align_normalize(matches,match_trace,search_parameters);
  match_align_swg_compute_alignment_type(
      matches,match_trace,pattern,search_parameters);
  if (match_trace->type == match_type_local) {
    // Compute Local Alignment
    match_align_swg_local_alignment(
        matches,search_parameters,pattern,
        text_trace,match_scaffold,match_trace);
  } else {
    // Add clipping to CIGAR
    match_aling_add_clipping(
        match_trace,matches->cigar_vector,
        pattern->clip_left,pattern->clip_right);
  }
  PROFILE_STOP(GP_MATCHES_ALIGN_SWG,PROFILE_LEVEL);
}

