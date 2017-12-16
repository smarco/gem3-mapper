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
 * DESCRIPTION:
 *   Filtering region module provides functions to produce a full-alignment
 *   of a filtering-region against its corresponding text-region
 */

#include "align/alignment.h"
#include "filtering/region/filtering_region_align.h"
#include "matches/align/match_align.h"
#include "io/output_map.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION  GEM_DEEP_DEBUG

/*
 * Debug
 */
void filtering_region_align_debug_prologue(
    filtering_region_t* const filtering_region) {
  PROF_INC_COUNTER(GP_ALIGNED_REGIONS);
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (region_align)\n");
    tab_global_inc();
    filtering_region_print(gem_log_get_stream(),filtering_region,false,true,true);
  }
}
void filtering_region_align_debug_epilogue(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches,
    match_trace_t* const match_trace) {
  // Check alignment result
  if (match_trace->swg_score == SWG_SCORE_MIN) {
    PROF_INC_COUNTER(GP_ALIGNED_DISCARDED);
  } else {
    PROF_INC_COUNTER(GP_ALIGNED_ACCEPTED);
  }
  // Deep debug
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    if (match_trace->swg_score == SWG_SCORE_MIN) {
      tab_fprintf(gem_log_get_stream(),
          "=> Region NOT-ALIGNED (distance=%lu,swg_score=%ld)\n",
          match_trace->event_distance,match_trace->swg_score);
      tab_global_dec();
    } else {
      // Text Candidate
      const text_trace_t* const text_trace = &filtering_region->text_trace;
      uint8_t* const text = text_trace->text;
      // Print debug info
      tab_fprintf(gem_log_get_stream(),"=> Region ALIGNED (distance=%lu,swg_score=%ld)\n",
          match_trace->event_distance,match_trace->swg_score);
      tab_global_inc();
      output_map_alignment_pretty(gem_log_get_stream(),match_trace,matches,pattern->key,pattern->key_length,
          text+(match_trace->match_alignment.match_position - filtering_region->text_begin_position),
          match_trace->match_alignment.effective_length,filtering_candidates->mm_allocator);
      tab_global_dec();
      tab_global_dec();
    }
  }
}
/*
 * Region (Re)Align by clone previous
 */
void filtering_region_align_clone(
    match_trace_t* const match_trace_src,
    match_trace_t* const match_trace_dst,
    filtering_region_t* const filtering_region_dst) {
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (region_align_clone)\n");
    tab_global_inc();
  }
  match_trace_dst->type = match_trace_src->type;
  // Clone match-trace (Match Text (Reference)
  match_trace_dst->text_trace = match_trace_src->text_trace;
  match_trace_dst->text = match_trace_src->text;
  match_trace_dst->text_length = match_trace_src->text_length;
  // Clone match-trace (Match)
  match_trace_dst->sequence_name = NULL;
  match_trace_dst->text_position = UINT64_MAX;
  // Clone match-trace (Score)
  match_trace_dst->event_distance = match_trace_src->event_distance;
  match_trace_dst->edit_distance = match_trace_src->edit_distance;
  match_trace_dst->swg_score = match_trace_src->swg_score;
  match_trace_dst->error_quality = match_trace_src->error_quality;
  // Clone match-alignment
  match_alignment_t* const match_alignment_dst = &match_trace_dst->match_alignment;
  match_alignment_t* const match_alignment_src = &match_trace_src->match_alignment;
  match_alignment_dst->match_text_offset = match_alignment_src->match_text_offset;
  match_alignment_dst->match_position =
      filtering_region_dst->text_begin_position + match_alignment_src->match_text_offset;
  match_alignment_dst->cigar_offset = match_alignment_src->cigar_offset;
  match_alignment_dst->cigar_length = match_alignment_src->cigar_length;
  match_alignment_dst->effective_length = match_alignment_src->effective_length;
  match_alignment_dst->score = match_alignment_src->score;
  match_trace_dst->match_scaffold = match_trace_src->match_scaffold; // Supporting Scaffolding
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"=> Region CLONED (distance=%lu,swg_score=%ld)\n",
        match_trace_dst->event_distance,match_trace_dst->swg_score);
    tab_global_dec();
  }
}
/*
 * Region (Re)Align
 */
bool filtering_region_align(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches,
    match_trace_t* const match_trace) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  const match_alignment_model_t match_alignment_model = search_parameters->match_alignment_model;
  text_trace_t* const text_trace = &filtering_region->text_trace;
  alignment_t* const alignment = &filtering_region->alignment;
  mm_allocator_t* const mm_allocator = filtering_candidates->mm_allocator;
  // DEBUG
  filtering_region_align_debug_prologue(filtering_region);
  // Select Model
  if (filtering_region->alignment.distance_min_bound > 0 ||
      filtering_region->alignment.num_tiles > 1 ||
      filtering_region->key_trimmed ||
      pattern->run_length) {
    PROF_INC_COUNTER(GP_ALIGNED_INEXACT);
    // Select alignment model
    switch (match_alignment_model) {
      case match_alignment_model_none:
        // Hamming Align
        match_align_pseudoalignment(
            matches,search_parameters,pattern,
            text_trace,alignment,match_trace);
        break;
      case match_alignment_model_hamming: {
        // Hamming Align
        match_align_hamming(
            matches,search_parameters,pattern,
            text_trace,alignment,match_trace);
        break;
      }
      case match_alignment_model_levenshtein: {
        // Levenshtein Align
        match_align_levenshtein(
            matches,search_parameters,pattern,
            text_trace,alignment,match_trace,mm_allocator);
        break;
      }
      case match_alignment_model_gap_affine: {
        // Gap-affine Align
        match_align_smith_waterman_gotoh(
            matches,search_parameters,pattern,
            &filtering_region->text_trace,
            &filtering_region->alignment,
            &filtering_region->match_scaffold,
            match_trace,mm_allocator);
        filtering_candidates->total_candidates_realigned_swg++;
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
    PROF_ADD_COUNTER(GP_ALIGNED_REGIONS_LENGTH,text_trace->text_length);
  } else {
    PROF_INC_COUNTER(GP_ALIGNED_EXACT);
    // Add exact match
    match_align_exact(
        matches,search_parameters,pattern,
        text_trace,alignment,match_trace);
    // Add input-clipping to CIGAR
    match_aling_add_clipping(
        match_trace,matches->cigar_vector,
        pattern->clip_left,pattern->clip_right);
  }
  // DEBUG
  filtering_region_align_debug_epilogue(filtering_candidates,
      filtering_region,pattern,matches,match_trace);
  // Return alignment result
  return (match_trace->swg_score != SWG_SCORE_MIN);
}
bool filtering_region_align_local(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches,
    match_trace_t* const match_trace) {
  // DEBUG
  filtering_region_align_debug_prologue(filtering_region);
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  text_trace_t* const text_trace = &filtering_region->text_trace;
  alignment_t* const alignment = &filtering_region->alignment;
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  mm_allocator_t* const mm_allocator = filtering_candidates->mm_allocator;
  // Gap-affine Align
  match_align_smith_waterman_gotoh_local(
      matches,search_parameters,pattern,
      text_trace,alignment,match_scaffold,
      match_trace,mm_allocator);
  // DEBUG
  filtering_region_align_debug_epilogue(filtering_candidates,
      filtering_region,pattern,matches,match_trace);
  // Return alignment result
  return (match_trace->swg_score != SWG_SCORE_MIN);
}
