/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region_align.h"
#include "align.h"
#include "match_align_local.h"
#include "output_map.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION  GEM_DEEP_DEBUG

/*
 * Region (Re)Align by clone previous
 */
void filtering_region_align_clone(
    match_trace_t* const match_trace_src,match_trace_t* const match_trace_dst,
    filtering_region_t* const filtering_region_dst,text_collection_t* const text_collection) {
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (region_align)\n");
    tab_global_inc();
  }
  // Clone match-trace (Match Text (Reference)
  match_trace_dst->trace_offset = match_trace_src->trace_offset;
  match_trace_dst->text = match_trace_src->text;
  match_trace_dst->text_length = match_trace_src->text_length;
  // Clone match-trace (Match)
  match_trace_dst->sequence_name = NULL;
  match_trace_dst->text_position = UINT64_MAX;
  match_trace_dst->emulated_rc_search = match_trace_src->emulated_rc_search;
  // Clone match-trace (Score)
  match_trace_dst->distance = match_trace_src->distance;
  match_trace_dst->edit_distance = match_trace_src->edit_distance;
  match_trace_dst->swg_score = match_trace_src->swg_score;
  // Clone match-alignment
  match_alignment_t* const match_alignment_dst = &match_trace_dst->match_alignment;
  match_alignment_t* const match_alignment_src = &match_trace_src->match_alignment;
  match_alignment_dst->match_position = filtering_region_dst->begin_position + match_alignment_src->match_offet;
  match_alignment_dst->cigar_offset = match_alignment_src->cigar_offset;
  match_alignment_dst->cigar_length = match_alignment_src->cigar_length;
  match_alignment_dst->effective_length = match_alignment_src->effective_length;
  match_alignment_dst->score = match_alignment_src->score;
#ifdef GEM_DEBUG
  match_trace_dst->match_scaffold = match_trace_src->match_scaffold;    // Supporting Scaffolding
#endif
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"=> Region ALIGNED-CACHED (distance=%lu,swg_score=%ld)\n",
        match_trace_dst->distance,match_trace_dst->swg_score);
    tab_global_dec();
  }
}

/*
 * Region (Re)Align
 */
GEM_INLINE bool filtering_region_align(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const text_collection_t* const text_collection,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,pattern_t* const pattern,matches_t* const matches,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_ALIGNED_REGIONS);
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (region_align)\n");
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  bool* const allowed_enc = search_parameters->allowed_enc;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Text
  uint64_t text_length;
  uint8_t* text;
  // Select Model
  if (filtering_region->align_distance==0 || alignment_model==alignment_model_none) {
    // Add exact match
    PROF_INC_COUNTER(GP_ALIGNED_EXACT);
    match_align_input_t align_input = {
        .key_length = key_length,
        .text_position = filtering_region->begin_position,
        .text_trace_offset = filtering_region->text_trace_offset,
        .text = NULL, // Explicitly nullify (no need to use it afterwards)
        .text_offset_begin = filtering_region->align_match_end_column+1 - key_length,
        .text_offset_end = filtering_region->align_match_end_column+1,
    };
    match_align_parameters_t align_parameters = {
        .emulated_rc_search = emulated_rc_search,
        .swg_penalties = swg_penalties,
    };
    match_trace->match_alignment.score = filtering_region->align_distance;
    match_align_exact(matches,match_trace,&align_input,&align_parameters);
    filtering_region->status = filtering_region_aligned; // Set status
    return true; // OK
  } else {
    PROF_INC_COUNTER(GP_ALIGNED_INEXACT);
    // Text Candidate
    const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
    text_length = filtering_region->end_position-filtering_region->begin_position;
    text = text_trace->text;
    // Select alignment model
    switch (alignment_model) {
      case alignment_model_hamming: {
        match_align_input_t align_input = {
            .key = key,
            .key_length = key_length,
            .text_trace_offset = filtering_region->text_trace_offset,
            .text_position = filtering_region->begin_position + filtering_region->base_position_offset, // Base position
            .text = text,
            .text_offset_begin = filtering_region->align_match_begin_column,
            .text_offset_end = filtering_region->align_match_begin_column + key_length,
        };
        match_align_parameters_t align_parameters = {
            .emulated_rc_search = emulated_rc_search,
            .allowed_enc = allowed_enc,
        };
        match_align_hamming(matches,match_trace,&align_input,&align_parameters);
        filtering_region->status = filtering_region_aligned;
        return true; // OK
      }
      case alignment_model_levenshtein: {
        const uint64_t match_effective_range = key_length + 2*filtering_region->align_distance;
        const uint64_t match_begin_column = BOUNDED_SUBTRACTION(filtering_region->align_match_end_column,match_effective_range,0);
        const uint64_t match_end_column = BOUNDED_ADDITION(filtering_region->align_match_end_column,filtering_region->align_distance,text_length);
        PROF_ADD_COUNTER(GP_ALIGNED_REGIONS_LENGTH,(match_end_column-match_end_column));
        match_align_input_t align_input = {
            .key = key,
            .key_length = key_length,
            .bpm_pattern = &pattern->bpm_pattern,
            .text_trace_offset = filtering_region->text_trace_offset,
            .text_position = filtering_region->begin_position,
            .text = text,
            .text_offset_begin = match_begin_column,
            .text_offset_end = match_end_column,
            .text_length = text_length,
        };
        match_align_parameters_t align_parameters = {
            .emulated_rc_search = emulated_rc_search,
            .max_error = filtering_region->align_distance,
            .left_gap_alignment = archive_text_get_position_strand(archive_text,align_input.text_position)==Forward,
            .swg_penalties = swg_penalties,
        };
        match_align_levenshtein(matches,match_trace,&align_input,&align_parameters,mm_stack);
        break;
      }
      case alignment_model_gap_affine: {
        // Adjust alignment boundaries (to allow optimization)
        const uint64_t align_distance_bound = filtering_region->align_distance;
        const uint64_t align_margin = key_length + filtering_region->align_distance;
        const uint64_t align_match_begin_column = BOUNDED_SUBTRACTION(filtering_region->align_match_end_column,align_margin,0);
        const uint64_t align_match_end_column = filtering_region->align_match_end_column;
        const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
        const uint64_t match_effective_range = key_length + 2*filtering_region->align_distance + max_bandwidth;
        const uint64_t text_offset_begin = BOUNDED_SUBTRACTION(filtering_region->align_match_end_column,match_effective_range,0);
        const uint64_t text_offset_end = BOUNDED_ADDITION(filtering_region->align_match_end_column,max_bandwidth,text_length-1);
        PROF_ADD_COUNTER(GP_ALIGNED_REGIONS_LENGTH,text_offset_end-text_offset_begin);
        match_align_input_t align_input = {
            .key = key,
            .key_length = key_length,
            .bpm_pattern = &pattern->bpm_pattern,
            .text_trace_offset = filtering_region->text_trace_offset,
            .text_position = filtering_region->begin_position,
            .text = text,
            .text_length = text_length,
            .text_offset_begin = text_offset_begin,
            .text_offset_end = text_offset_end,
            .align_distance_bound = align_distance_bound,
            .align_match_begin_column = align_match_begin_column,
            .align_match_end_column = align_match_end_column,
        };
        match_align_parameters_t align_parameters = {
            .emulated_rc_search = emulated_rc_search,
            .max_error = pattern->max_effective_filtering_error,
            .max_bandwidth = max_bandwidth,
            .left_gap_alignment = archive_text_get_position_strand(archive_text,align_input.text_position)==Forward,
            .min_identity = as_parameters->alignment_min_identity_nominal,
            .scaffolding = search_parameters->alignment_scaffolding,
            .scaffolding_min_coverage = as_parameters->alignment_scaffolding_min_coverage_nominal,
            .scaffolding_matching_min_length = as_parameters->alignment_scaffolding_min_matching_length_nominal,
            .scaffolding_homopolymer_min_context = as_parameters->alignment_scaffolding_homopolymer_min_context_nominal,
            .allowed_enc = allowed_enc,
            .swg_penalties = swg_penalties,
            .swg_threshold = as_parameters->swg_threshold_nominal,
            .cigar_curation = search_parameters->cigar_curation,
            .cigar_curation_min_end_context = as_parameters->cigar_curation_min_end_context_nominal,
        };
        // Smith-Waterman-Gotoh Alignment (Gap-affine)
        match_align_smith_waterman_gotoh(matches,match_trace,
            &align_input,&align_parameters,&filtering_region->match_scaffold,mm_stack);
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Check (re)alignment result
  if (match_trace->distance!=ALIGN_DISTANCE_INF && match_trace->swg_score >= 0) { // TODO Check
    // Assign Aligned-status & store offset
    filtering_region->status = filtering_region_aligned;
    // DEBUG
    gem_cond_debug_block(DEBUG_FILTERING_REGION) {
      tab_fprintf(gem_log_get_stream(),"=> Region ALIGNED (distance=%lu,swg_score=%ld)\n",
          match_trace->distance,match_trace->swg_score);
      tab_global_inc();
      output_map_alignment_pretty(gem_log_get_stream(),match_trace,matches,key,key_length,
          text+(match_trace->match_alignment.match_position - filtering_region->begin_position),
          match_trace->match_alignment.effective_length,mm_stack);
      tab_global_dec();
      tab_global_dec();
    }
    return true; // OK
  } else {
    filtering_region->status = filtering_region_aligned_subdominant;
    // DEBUG
    gem_cond_debug_block(DEBUG_FILTERING_REGION) {
      tab_fprintf(gem_log_get_stream(),"=> Region SUBDOMINANT (distance=%lu,swg_score=%ld)\n",
          match_trace->distance,match_trace->swg_score);
      tab_global_dec();
    }
    return false; // Discarded
  }
}
GEM_INLINE bool filtering_region_align_unbounded(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const text_collection_t* const text_collection,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,pattern_t* const pattern,matches_t* const matches,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_ALIGNED_REGIONS);
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (align_unbounded)\n");
    tab_global_inc();
    filtering_region_print(gem_log_get_stream(),filtering_region,text_collection,false);
  }
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  bool* const allowed_enc = search_parameters->allowed_enc;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Text
  uint64_t text_length;
  uint8_t* text;
  // Select alignment model
  switch (alignment_model) {
    case alignment_model_hamming:
    case alignment_model_levenshtein:
      // TODO
      return false;
      break;
    case alignment_model_gap_affine: {
      // Text Candidate
      const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
      text_length = filtering_region->end_position-filtering_region->begin_position;
      text = text_trace->text;
      // Adjust alignment boundaries (to allow optimization)
      const uint64_t text_offset_begin = 0;
      const uint64_t text_offset_end = text_length-1;
      PROF_ADD_COUNTER(GP_ALIGNED_REGIONS_LENGTH,text_offset_end-text_offset_begin);
      match_align_input_t align_input = {
          .key = key,
          .key_length = key_length,
          .bpm_pattern = &pattern->bpm_pattern,
          .text_trace_offset = filtering_region->text_trace_offset,
          .text_position = filtering_region->begin_position,
          .text = text,
          .text_length = text_length,
          .text_offset_begin = text_offset_begin,
          .text_offset_end = text_offset_end,
      };
      match_align_parameters_t align_parameters = {
          .emulated_rc_search = emulated_rc_search,
          .max_error = BOUNDED_SUBTRACTION(key_length,as_parameters->alignment_min_identity_nominal,0),
          .max_bandwidth = pattern->max_effective_bandwidth,
          .left_gap_alignment = archive_text_get_position_strand(archive_text,align_input.text_position)==Forward,
          .min_identity = as_parameters->alignment_min_identity_nominal,
          .scaffolding = true,
          .scaffolding_min_coverage = as_parameters->alignment_scaffolding_min_coverage_nominal,
          .scaffolding_matching_min_length = as_parameters->alignment_scaffolding_min_matching_length_nominal,
          .scaffolding_homopolymer_min_context = as_parameters->alignment_scaffolding_homopolymer_min_context_nominal,
          .allowed_enc = allowed_enc,
          .swg_penalties = swg_penalties,
          .swg_threshold = as_parameters->swg_threshold_nominal,
          .cigar_curation = search_parameters->cigar_curation,
          .cigar_curation_min_end_context = as_parameters->cigar_curation_min_end_context_nominal,
      };
      // Smith-Waterman-Gotoh Alignment (Gap-affine)
      match_align_local_smith_waterman_gotoh(matches,match_trace,
          &align_input,&align_parameters,&filtering_region->match_scaffold,mm_stack);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Check (re)alignment result
  if (match_trace->distance!=ALIGN_DISTANCE_INF) {
    // DEBUG
    gem_cond_debug_block(DEBUG_FILTERING_REGION) {
      tab_fprintf(gem_log_get_stream(),"=> Region ALIGNED (distance=%lu,swg_score=%ld)\n",
          match_trace->distance,match_trace->swg_score);
      tab_global_inc();
      output_map_alignment_pretty(gem_log_get_stream(),match_trace,matches,key,key_length,
          text+(match_trace->match_alignment.match_position - filtering_region->begin_position),
          match_trace->match_alignment.effective_length,mm_stack);
      tab_global_dec();
      tab_global_dec();
    }
    return true;
  } else {
    // DEBUG
    gem_cond_debug_block(DEBUG_FILTERING_REGION) {
      tab_fprintf(gem_log_get_stream(),"=> Region SUBDOMINANT (distance=%lu,swg_score=%ld)\n",
          match_trace->distance,match_trace->swg_score);
      tab_global_dec();
    }
    return false;
  }
}
