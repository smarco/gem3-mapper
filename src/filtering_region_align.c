/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region_align.h"
#include "filtering_region_align_configure.h"
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
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (region_align_clone)\n");
    tab_global_inc();
  }
  // Clone match-trace (Match Text (Reference)
  match_trace_dst->text_trace_offset = match_trace_src->text_trace_offset;
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
  match_alignment_dst->match_text_offset = match_alignment_src->match_text_offset;
  match_alignment_dst->match_position = filtering_region_dst->begin_position + match_alignment_dst->match_text_offset;
  match_alignment_dst->cigar_offset = match_alignment_src->cigar_offset;
  match_alignment_dst->cigar_length = match_alignment_src->cigar_length;
  match_alignment_dst->effective_length = match_alignment_src->effective_length;
  match_alignment_dst->score = match_alignment_src->score;
#ifdef GEM_DEBUG
  match_trace_dst->match_scaffold = match_trace_src->match_scaffold; // Supporting Scaffolding
#endif
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"=> Region CLONED (distance=%lu,swg_score=%ld)\n",
        match_trace_dst->distance,match_trace_dst->swg_score);
    tab_global_dec();
  }
}
/*
 * Adjust distance bound by scaffolding
 */
void filtering_region_align_adjust_distance_by_scaffolding(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const as_parameters_t* const as_parameters,pattern_t* const pattern,
    matches_t* const matches,text_collection_t* const text_collection,
    mm_stack_t* const mm_stack) {
  // Retrieve Candidate (if needed)
  filtering_region_retrieve_text(filtering_region,archive_text,text_collection,mm_stack);
  // Text Candidate
  const text_trace_t* const text_trace =
      text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint64_t text_length = filtering_region->end_position-filtering_region->begin_position;
  uint8_t* const text = text_trace->text;
  // Compile pattern (if trimmed)
  if (filtering_region->key_trimmed && filtering_region->bpm_pattern_trimmed==NULL) {
    pattern_trimmed_init(pattern,
        &filtering_region->bpm_pattern_trimmed,&filtering_region->bpm_pattern_trimmed_tiles,
        filtering_region->key_trim_left,filtering_region->key_trim_right,mm_stack);
  }
  // Configure Alignment
  match_align_input_t align_input;
  match_align_parameters_t align_parameters;
  const strand_t position_strand = archive_text_get_position_strand(archive_text,filtering_region->begin_position);
  const bool left_gap_alignment = (position_strand==Forward);
  filtering_region_align_configure_scaffold(&align_input,&align_parameters,
      filtering_region,as_parameters,pattern,text,text_length,left_gap_alignment);
  // Scaffold alignment & re-check distance limits
  match_scaffold_adaptive(matches,&align_input,&align_parameters,&filtering_region->match_scaffold,mm_stack);
  // Adjust distance using scaffold
  filtering_region->align_distance_min_bound = filtering_region->match_scaffold.match_alignment.score;
}
/*
 * Region (Re)Align
 */
void filtering_region_align_exact(
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,pattern_t* const pattern,
    matches_t* const matches,match_trace_t* const match_trace) {
  PROF_INC_COUNTER(GP_ALIGNED_EXACT);
  // Parameters
  match_align_input_t align_input;
  match_align_parameters_t align_parameters;
  // Configure Alignment
  filtering_region_align_configure_exact(&align_input,&align_parameters,
      filtering_region,as_parameters,pattern,emulated_rc_search);
  // Add exact match
  match_trace->match_alignment.score = filtering_region->align_distance;
  match_align_exact(matches,match_trace,&align_input,&align_parameters);
  filtering_region->status = filtering_region_aligned; // Set status
}
void filtering_region_align_inexact(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const text_collection_t* const text_collection,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,pattern_t* const pattern,
    matches_t* const matches,match_trace_t* const match_trace,
    mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_ALIGNED_INEXACT);
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  match_align_input_t align_input;
  match_align_parameters_t align_parameters;
  // Text Candidate
  const text_trace_t* const text_trace =
      text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint64_t text_length = filtering_region->end_position-filtering_region->begin_position;
  uint8_t* const text = text_trace->text;
  // Select alignment model
  switch (alignment_model) {
    case alignment_model_hamming: {
      // Configure Alignment
      filtering_region_align_configure_hamming(&align_input,&align_parameters,
          filtering_region,as_parameters,pattern,text,text_length,emulated_rc_search);
      // Hamming Align
      match_align_hamming(matches,match_trace,&align_input,&align_parameters);
      filtering_region->status = filtering_region_aligned;
      break;
    }
    case alignment_model_levenshtein: {
      // Compile pattern (if trimmed)
      if (filtering_region->key_trimmed && filtering_region->bpm_pattern_trimmed==NULL) {
        pattern_trimmed_init(pattern,
            &filtering_region->bpm_pattern_trimmed,&filtering_region->bpm_pattern_trimmed_tiles,
            filtering_region->key_trim_left,filtering_region->key_trim_right,mm_stack);
      }
      // Configure Alignment
      const strand_t position_strand = archive_text_get_position_strand(archive_text,filtering_region->begin_position);
      const bool left_gap_alignment = (position_strand==Forward);
      filtering_region_align_configure_levenshtein(
          &align_input,&align_parameters,filtering_region,as_parameters,
          pattern,text,text_length,emulated_rc_search,left_gap_alignment);
      // Levenshtein Align
      match_align_levenshtein(matches,match_trace,&align_input,&align_parameters,mm_stack);
      break;
    }
    case alignment_model_gap_affine: {
      // Compile pattern (if trimmed)
      if (filtering_region->key_trimmed && filtering_region->bpm_pattern_trimmed==NULL) {
        pattern_trimmed_init(pattern,
            &filtering_region->bpm_pattern_trimmed,&filtering_region->bpm_pattern_trimmed_tiles,
            filtering_region->key_trim_left,filtering_region->key_trim_right,mm_stack);
      }
      // Configure Alignment
      const strand_t position_strand = archive_text_get_position_strand(archive_text,filtering_region->begin_position);
      const bool left_gap_alignment = (position_strand==Forward);
      filtering_region_align_configure_swg(
          &align_input,&align_parameters,filtering_region,as_parameters,
          pattern,text,text_length,emulated_rc_search,left_gap_alignment);
      // Gap-affine Align
      match_align_smith_waterman_gotoh(matches,match_trace,
          &align_input,&align_parameters,&filtering_region->match_scaffold,mm_stack);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_ADD_COUNTER(GP_ALIGNED_REGIONS_LENGTH,align_input.text_offset_end-align_input.text_offset_begin);
}
bool filtering_region_align(
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
  // Select Model
  if (!filtering_region->key_trimmed &&
      (filtering_region->align_distance==0 || alignment_model==alignment_model_none)) {
    filtering_region_align_exact(filtering_region,as_parameters,
        emulated_rc_search,pattern,matches,match_trace);
  } else {
    filtering_region_align_inexact(filtering_region,archive_text,text_collection,
        as_parameters,emulated_rc_search,pattern,matches,match_trace,mm_stack);
  }
  // Check (re)alignment result
  if (match_trace->distance!=ALIGN_DISTANCE_INF && match_trace->swg_score >= 0) {
    // Assign Aligned-status & store offset
    filtering_region->status = filtering_region_aligned;
    // DEBUG
    gem_cond_debug_block(DEBUG_FILTERING_REGION) {
      // Text Candidate
      const text_trace_t* const text_trace =
          text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
      uint8_t* const text = text_trace->text;
      // Print debug info
      tab_fprintf(gem_log_get_stream(),"=> Region ALIGNED (distance=%lu,swg_score=%ld)\n",
          match_trace->distance,match_trace->swg_score);
      tab_global_inc();
      output_map_alignment_pretty(gem_log_get_stream(),match_trace,matches,pattern->key,pattern->key_length,
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
bool filtering_region_align_unbounded(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const text_collection_t* const text_collection,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,pattern_t* const pattern,matches_t* const matches,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_ALIGNED_REGIONS);
//  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
//    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (align_unbounded)\n");
//    tab_global_inc();
//    filtering_region_print(gem_log_get_stream(),filtering_region,text_collection,false);
//  }
//  // Parameters
//  search_parameters_t* const search_parameters = as_parameters->search_parameters;
//  const alignment_model_t alignment_model = search_parameters->alignment_model;
//  uint8_t* const key = pattern->key;
//  const uint64_t key_length = pattern->key_length;
//  // Text Candidate
//  const text_trace_t* const text_trace =
//      text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
//  uint64_t text_length = filtering_region->end_position-filtering_region->begin_position;
//  uint8_t* text = text_trace->text;
//  // Select alignment model
//  switch (alignment_model) {
//    case alignment_model_hamming:
//    case alignment_model_levenshtein:
//      return false; // TODO
//      break;
//    case alignment_model_gap_affine: {
//      // Configure Alignment
//      match_align_input_t align_input;
//      match_align_parameters_t align_parameters;
//      const strand_t strand = archive_text_get_position_strand(archive_text,filtering_region->begin_position);
//      const bool left_gap_alignment = (strand==Forward);
//      filtering_region_align_configure_local_swg(
//          &align_input,&align_parameters,filtering_region,as_parameters,
//          pattern,text,text_length,emulated_rc_search,left_gap_alignment);
//      PROF_ADD_COUNTER(GP_ALIGNED_REGIONS_LENGTH,align_input.text_offset_end-align_input.text_offset_begin);
//      // Smith-Waterman-Gotoh Alignment (Gap-affine)
//      match_align_local_smith_waterman_gotoh(matches,match_trace,
//          &align_input,&align_parameters,&filtering_region->match_scaffold,mm_stack);
//      break;
//    }
//    default:
//      GEM_INVALID_CASE();
//      break;
//  }
//  // Check (re)alignment result
//  if (match_trace->distance!=ALIGN_DISTANCE_INF) {
//    // DEBUG
//    gem_cond_debug_block(DEBUG_FILTERING_REGION) {
//      tab_fprintf(gem_log_get_stream(),"=> Region ALIGNED (distance=%lu,swg_score=%ld)\n",
//          match_trace->distance,match_trace->swg_score);
//      tab_global_inc();
//      output_map_alignment_pretty(gem_log_get_stream(),match_trace,matches,key,key_length,
//          text+(match_trace->match_alignment.match_position - filtering_region->begin_position),
//          match_trace->match_alignment.effective_length,mm_stack);
//      tab_global_dec();
//      tab_global_dec();
//    }
//    return true;
//  } else {
//    // DEBUG
//    gem_cond_debug_block(DEBUG_FILTERING_REGION) {
//      tab_fprintf(gem_log_get_stream(),
//          "=> Region SUBDOMINANT (distance=%lu,swg_score=%ld)\n",
//          match_trace->distance,match_trace->swg_score);
//      tab_global_dec();
//    }
//    return false;
//  }
  return false;
}
