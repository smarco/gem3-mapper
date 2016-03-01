/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_region_align.h"
#include "filtering/filtering_region_align_configure.h"
#include "align/align.h"
#include "io/output_map.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION  GEM_DEEP_DEBUG

/*
 * Debug
 */
void filtering_region_align_debug(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const bool aligned,
    matches_t* const matches,
    match_trace_t* const match_trace) {
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    if (!aligned) {
      tab_fprintf(gem_log_get_stream(),
          "=> Region NOT-ALIGNED (distance=%lu,swg_score=%ld)\n",
          match_trace->distance,match_trace->swg_score);
      tab_global_dec();
    } else {
      // Text Candidate
      const text_trace_t* const text_trace = text_collection_get_trace(
          filtering_candidates->text_collection,filtering_region->text_trace_offset);
      uint8_t* const text = text_trace->text;
      // Print debug info
      tab_fprintf(gem_log_get_stream(),"=> Region ALIGNED (distance=%lu,swg_score=%ld)\n",
          match_trace->distance,match_trace->swg_score);
      tab_global_inc();
      output_map_alignment_pretty(gem_log_get_stream(),match_trace,matches,pattern->key,pattern->key_length,
          text+(match_trace->match_alignment.match_position - filtering_region->text_begin_position),
          match_trace->match_alignment.effective_length,filtering_candidates->mm_stack);
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
    filtering_region_t* const filtering_region_dst,
    const uint64_t run_length) {
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (region_align_clone)\n");
    tab_global_inc();
  }
  match_trace_dst->type = match_trace_src->type;
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
  match_alignment_dst->match_position =
      filtering_region_dst->text_begin_position + match_alignment_dst->match_text_offset;
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
 * Region (Re)Align
 */
void filtering_region_align_exact(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    matches_t* const matches,
    match_trace_t* const match_trace) {
  PROF_INC_COUNTER(GP_ALIGNED_EXACT);
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  match_align_input_t align_input;
  match_align_parameters_t align_parameters;
  // Configure Alignment
  filtering_region_align_configure_exact(
      &align_input,&align_parameters,filtering_region,
      search_parameters,pattern,emulated_rc_search);
  // Add exact match
  match_align_exact(matches,match_trace,&align_input,&align_parameters);
  filtering_region->status = filtering_region_aligned; // Set status
}
void filtering_region_align_inexact(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    const bool local_alignment,
    matches_t* const matches,
    match_trace_t* const match_trace) {
  PROF_INC_COUNTER(GP_ALIGNED_INEXACT);
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  match_align_input_t align_input;
  match_align_parameters_t align_parameters;
  // Select alignment model
  switch (alignment_model) {
    case alignment_model_none:
      GEM_NOT_IMPLEMENTED();
      break;
    case alignment_model_hamming: {
      // Configure Alignment
      filtering_region_align_configure_hamming(
          &align_input,&align_parameters,filtering_region,
          search_parameters,pattern,text_trace,emulated_rc_search);
      // Hamming Align
      match_align_hamming(matches,match_trace,&align_input,&align_parameters);
      filtering_region->status = filtering_region_aligned;
      break;
    }
    case alignment_model_levenshtein: {
      // Configure Alignment
      const strand_t position_strand =
          archive_text_get_position_strand(archive_text,filtering_region->text_begin_position);
      const bool left_gap_alignment = (position_strand==Forward);
      filtering_region_align_configure_levenshtein(
          &align_input,&align_parameters,filtering_region,search_parameters,
          pattern,text_trace,emulated_rc_search,left_gap_alignment,mm_stack);
      // Levenshtein Align
      match_align_levenshtein(matches,match_trace,&align_input,&align_parameters,mm_stack);
      break;
    }
    case alignment_model_gap_affine: {
      // Configure Alignment
      const strand_t position_strand =
          archive_text_get_position_strand(archive_text,filtering_region->text_begin_position);
      const bool left_gap_alignment = (position_strand==Forward);
      filtering_region_align_configure_swg(
          &align_input,&align_parameters,filtering_region,search_parameters,pattern,
          text_trace,emulated_rc_search,left_gap_alignment,local_alignment,mm_stack);
      // Gap-affine Align
      match_align_smith_waterman_gotoh(matches,match_trace,
          &align_input,&align_parameters,&filtering_region->match_scaffold,mm_stack);
      //      // DEBUG
      //      fprintf(stderr,"[GEM]>Text-Alignment\n");
      //      match_alignment_print_pretty(stderr,&match_trace->match_alignment,
      //          matches->cigar_vector,pattern->key,pattern->key_length,
      //          text_trace->text+match_trace->match_alignment.match_text_offset,
      //          text_trace->text_length,mm_stack);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_ADD_COUNTER(GP_ALIGNED_REGIONS_LENGTH,align_input.text_length);
}
bool filtering_region_align(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    const bool local_alignment,
    matches_t* const matches,
    match_trace_t* const match_trace) {
  // DEBUG
  PROF_INC_COUNTER(GP_ALIGNED_REGIONS);
  gem_cond_debug_block(DEBUG_FILTERING_REGION) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Region (region_align)\n");
    tab_global_inc();
  }
  // Select Model
  if (pattern->run_length ||
      filtering_region->key_trimmed ||
      filtering_region->region_alignment.num_tiles>1 ||
      filtering_region->region_alignment.distance_min_bound > 0) {
    filtering_region_align_inexact(filtering_candidates,filtering_region,
        pattern,emulated_rc_search,local_alignment,matches,match_trace);
  } else {
    filtering_region_align_exact(filtering_candidates,filtering_region,
        pattern,emulated_rc_search,matches,match_trace);
  }
  // Check (re)alignment result
  if (match_trace->distance==ALIGN_DISTANCE_INF || match_trace->swg_score < 0) {
    filtering_region_align_debug(filtering_candidates,
        filtering_region,pattern,false,matches,match_trace); // DEBUG
    return false; // Discarded
  } else {
    // Assign Aligned-status & store offset
    filtering_region->status = filtering_region_aligned;
    filtering_region_align_debug(filtering_candidates,
        filtering_region,pattern,true,matches,match_trace); // DEBUG
    return true; // OK
  }
}
