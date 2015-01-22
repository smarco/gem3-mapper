/*
 * PROJECT: GEMMapper
 * FILE: matches_align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches_align.h"

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// archive_select_curate_match(select_parameters,matches,match_trace);
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

GEM_INLINE void matches_align_exact(
    matches_t* const matches,match_trace_t* const match_trace,
    const strand_t strand,const swg_penalties_t* const swg_penalties,
    const uint64_t key_length,const uint64_t text_trace_offset,
    uint64_t match_position,const uint64_t match_distance,const uint64_t match_length) {
  PROF_START(GP_MATCHES_ALIGN_EXACT);
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->position = (match_distance==0) ? // Adjust position (if needed)
      match_position + (match_length - key_length) : match_position;
  match_trace->distance = match_distance;
  match_trace->swg_score = swg_penalties->generic_match_score*key_length;
  match_trace->strand = strand;
  // Exact match
  match_trace->cigar_buffer_offset = vector_get_used(matches->cigar_buffer);
  match_trace->cigar_length = 1;
  match_trace->effective_length = key_length;
  // Insert all-matching CIGAR
  cigar_element_t cigar_element;
  cigar_element.type = cigar_match;
  cigar_element.length = key_length;
  vector_insert(matches->cigar_buffer,cigar_element,cigar_element_t);
  PROF_STOP(GP_MATCHES_ALIGN_EXACT);
}
GEM_INLINE void matches_align_hamming(
    matches_t* const matches,match_trace_t* const match_trace,
    const strand_t strand,const bool* const allowed_enc,const uint8_t* const key,const uint64_t key_length,
    const uint64_t text_trace_offset,const uint64_t match_position,const uint8_t* const text) {
  PROF_START(GP_MATCHES_ALIGN_HAMMING);
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->position = match_position;
  match_trace->strand = strand;
  // Hamming Check
  match_trace->cigar_buffer_offset = vector_get_used(matches->cigar_buffer);
  vector_reserve_additional(matches->cigar_buffer,key_length);
  cigar_element_t* cigar_element = vector_get_free_elm(matches->cigar_buffer,cigar_element_t);
  uint64_t i, mismatches;
  for (i=0,mismatches=0;i<key_length;++i) {
    const uint8_t candidate_enc = text[i];
    // Check Mismatch
    if (!allowed_enc[candidate_enc] || candidate_enc != key[i]) {
      ++mismatches;
      cigar_element->type = cigar_mismatch;
      cigar_element->mismatch = candidate_enc;
    }
  }
  match_trace->distance = mismatches;
  match_trace->cigar_length = mismatches;
  vector_update_used(matches->cigar_buffer,cigar_element);
  match_trace->effective_length = key_length;
  PROF_STOP(GP_MATCHES_ALIGN_HAMMING);
}
GEM_INLINE void matches_align_levenshtein(
    matches_t* const matches,match_trace_t* const match_trace,
    const strand_t strand,const uint8_t* const key,const bpm_pattern_t* const bpm_pattern,
    const uint64_t text_trace_offset,const uint64_t match_position,const uint64_t max_distance,
    const uint8_t* const text,const uint64_t text_length,const region_matching_t* const regions_matching,
    const uint64_t num_regions_matching,mm_stack_t* const mm_stack) {
  PROF_START(GP_MATCHES_ALIGN_LEVENSHTEIN);
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->position = match_position;
  match_trace->strand = strand;
  // Levenshtein Align
  bpm_align_match(key,bpm_pattern,
      &match_trace->position,text,text_length,max_distance,
      matches->cigar_buffer,&match_trace->cigar_buffer_offset,&match_trace->cigar_length,
      &match_trace->distance,&match_trace->effective_length,mm_stack);
  PROF_STOP(GP_MATCHES_ALIGN_LEVENSHTEIN);
}
GEM_INLINE void matches_align_smith_waterman_gotoh(
    matches_t* const matches,match_trace_t* const match_trace,const strand_t strand,const bool* const allowed_enc,
    const swg_penalties_t* const swg_penalties,const uint8_t* const key,const uint64_t key_length,
    const uint64_t text_trace_offset,const uint64_t match_position,const uint64_t max_distance,
    const uint8_t* const text,const uint64_t text_length,const region_matching_t* const regions_matching,
    const uint64_t num_regions_matching,mm_stack_t* const mm_stack) {
  PROF_START(GP_MATCHES_ALIGN_SWG);
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->position = match_position;
  match_trace->strand = strand;
  match_trace->cigar_buffer_offset = vector_get_used(matches->cigar_buffer);
  match_trace->cigar_length = 0;
  match_trace->effective_length = 0;

//  swg_align_match_full_32b(key,key_length,swg_penalties,&match_trace->position,text,text_length,
//      matches->cigar_buffer,&match_trace->cigar_length,&match_trace->effective_length,&match_trace->score,mm_stack);
//  swg_align_match_banded_32b(key,key_length,swg_penalties,&match_trace->position,
//          text,text_length,max_distance+1,true,true,matches->cigar_buffer,
//          &match_trace->cigar_length,&match_trace->effective_length,&match_trace->score,mm_stack);
//  swg_align_match_banded_32b_opt(key,key_length,swg_penalties,&match_trace->position,
//          text,text_length,max_distance+1,true,true,matches->cigar_buffer,
//          &match_trace->cigar_length,&match_trace->effective_length,&match_trace->score,mm_stack);

  // SWG Align
  int32_t score = 0; // Init
  if (num_regions_matching > 0) {
    // Initialize
    const uint32_t generic_match_score = swg_penalties->generic_match_score; // TODO Weighted Matrices
    uint64_t dummy;
    /*
     * Chain to FIRST matching region
     */
    const region_matching_t* first_region_matching = regions_matching;
    const uint64_t key_chunk_length = first_region_matching->key_begin;
    const uint64_t text_chunk_length = first_region_matching->text_begin;
    swg_align_match(key,key_chunk_length,allowed_enc,swg_penalties,&match_trace->position,
        text,text_chunk_length,max_distance+1,true,false,matches->cigar_buffer,
        &match_trace->cigar_length,&match_trace->effective_length,&score,mm_stack);
    PROF_ADD_COUNTER(GP_MATCHES_ALIGN_SWG_CELLS,key_chunk_length*text_chunk_length);
    /*
     * Chain matching regions
     */
    // Add matching region to CIGAR
    const uint64_t region_matching_length = first_region_matching->key_end-first_region_matching->key_begin;
    matches_cigar_buffer_append_match(matches->cigar_buffer,&match_trace->cigar_length,region_matching_length);
    match_trace->effective_length += region_matching_length;
    score += generic_match_score*region_matching_length;
    uint64_t i;
    for (i=1;i<num_regions_matching;++i) {
      const region_matching_t* prev_region_matching = regions_matching + (i-1);
      const region_matching_t* current_region_matching = prev_region_matching + 1;
      // Align gap between regions matching
      const uint64_t key_gap_length = current_region_matching->key_begin - prev_region_matching->key_end;
      const uint64_t text_gap_length = current_region_matching->text_begin - prev_region_matching->text_end;
      const uint64_t key_chunk_offset = prev_region_matching->key_end;
      const uint64_t text_chunk_offset = prev_region_matching->text_end;
      swg_align_match(key+key_chunk_offset,key_gap_length,allowed_enc,swg_penalties,&dummy,
          text+text_chunk_offset,text_gap_length,max_distance+1,false,false,matches->cigar_buffer,
          &match_trace->cigar_length,&match_trace->effective_length,&score,mm_stack);
      PROF_ADD_COUNTER(GP_MATCHES_ALIGN_SWG_CELLS,key_gap_length*text_gap_length);
      // Add matching region to CIGAR
      const uint64_t region_matching_length = regions_matching[i].key_end-regions_matching[i].key_begin;
      matches_cigar_buffer_append_match(matches->cigar_buffer,&match_trace->cigar_length,region_matching_length);
      match_trace->effective_length += region_matching_length;
      score += generic_match_score*region_matching_length;
    }
    /*
     * Chain from LAST matching region
     */
    const region_matching_t* last_region_matching = regions_matching + (num_regions_matching-1);
    const uint64_t key_gap_length = key_length - last_region_matching->key_end;
    const uint64_t text_gap_length = text_length - last_region_matching->text_end;
    const uint64_t key_gap_offset = last_region_matching->key_end;
    const uint64_t text_gap_offset = last_region_matching->text_end;
    swg_align_match(key+key_gap_offset,key_gap_length,allowed_enc,swg_penalties,&dummy,
        text+text_gap_offset,text_gap_length,max_distance+1,false,true,matches->cigar_buffer,
        &match_trace->cigar_length,&match_trace->effective_length,&score,mm_stack);
    PROF_ADD_COUNTER(GP_MATCHES_ALIGN_SWG_CELLS,key_gap_length*text_gap_length);
  } else {
    swg_align_match(key,key_length,allowed_enc,swg_penalties,&match_trace->position,
        text,text_length,max_distance+1,true,true,matches->cigar_buffer,
        &match_trace->cigar_length,&match_trace->effective_length,&score,mm_stack);
    PROF_ADD_COUNTER(GP_MATCHES_ALIGN_SWG_CELLS,key_length*text_length);
  }
  // Update distance (SWG score)
  if (score >= 0) {
    match_trace->swg_score = score;
    match_trace->distance = matches_cigar_calculate_edit_distance(
        matches,match_trace->cigar_buffer_offset,match_trace->cigar_length); // Event distance
  } else {
    match_trace->swg_score = ALIGN_DISTANCE_INF;
    match_trace->distance = ALIGN_DISTANCE_INF;
  }
  PROF_STOP(GP_MATCHES_ALIGN_SWG);
  PROF_ADD_COUNTER(GP_MATCHES_ALIGN_SWG_CELLS_POTENTIAL,key_length*text_length);
}


