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
    const strand_t strand,const uint64_t key_length,
    const uint64_t text_trace_offset,const uint64_t match_position,
    const uint64_t match_distance,const uint64_t match_length) {
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->position = (match_distance==0) ? // Adjust position (if needed)
      match_position + (match_length - key_length) : match_position;
  match_trace->distance = match_distance;
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
}
GEM_INLINE void matches_align_hamming(
    matches_t* const matches,match_trace_t* const match_trace,
    const strand_t strand,const bool* const allowed_enc,const uint8_t* const key,const uint64_t key_length,
    const uint64_t text_trace_offset,const uint64_t match_position,const uint8_t* const text) {
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
}
GEM_INLINE void matches_align_levenshtein(
    matches_t* const matches,match_trace_t* const match_trace,
    const strand_t strand,const uint8_t* const key,const bpm_pattern_t* const bpm_pattern,
    const uint64_t text_trace_offset,const uint64_t match_position,const uint64_t max_distance,
    const uint8_t* const text,const uint64_t text_length,const region_matching_t* const regions_matching,
    const uint64_t num_regions_matching,mm_stack_t* const mm_stack) {
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->position = match_position;
  match_trace->strand = strand;
  // Levenshtein Align
  bpm_align_match(key,bpm_pattern,
      &match_trace->position,text,text_length,max_distance,
      matches->cigar_buffer,&match_trace->cigar_buffer_offset,&match_trace->cigar_length,
      &match_trace->distance,&match_trace->effective_length,mm_stack);
}
GEM_INLINE void matches_align_smith_waterman_gotoh(
    matches_t* const matches,match_trace_t* const match_trace,
    const strand_t strand,const uint8_t* const key,const uint64_t key_length,
    const uint64_t text_trace_offset,const uint64_t match_position,const uint8_t* const text,const uint64_t text_length,
    const region_matching_t* const regions_matching,const uint64_t num_regions_matching,mm_stack_t* const mm_stack) {
  GEM_NOT_IMPLEMENTED(); // TODO
}


// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
//const uint64_t match_length = matches_get_effective_length(
//    matches,seq_length,match_trace.cigar_buffer_offset,match_trace.cigar_length);
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute
// FIXME match_length should be a match attribute







