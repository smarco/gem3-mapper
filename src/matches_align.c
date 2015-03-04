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

/*
 * Exact
 */
GEM_INLINE void matches_align_exact(
    matches_t* const matches,match_trace_t* const match_trace,
    const strand_t strand,const swg_penalties_t* const swg_penalties,
    const uint64_t key_length,const uint64_t text_trace_offset,const uint64_t match_position,
    const uint64_t text_begin_offset,const uint64_t text_end_offset,const uint64_t match_distance) {
  PROF_START(GP_MATCHES_ALIGN_EXACT);
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->index_position = (match_distance==0) ? // Adjust position (if needed)
      match_position + (text_end_offset - key_length) : match_position;
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
  cigar_element.match_length = key_length;
  vector_insert(matches->cigar_buffer,cigar_element,cigar_element_t);
  PROF_STOP(GP_MATCHES_ALIGN_EXACT);
}
/*
 * Hamming
 */
GEM_INLINE void matches_align_hamming(
    matches_t* const matches,match_trace_t* const match_trace,
    const strand_t strand,const bool* const allowed_enc,
    const uint8_t* const key,const uint64_t key_length,
    const uint64_t text_trace_offset,const uint64_t match_position,
    const uint8_t* const text,const uint64_t text_begin_offset) {
  PROF_START(GP_MATCHES_ALIGN_HAMMING);
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->index_position = match_position;
  match_trace->strand = strand;
  // Hamming Check
  const uint64_t text_end_offset = text_begin_offset + key_length;
  match_trace->cigar_buffer_offset = vector_get_used(matches->cigar_buffer);
  vector_reserve_additional(matches->cigar_buffer,key_length);
  cigar_element_t* cigar_element = vector_get_free_elm(matches->cigar_buffer,cigar_element_t);
  uint64_t i, mismatches;
  for (i=text_begin_offset,mismatches=0;i<text_end_offset;++i) {
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
/*
 * Levenshtein
 */
GEM_INLINE void matches_align_levenshtein(
    matches_t* const matches,match_trace_t* const match_trace,const strand_t strand,
    const uint8_t* const key,const bpm_pattern_t* const bpm_pattern,
    const uint64_t text_trace_offset,const uint64_t match_position,
    uint8_t* const text,const uint64_t text_begin_offset,const uint64_t text_end_offset,
    const uint64_t max_distance,const region_matching_t* const regions_matching,
    const uint64_t num_regions_matching,mm_stack_t* const mm_stack) {
  PROF_START(GP_MATCHES_ALIGN_LEVENSHTEIN);
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->index_position = match_position + text_begin_offset; // Adjust position
  match_trace->strand = strand;
  // Levenshtein Align
  bpm_align_match(key,bpm_pattern,
      &match_trace->index_position,text+text_begin_offset,text_end_offset-text_begin_offset,max_distance,
      matches->cigar_buffer,&match_trace->cigar_buffer_offset,&match_trace->cigar_length,
      &match_trace->distance,&match_trace->effective_length,mm_stack);
  PROF_STOP(GP_MATCHES_ALIGN_LEVENSHTEIN);
}
/*
 * Smith-Waterman-Gotoh
 */
GEM_INLINE void matches_align_add_cigar_element(
    const swg_penalties_t* swg_penalties,vector_t* const cigar_vector,uint64_t* const cigar_length,
    int64_t* const effective_length,int32_t* const score,cigar_element_t* const cigar_element) {
  switch (cigar_element->type) {
    case cigar_match: {
      const int32_t match_length = cigar_element->match_length;
      matches_cigar_vector_append_match(cigar_vector,cigar_length,match_length);
      *score += swg_score_match(swg_penalties,match_length);
      *effective_length += match_length;
      break;
    }
    case cigar_mismatch:
      matches_cigar_vector_append_mismatch(cigar_vector,cigar_length,cigar_element->mismatch);
      *score += swg_score_mismatch(swg_penalties);
      ++(*effective_length);
      break;
    case cigar_ins: {
      const int32_t indel_length = cigar_element->indel.indel_length;
      matches_cigar_vector_append_insertion(cigar_vector,cigar_length,indel_length,cigar_element->indel.indel_text);
      *score += swg_score_insertion(swg_penalties,indel_length);
      *effective_length += indel_length;
      break;
    }
    case cigar_del:{
      const int32_t indel_length = cigar_element->indel.indel_length;
      matches_cigar_vector_append_deletion(cigar_vector,cigar_length,indel_length);
      *score += swg_score_deletion(swg_penalties,indel_length);
      *effective_length -= indel_length;
      break;
    }
    case cigar_soft_trim:
    case cigar_null:
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE void matches_align_matching_region(
    const region_matching_t* const region_matching,const bool* const allowed_enc,
    const swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,uint8_t* const text,vector_t* const cigar_vector,
    uint64_t* const cigar_length,int64_t* const effective_length,int32_t* const score,
    mm_stack_t* const mm_stack) {
  const uint64_t key_matching_length = region_matching->key_end-region_matching->key_begin;
  const uint64_t text_matching_length = region_matching->text_end-region_matching->text_begin;
  if (region_matching->error==0) {
    matches_cigar_vector_append_match(cigar_vector,cigar_length,key_matching_length);
    *effective_length += key_matching_length;
    *score += swg_score_match(swg_penalties,key_matching_length);
  } else if (region_matching->cigar_length > 0) {
    cigar_element_t* const scaffolding_buffer = vector_get_elm(cigar_vector,region_matching->cigar_buffer_offset,cigar_element_t);
    uint64_t i;
    for (i=0;i<region_matching->cigar_length;++i) {
      matches_align_add_cigar_element(swg_penalties,cigar_vector,cigar_length,effective_length,score,scaffolding_buffer+i);
    }
  } else {
    uint64_t dummy;
    swg_align_match(key+region_matching->key_begin,key_matching_length,allowed_enc,
        swg_query_profile,swg_penalties,&dummy,text+region_matching->text_begin,text_matching_length,
        region_matching->error,false,false,cigar_vector,cigar_length,effective_length,score,mm_stack);
  }
}
GEM_INLINE void matches_align_smith_waterman_gotoh_chained(
    matches_t* const matches,match_trace_t* const match_trace,const bool* const allowed_enc,
    const swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,const uint64_t key_length,uint8_t* const text,const uint64_t text_length,
    const uint64_t max_bandwidth,const region_matching_t* const regions_matching,
    const uint64_t num_regions_matching,mm_stack_t* const mm_stack) {
  // Initialize
  uint64_t key_chunk_begin_offset;
  uint64_t key_chunk_length;
  uint64_t text_chunk_end_offset;
  uint64_t text_chunk_begin_offset;
  uint64_t text_chunk_length;
  // Chain to FIRST matching region
  const region_matching_t* const first_region_matching = regions_matching;
  key_chunk_length = first_region_matching->key_begin;
  text_chunk_length = BOUNDED_ADDITION(key_chunk_length,max_bandwidth,first_region_matching->text_begin);
  text_chunk_begin_offset = first_region_matching->text_begin - text_chunk_length;
  uint8_t* const effective_text = text + text_chunk_begin_offset;
  match_trace->index_position += text_chunk_begin_offset;
  swg_align_match(key,key_chunk_length,allowed_enc,swg_query_profile,swg_penalties,
      &match_trace->index_position,effective_text,text_chunk_length,max_bandwidth,true,false,
      matches->cigar_buffer,&match_trace->cigar_length,&match_trace->effective_length,
      &match_trace->swg_score,mm_stack);
  // Add matching region to CIGAR
  matches_align_matching_region(first_region_matching,allowed_enc,
      swg_query_profile,swg_penalties,key,text,matches->cigar_buffer,
      &match_trace->cigar_length,&match_trace->effective_length,&match_trace->swg_score,mm_stack);
  // Chain matching regions
  uint64_t i, dummy;
  for (i=1;i<num_regions_matching;++i) {
    // Regions
    const region_matching_t* const prev_region_matching = regions_matching + (i-1);
    const region_matching_t* const current_region_matching = regions_matching + i;
    // Align gap between regions matching
    const uint64_t key_chunk_begin_offset = prev_region_matching->key_end;
    const uint64_t key_chunk_length = current_region_matching->key_begin - key_chunk_begin_offset;
    const uint64_t text_chunk_begin_offset = prev_region_matching->text_end;
    const uint64_t text_chunk_length = current_region_matching->text_begin - text_chunk_begin_offset;
    swg_align_match(key+key_chunk_begin_offset,key_chunk_length,allowed_enc,swg_query_profile,swg_penalties,
        &dummy,text+text_chunk_begin_offset,text_chunk_length,max_bandwidth,false,false,
        matches->cigar_buffer,&match_trace->cigar_length,&match_trace->effective_length,
        &match_trace->swg_score,mm_stack);
    // Add matching region to CIGAR
    matches_align_matching_region(current_region_matching,allowed_enc,
        swg_query_profile,swg_penalties,key,text,matches->cigar_buffer,
        &match_trace->cigar_length,&match_trace->effective_length,&match_trace->swg_score,mm_stack);
  }
  // Chain from LAST matching region
  const region_matching_t* const last_region_matching = regions_matching + (num_regions_matching-1);
  key_chunk_begin_offset = last_region_matching->key_end;
  key_chunk_length = key_length - last_region_matching->key_end;
  text_chunk_end_offset = BOUNDED_ADDITION(last_region_matching->text_end,key_chunk_length+max_bandwidth,text_length);
  text_chunk_begin_offset = last_region_matching->text_end;
  text_chunk_length = text_chunk_end_offset-last_region_matching->text_end;
  swg_align_match(key+key_chunk_begin_offset,key_chunk_length,allowed_enc,swg_query_profile,swg_penalties,
      &dummy,text+text_chunk_begin_offset,text_chunk_length,max_bandwidth,false,true,matches->cigar_buffer,
      &match_trace->cigar_length,&match_trace->effective_length,&match_trace->swg_score,mm_stack);
}
GEM_INLINE void matches_align_smith_waterman_gotoh(
    matches_t* const matches,match_trace_t* const match_trace,const strand_t strand,const bool* const allowed_enc,
    const swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,const uint64_t key_length,const uint64_t text_trace_offset,
    const uint64_t match_position,uint8_t* const text,const uint64_t text_length,
    const uint64_t text_begin_offset,const uint64_t text_end_offset,
    const uint64_t max_bandwidth,const region_matching_t* const regions_matching,
    const uint64_t num_regions_matching,mm_stack_t* const mm_stack) {
  PROF_START(GP_MATCHES_ALIGN_SWG);
  // Configure match-trace
  match_trace->trace_offset = text_trace_offset;
  match_trace->index_position = match_position;
  match_trace->strand = strand;
  match_trace->cigar_buffer_offset = vector_get_used(matches->cigar_buffer);
  match_trace->cigar_length = 0;
  match_trace->effective_length = 0;
  match_trace->swg_score = 0;
  // Check the number of matching regions
  if (num_regions_matching > 0) {
    // Chain matching regions and align gaps (SWG)
    matches_align_smith_waterman_gotoh_chained(matches,match_trace,allowed_enc,
        swg_query_profile,swg_penalties,key,key_length,text,text_length,
        max_bandwidth,regions_matching,num_regions_matching,mm_stack);
  } else {
    // Full SWG
    uint8_t* const effective_text = text + text_begin_offset;
    const uint64_t effective_text_length = text_end_offset - text_begin_offset;
    match_trace->index_position += text_begin_offset;
    swg_align_match(key,key_length,allowed_enc,swg_query_profile,swg_penalties,&match_trace->index_position,
        effective_text,effective_text_length,max_bandwidth,true,true,matches->cigar_buffer,
        &match_trace->cigar_length,&match_trace->effective_length,&match_trace->swg_score,mm_stack);
  }
  // Update distance (SWG score)
  if (match_trace->swg_score >= 0) {
    match_trace->distance = matches_cigar_compute_edit_distance(
        matches,match_trace->cigar_buffer_offset,match_trace->cigar_length);
  } else {
    match_trace->swg_score = ALIGN_DISTANCE_INF;
    match_trace->distance = ALIGN_DISTANCE_INF;
  }
  PROF_STOP(GP_MATCHES_ALIGN_SWG);
}
