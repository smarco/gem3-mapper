/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_verify.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_region_verify.h"
#include "align/align.h"
#include "align/align_bpm_distance.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION  GEM_DEEP_DEBUG
//#define FILTERING_REGION_VERIFY_CHECK_KMER_FILTER

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Verify Hamming
 */
void filtering_region_verify_hamming_text(
    filtering_region_t* const filtering_region,
    const uint8_t* const text,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool* const allowed_enc,
    const uint64_t max_mismatches) {
  // Check candidate
  uint64_t i, mismatches;
  for (i=0,mismatches=0;i<key_length;++i) {
    const uint8_t candidate_enc = text[i];
    // Check Mismatch
    if (!allowed_enc[candidate_enc] || candidate_enc != key[i]) {
      // Check Real Mismatch
      if (++mismatches > max_mismatches) {
        filtering_region->align_distance = ALIGN_DISTANCE_INF;
        filtering_region->align_match_begin_column = ALIGN_COLUMN_INF;
        filtering_region->align_match_end_column = ALIGN_COLUMN_INF;
        return;
      }
    }
  }
  filtering_region->align_distance = mismatches;
  filtering_region->align_match_begin_column = 0;
  filtering_region->align_match_end_column = key_length-1;
}
void filtering_region_verify_hamming(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool* const allowed_enc,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint8_t* const text = text_trace->text + filtering_region->base_begin_position_offset;
  const uint64_t text_length =
      filtering_region->base_end_position_offset -
      filtering_region->base_begin_position_offset;
  // Check length
  if (text_length >= key_length) {
    // Verify Hamming
    filtering_region_verify_hamming_text(filtering_region,text,key,key_length,allowed_enc,max_error);
    if (filtering_region->align_distance != ALIGN_DISTANCE_INF) {
      filtering_region->align_match_begin_column = filtering_region->base_begin_position_offset;
      filtering_region->align_match_end_column = filtering_region->base_end_position_offset;
    }
  } else {
    filtering_region->align_distance = ALIGN_DISTANCE_INF;
  }
}
/*
 * Verify Leveshtein
 */
void filtering_region_verify_levenshtein(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = filtering_region->end_position - filtering_region->begin_position;
  /*
   * Generalized Kmer-Counting filter [Prefilter]
   */
  if (!filtering_region->key_trimmed) {
    const uint64_t test_positive = kmer_counting_filter(&pattern->kmer_counting,text,text_length);
    if (test_positive==ALIGN_DISTANCE_INF) {
      filtering_region->align_distance = ALIGN_DISTANCE_INF;
      PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
#ifdef FILTERING_REGION_VERIFY_CHECK_KMER_FILTER
      // DEBUG
      const bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
      bpm_compute_edit_distance(bpm_pattern,
          text,eff_text_length,&filtering_region->align_distance,
          &filtering_region->align_match_end_column,max_error,true);
      gem_cond_error_msg(filtering_region->align_distance != ALIGN_DISTANCE_INF,
          "Filtering.Region.Verify: K-mer filtering wrong discarding (edit-distance=%lu)",
          filtering_region->align_distance);
#endif
      return;
    }
  } else {
    // Compile BPM-Pattern trimmed
    if (filtering_region->bpm_pattern_trimmed==NULL) {
      pattern_trimmed_init(pattern,&filtering_region->bpm_pattern_trimmed,
          &filtering_region->bpm_pattern_trimmed_tiles,filtering_region->key_trim_left,
          filtering_region->key_trim_right,mm_stack);
    }
  }
  PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
  /*
   * Myers's BPM algorithm [EditFilter]
   */
  bpm_pattern_t* const bpm_pattern = (!filtering_region->key_trimmed) ?
      pattern->bpm_pattern : filtering_region->bpm_pattern_trimmed;
  bpm_compute_edit_distance(bpm_pattern,
      text,text_length,&filtering_region->align_distance,
      &filtering_region->align_match_end_column,max_error,true);
}
/*
 * Verify (Switch)
 */
bool filtering_region_verify(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  PROF_START(GP_FC_VERIFY_CANDIDATES_REGION);
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  // Check align-distance (already known or verified)
  if (filtering_region->align_distance!=ALIGN_DISTANCE_INF) {
    filtering_region->status = filtering_region_accepted;
    filtering_region->align_match_end_column = filtering_region->base_end_position_offset-1;
    PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
    PROF_STOP(GP_FC_VERIFY_CANDIDATES_REGION);
    return true;
  }
  // Retrieve text-candidate
  filtering_region_retrieve_text(filtering_region,archive_text,text_collection,mm_stack);
  text_trace_t* const text_trace =
      text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  // Select alignment model
  switch (search_parameters->alignment_model) {
    case alignment_model_hamming: {
      // Verify Hamming
      filtering_region_verify_hamming(filtering_region,
          pattern,text_trace,search_parameters->allowed_enc,mm_stack);
      break;
    }
    case alignment_model_levenshtein:
    case alignment_model_gap_affine:
      // Verify Levenshtein
      filtering_region_verify_levenshtein(filtering_region,pattern,text_trace,mm_stack);
      break;
    case alignment_model_none:
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Check distance
  if (filtering_region->align_distance != ALIGN_DISTANCE_INF) {
    filtering_region->status = filtering_region_accepted;
    PROF_STOP(GP_FC_VERIFY_CANDIDATES_REGION);
    PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
    return true;
  } else {
    filtering_region->status = filtering_region_verified_discarded;
    PROF_STOP(GP_FC_VERIFY_CANDIDATES_REGION);
    PROF_INC_COUNTER(GP_DISCARDED_REGIONS);
    return false;
  }
  // Return fail
  PROF_STOP(GP_FC_VERIFY_CANDIDATES_REGION);
  return false;
}
uint64_t filtering_region_verify_multiple_hits(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  // Text (candidate)
  text_trace_t* const text_trace = text_collection_get_trace(
      filtering_candidates->text_collection,filtering_region->text_trace_offset);
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->text_length;
  // Pattern
  const bpm_pattern_t* const bpm_pattern = pattern->bpm_pattern;
  const uint64_t max_filtering_error = pattern->max_effective_filtering_error;
  // Select alignment model
  uint64_t num_matches_found;
  switch (filtering_candidates->search_parameters->alignment_model) {
    case alignment_model_levenshtein:
    case alignment_model_gap_affine:
      // 3. Myers's BPM algorithm
      num_matches_found = bpm_compute_edit_distance_all(bpm_pattern,
          filtering_candidates->filtering_regions,filtering_region->text_trace_offset,
          filtering_region->begin_position,text,text_length,max_filtering_error);
      break;
    case alignment_model_none:
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Return number of filtering regions added (accepted)
  PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS,num_matches_found);
  return num_matches_found;
}
uint64_t filtering_region_verify_extension(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t text_trace_offset,
    const uint64_t index_position,
    const pattern_t* const pattern) {
  PROFILE_START(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Text (candidate)
  text_trace_t* const text_trace = text_collection_get_trace(
      filtering_candidates->text_collection,text_trace_offset);
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->text_length;
  // Pattern
  const uint64_t max_filtering_error = pattern->max_effective_filtering_error;
  // 1. Hamming switch
  // TODO if (alignment_model==alignment_model_hamming) { }
  // 2. Generalized Counting filter
  // TODO
  // 3. Myers's BPM algorithm
  const bpm_pattern_t* const bpm_pattern = pattern->bpm_pattern;
  const uint64_t num_matches_found = bpm_compute_edit_distance_all(
      bpm_pattern,filtering_candidates->filtering_regions,
      text_trace_offset,index_position,text,text_length,max_filtering_error);
  PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS,num_matches_found);
  PROF_ADD_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_LENGTH,text_length);
  // Filter out already verified regions
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // Add to verified regions
  verified_region_t* verified_region;
  vector_alloc_new(filtering_candidates->verified_regions,verified_region_t,verified_region);
  verified_region->begin_position = index_position;
  verified_region->end_position = index_position + text_length;
  // Return number of filtering regions added (accepted)
  PROF_ADD_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_FOUND,num_matches_found);
  PROFILE_STOP(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  return num_matches_found;
}
