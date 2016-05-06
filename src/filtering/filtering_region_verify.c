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
    region_alignment_t* const region_alignment,
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
        region_alignment->distance_min_bound = ALIGN_DISTANCE_INF;
        return;
      }
    }
  }
  region_alignment->distance_min_bound = mismatches;
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
  const uint64_t max_error = filtering_region->max_error;
  const uint64_t text_base_offset = filtering_region->text_source_region_offset - filtering_region->key_source_region_offset;
  const uint8_t* const text = text_trace->text + text_base_offset;
  const uint64_t text_length = key_length;
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  // Check length
  if (text_length >= key_length) {
    // Verify Hamming
    filtering_region_verify_hamming_text(region_alignment,text,key,key_length,allowed_enc,max_error);
    if (region_alignment->distance_min_bound != ALIGN_DISTANCE_INF) {
      region_alignment->num_tiles = 1;
      region_alignment_tile_t* const alignment_tiles = mm_stack_calloc(mm_stack,1,region_alignment_tile_t,false);
      region_alignment->alignment_tiles = alignment_tiles;
      alignment_tiles->match_distance = region_alignment->distance_min_bound;
      alignment_tiles->text_begin_offset = text_base_offset;
      alignment_tiles->text_end_offset = text_base_offset + key_length;
    }
  } else {
    region_alignment->distance_min_bound = ALIGN_DISTANCE_INF;
  }
}
/*
 * Verify Leveshtein
 */
bool filtering_region_verify_levenshtein_kmer_filter(
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->text_length;
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  // Compile kmer-filter trimmed
  kmer_counting_t* kmer_counting;
  if (filtering_region->key_trimmed) {
    mm_stack_push_state(mm_stack);
    kmer_counting = mm_stack_alloc(mm_stack,kmer_counting_t);
    kmer_counting_compile(
        kmer_counting,pattern->key+filtering_region->key_trim_left,
        filtering_region->key_trimmed_length,filtering_region->max_error,mm_stack);
  } else {
    kmer_counting = &pattern->kmer_counting;
  }
  // Kmer Filter
  const uint64_t test_positive = kmer_counting_filter(kmer_counting,text,text_length);
  if (filtering_region->key_trimmed) mm_stack_pop_state(mm_stack);
  if (test_positive==ALIGN_DISTANCE_INF) {
    region_alignment->distance_min_bound = ALIGN_DISTANCE_INF;
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
    #ifdef FILTERING_REGION_VERIFY_CHECK_KMER_FILTER
    // DEBUG
    const bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
    uint64_t distance, match_column;
    bpm_compute_edit_distance(bpm_pattern,text,text_length,&distance,match_column,max_error,true);
    gem_cond_error_msg(distance != ALIGN_DISTANCE_INF,
        "Filtering.Region.Verify: K-mer filtering wrong discarding (edit-distance=%lu)",distance);
    #endif
    return false;
  }
  PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
  return true;
}
void filtering_region_verify_levenshtein_bpm(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    text_trace_t* const text_trace,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint8_t* const text = text_trace->text;
  const uint64_t num_pattern_tiles = bpm_pattern_tiles->num_pattern_tiles;
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  region_alignment_tile_t* const alignment_tiles = region_alignment->alignment_tiles;
  const uint64_t max_error = MIN(filtering_region->max_error,bpm_pattern->pattern_length);
  // Align tiles
  uint64_t max_remaining_error = filtering_region->max_error;
  uint64_t tile_pos, global_distance=0;
  PROF_ADD_COUNTER(GP_BMP_DISTANCE_NUM_TILES,num_pattern_tiles);
  for (tile_pos=0;tile_pos<num_pattern_tiles;++tile_pos) {
    region_alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
    if (alignment_tile->match_distance!=ALIGN_DISTANCE_INF) continue; // Next, already computed!
    if (global_distance <= max_error) {
      bpm_pattern_t* const bpm_pattern_tile = bpm_pattern_tiles+tile_pos;
      const uint64_t max_tile_error = MIN(max_remaining_error,bpm_pattern->pattern_length);
      const uint64_t tile_offset = alignment_tile->text_begin_offset;
      const uint64_t tile_wide = alignment_tile->text_end_offset-alignment_tile->text_begin_offset;
      uint64_t tile_distance, tile_match_column;
      bpm_compute_edit_distance(bpm_pattern_tile,text+tile_offset,
          tile_wide,&tile_distance,&tile_match_column,max_tile_error,false);
      PROF_INC_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED);
      // Store tile alignment
      if (tile_distance!=ALIGN_DISTANCE_INF) {
        const uint64_t tile_end_offset = tile_match_column+1;
        const uint64_t tile_tall = bpm_pattern_tile->pattern_length;
        const uint64_t tile_begin_offset = BOUNDED_SUBTRACTION(tile_end_offset,tile_tall+tile_distance,0);
        alignment_tile->match_distance = tile_distance;
        alignment_tile->text_end_offset = tile_offset + tile_end_offset;
        alignment_tile->text_begin_offset = tile_offset + tile_begin_offset;
        // Update distance
        max_remaining_error = BOUNDED_SUBTRACTION(max_remaining_error,tile_distance,0);
        global_distance += tile_distance;
        continue;
      }
    }
    // Store tile alignment
    global_distance = ALIGN_DISTANCE_INF;
  }
  // Setup alignment result
  region_alignment->distance_min_bound = global_distance;
}
void filtering_region_verify_levenshtein(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool kmer_filter,
    mm_stack_t* const mm_stack) {
  // Generalized Kmer-Counting filter [Prefilter]
  if (kmer_filter) {
    if (!filtering_region_verify_levenshtein_kmer_filter(
        filtering_region,pattern,text_trace,mm_stack)) return;
  }
  // Get BPM-Pattern
  bpm_pattern_t* bpm_pattern, *bpm_pattern_tiles;
  filtering_region_bpm_pattern_select(filtering_region,
      pattern,&bpm_pattern,&bpm_pattern_tiles,mm_stack);
  // Prepare Alignment
  filtering_region_alignment_prepare(filtering_region,
      bpm_pattern,bpm_pattern_tiles,mm_stack);
  // Myers's BPM algorithm [EditFilter]
  filtering_region_verify_levenshtein_bpm(filtering_candidates,
      filtering_region,bpm_pattern,bpm_pattern_tiles,text_trace,mm_stack);
}
/*
 * Verify (Switch)
 */
bool filtering_region_verify(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const bool kmer_filter) {
  PROF_START(GP_FC_VERIFY_CANDIDATES_REGION);
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  text_collection_t* const text_collection = filtering_candidates->text_collection;
  mm_stack_t* const mm_stack = filtering_candidates->mm_stack;
  region_alignment_t* const region_alignment = &filtering_region->region_alignment;
  // Check align-distance (already known or verified)
  if (region_alignment->distance_min_bound == ALIGN_DISTANCE_INF) {
    // Retrieve text-candidate
    filtering_region_retrieve_text(filtering_region,pattern,archive_text,text_collection,mm_stack);
    text_trace_t* const text_trace = text_collection_get_trace(
        text_collection,filtering_region->text_trace_offset);
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
        filtering_region_verify_levenshtein(
            filtering_candidates,filtering_region,
            pattern,text_trace,kmer_filter,mm_stack);
        break;
      case alignment_model_none:
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Check distance
  if (region_alignment->distance_min_bound != ALIGN_DISTANCE_INF) {
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
  const uint64_t max_filtering_error = filtering_region->max_error;
  // Select alignment model
  uint64_t num_matches_found;
  switch (filtering_candidates->search_parameters->alignment_model) {
    case alignment_model_levenshtein:
    case alignment_model_gap_affine:
      // 3. Myers's BPM algorithm
      num_matches_found = bpm_compute_edit_distance_all(
          bpm_pattern,filtering_candidates->filtering_regions,
          filtering_region->text_trace_offset,filtering_region->text_begin_position,
          text,text_length,max_filtering_error,filtering_candidates->mm_stack);
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
      text_trace_offset,index_position,text,text_length,
      max_filtering_error,filtering_candidates->mm_stack);
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
