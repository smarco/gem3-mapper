/*
 * PROJECT: GEMMapper
 * FILE: filtering_region.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region.h"

#include "match_scaffold.h"
#include "match_align_dto.h"

/*
 * Sorting
 */
int filtering_region_locator_cmp_position(const filtering_region_locator_t* const a,const filtering_region_locator_t* const b) {
  return a->position - b->position;
}
GEM_INLINE void filtering_region_locator_sort_positions(vector_t* const filtering_region_locators) {
  void* array = vector_get_mem(filtering_region_locators,filtering_region_locator_t);
  const size_t count = vector_get_used(filtering_region_locators);
  qsort(array,count,sizeof(filtering_region_locator_t),(int (*)(const void *,const void *))filtering_region_locator_cmp_position);
}
/*
 * (Re)Align
 */
GEM_INLINE bool filtering_region_align(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const text_collection_t* const text_collection,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,const bool left_gap_alignment,
    pattern_t* const pattern,matches_t* const matches,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  bool* const allowed_enc = search_parameters->allowed_enc;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Select Model
  if (filtering_region->align_distance==0 || alignment_model==alignment_model_none) {
    // Add exact match
    PROF_INC_COUNTER(GP_ACCEPTED_EXACT);
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
    PROF_INC_COUNTER(GP_ACCEPTED_INEXACT);
    // Retrieve Candidate (fetch text-trace if needed)
    const uint64_t text_length = filtering_region->end_position-filtering_region->begin_position;
    if (filtering_region->text_trace_offset == UINT64_MAX) {
      filtering_region->text_trace_offset = archive_text_retrieve(archive_text,text_collection,
          filtering_region->begin_position,text_length,false,mm_stack);
    }
    const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
    uint8_t* const text = text_trace->text;
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
        const uint64_t match_region_length = match_end_column-match_end_column;
        PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS_LENGTH,match_region_length);
        match_align_input_t align_input = {
            .key = key,
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
            .left_gap_alignment = left_gap_alignment,
        };
        match_align_levenshtein(matches,match_trace,&align_input,&align_parameters,mm_stack);
        break;
      }
      case alignment_model_gap_affine: {
        // Adjust alignment boundaries (to allow optimization)
        const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
        const uint64_t match_effective_range = key_length + 2*filtering_region->align_distance + max_bandwidth;
        const uint64_t text_offset_begin = BOUNDED_SUBTRACTION(filtering_region->align_match_end_column,match_effective_range,0);
        const uint64_t text_offset_end = BOUNDED_ADDITION(filtering_region->align_match_end_column,max_bandwidth,text_length-1);
        PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS_LENGTH,text_offset_end-text_offset_begin);
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
            .max_error = pattern->max_effective_filtering_error,
            .max_bandwidth = max_bandwidth,
            .left_gap_alignment = false, // FIXME left_gap_alignment
            .local_min_identity = as_parameters->local_min_identity_nominal,
            .min_coverage = as_parameters->region_scaffolding_coverage_threshold_nominal,
            .min_matching_length = as_parameters->region_scaffolding_min_length_nominal,
            .min_context_length = as_parameters->region_scaffolding_min_context_length_nominal,
            .allowed_enc = allowed_enc,
            .swg_penalties = swg_penalties,
            .cigar_curation = search_parameters->allow_cigar_curation,
        };
        match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
        // Scaffold the alignment
        if (search_parameters->allow_region_chaining) {
          match_scaffold_alignment(matches,&align_input,&align_parameters,match_scaffold,mm_stack);
        }
        // Smith-Waterman-Gotoh Alignment (Gap-affine)
        match_align_smith_waterman_gotoh(matches,match_trace,&align_input,&align_parameters,match_scaffold,mm_stack);
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Check (re)alignment result
  if (match_trace->distance!=ALIGN_DISTANCE_INF && match_trace->swg_score >= 0) { // TODO Check
    filtering_region->status = filtering_region_aligned;
    return true; // OK
  } else {
    filtering_region->align_distance = ALIGN_DISTANCE_INF;
    filtering_region->status = filtering_region_aligned_subdominant;
    return false; // Discarded
  }
}
GEM_INLINE bool filtering_region_local_align(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const text_collection_t* const text_collection,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,const bool left_gap_alignment,
    pattern_t* const pattern,matches_t* const matches,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
  // Parameters
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  bool* const allowed_enc = search_parameters->allowed_enc;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  // Select alignment model
  switch (alignment_model) {
    case alignment_model_hamming:
    case alignment_model_levenshtein:
      // There is no such a thing as local-alignment in these strict models
      return false;
      break;
    case alignment_model_gap_affine: {
      // Retrieve Candidate (fetch text-trace if needed)
      const uint64_t text_length = filtering_region->end_position-filtering_region->begin_position;
      if (filtering_region->text_trace_offset == UINT64_MAX) {
        filtering_region->text_trace_offset = archive_text_retrieve(archive_text,text_collection,
            filtering_region->begin_position,text_length,false,mm_stack);
      }
      const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
      uint8_t* const text = text_trace->text;
      // Adjust alignment boundaries (to allow optimization)
      const uint64_t max_bandwidth = pattern->max_effective_bandwidth;
      const uint64_t text_offset_begin = 0;
      const uint64_t text_offset_end = text_length-1;
      PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS_LENGTH,text_offset_end-text_offset_begin);
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
          .max_error = BOUNDED_SUBTRACTION(key_length,as_parameters->local_min_identity_nominal,0),
          .max_bandwidth = max_bandwidth,
          .left_gap_alignment = left_gap_alignment,
          .local_min_identity = as_parameters->local_min_identity_nominal,
          .min_coverage = as_parameters->region_scaffolding_coverage_threshold_nominal,
          .min_matching_length = as_parameters->region_scaffolding_min_length_nominal,
          .min_context_length = as_parameters->region_scaffolding_min_context_length_nominal,
          .allowed_enc = allowed_enc,
          .swg_penalties = swg_penalties,
          .cigar_curation = search_parameters->allow_cigar_curation,
      };
      match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
      // Scaffold from Levenshtein-alignment
      PROF_INC_COUNTER(GP_MATCHING_REGIONS_SCAFFOLDED);
      PROF_START(GP_MATCHING_REGIONS_SCAFFOLD);
      const bool valid_scaffold = match_scaffold_levenshtein(matches,&align_input,&align_parameters,match_scaffold,mm_stack);
      PROF_STOP(GP_MATCHING_REGIONS_SCAFFOLD);
      PROF_ADD_COUNTER(GP_MATCHING_REGIONS_SCAFFOLD_COVERAGE,(100*match_scaffold->scaffolding_coverage)/key_length);
      if (!valid_scaffold) return false;
      // Smith-Waterman-Gotoh Alignment (Gap-affine)
      match_align_smith_waterman_gotoh(matches,match_trace,&align_input,&align_parameters,match_scaffold,mm_stack);
      return (match_trace->distance!=ALIGN_DISTANCE_INF); // (Re)alignment result
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}
/*
 * Verify
 */
GEM_INLINE void filtering_region_verify_hamming(
    filtering_region_t* const filtering_region,const uint8_t* const text,
    const uint8_t* const key,const uint64_t key_length,
    const bool* const allowed_enc,const uint64_t max_mismatches) {
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
GEM_INLINE int64_t filtering_region_verify_levenshtein(
    filtering_region_t* const candidate_region,
    const text_collection_t* const text_collection,
    const uint8_t* const key,const uint64_t key_length) {
  // USAGE DEBUG
  // gem_cond_debug_block(DEBUG_ALIGN_CANDIDATES) {
  //   filtering_region_verify_levenshtein(filtering_region,text_collection,key,key_length);
  // }
  // Candidate
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,candidate_region->text_trace_offset);
  const uint8_t* const text = text_trace->text;
  // Check
  const uint64_t base_position = candidate_region->begin_position+candidate_region->base_position_offset;
  gem_slog("[Checking position %lu]\n",base_position);
  gem_slog("\tRange [%lu,%lu]\n",base_position,candidate_region->end_position);
  gem_slog("\tEffective Range [%lu,%lu]\n",
      candidate_region->begin_position,candidate_region->end_position);
  const uint64_t eff_text_length =
      candidate_region->end_position - candidate_region->begin_position;
  uint64_t i, dp_position;
  const int64_t dp_distance = align_levenshtein_get_distance(
      (const char * const)key,key_length,(const char * const)text,eff_text_length,true,&dp_position);
  gem_slog("\tDP-Alignment (distance=%lu,position=%lu)\n",dp_distance,dp_position);
  gem_slog("\tPattern: ");
  for (i=0;i<key_length;++i) gem_slog("%c",dna_decode(key[i]));
  gem_slog("\n\tText: ");
  for (i=0;i<eff_text_length;++i) gem_slog("%c",dna_decode(text[i]));
  gem_slog("\n");
  // Return distance
  return dp_distance;
}
GEM_INLINE bool filtering_region_verify(
    filtering_region_t* const filtering_region,const text_collection_t* const text_collection,
    search_parameters_t* const search_parameters,const pattern_t* const pattern) {
  // Parameters
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint8_t* const text = text_trace->text; // Candidate
  const uint64_t max_filtering_error = pattern->max_effective_filtering_error;
  // Select alignment model
  switch (search_parameters->alignment_model) {
    case alignment_model_hamming: { // 1. Hamming switch
      const uint8_t* const key = pattern->key;
      const uint64_t key_length = pattern->key_length;
      const uint64_t text_length =
          filtering_region->end_position - filtering_region->begin_position - filtering_region->base_position_offset;
      const bool* const allowed_enc = search_parameters->allowed_enc;
      if (text_length >= key_length) {
        filtering_region_verify_hamming(filtering_region,
            text+filtering_region->base_position_offset,key,key_length,allowed_enc,max_filtering_error);
        if (filtering_region->align_distance != ALIGN_DISTANCE_INF) {
          // Adjust using text_offset
          filtering_region->align_match_begin_column += text_length;
          filtering_region->align_match_end_column += text_length;
          filtering_region->status = filtering_region_accepted;
          return true;
        }
      }
      filtering_region->status = filtering_region_discarded;
      return false;
      break;
    }
    case alignment_model_levenshtein:
    case alignment_model_gap_affine: {
      const uint64_t eff_text_length = filtering_region->end_position - filtering_region->begin_position;
      // 2. Generalized Counting filter
      // PROF_START(GP_FC_KMER_COUNTER_FILTER);
      // const uint64_t test_positive = kmer_counting_filter(&pattern->kmer_counting,text,eff_text_length,max_filtering_error);
      // PROF_STOP(GP_FC_KMER_COUNTER_FILTER);

//      // 3. Multiple Myers's BPM algorithm
//      const uint64_t num_matches_found = bpm_search_all(
//          (bpm_pattern_t* const)&pattern->bpm_pattern,filtering_regions,
//          text_trace_offset,index_position,text,text_length,max_filtering_error);

      // 3. Myers's BPM algorithm
      bpm_get_distance_cutoff_tiled(
          (bpm_pattern_t* const)&pattern->bpm_pattern,text,eff_text_length,&filtering_region->align_distance,
          &filtering_region->align_match_end_column,max_filtering_error);
      if (filtering_region->align_distance != ALIGN_DISTANCE_INF) {
        PROF_INC_COUNTER(GP_LEVENSHTEIN_ACCEPTED);
        filtering_region->status = filtering_region_accepted;
        return true;
      } else {
        filtering_region->status = filtering_region_discarded;
        return false;
      }
      break;
    }
    case alignment_model_none:
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Return fail
  return false;
}
GEM_INLINE uint64_t filtering_region_verify_multiple_hits(
    vector_t* const filtering_regions,filtering_region_t* const filtering_region,
    const text_collection_t* const text_collection,search_parameters_t* const search_parameters,
    const pattern_t* const pattern) {
  // Text (candidate)
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->length;
  // Pattern
  const bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
  const uint64_t max_filtering_error = pattern->max_effective_filtering_error;
  // Select alignment model
  uint64_t num_matches_found;
  switch (search_parameters->alignment_model) {
//    case alignment_model_hamming:
    case alignment_model_levenshtein:
    case alignment_model_gap_affine:
      // 3. Myers's BPM algorithm
      num_matches_found = bpm_search_all(bpm_pattern,filtering_regions,
          filtering_region->text_trace_offset,filtering_region->begin_position,
          text,text_length,max_filtering_error);
      break;
    case alignment_model_none:
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Return number of filtering regions added (accepted)
  return num_matches_found;
}
GEM_INLINE uint64_t filtering_region_verify_extension(
    vector_t* const filtering_regions,vector_t* const verified_regions,
    const text_collection_t* const text_collection,
    const uint64_t const text_trace_offset,const uint64_t index_position,
    search_parameters_t* const search_parameters,const pattern_t* const pattern) {
  // Text (candidate)
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->length;
  // Pattern
  const uint64_t max_filtering_error = pattern->max_effective_filtering_error;
  // 1. Hamming switch
  // TODO if (alignment_model==alignment_model_hamming) { }
  // 2. Generalized Counting filter
  // TODO
  // 3. Myers's BPM algorithm
  const uint64_t num_matches_found = bpm_search_all(
      (bpm_pattern_t* const)&pattern->bpm_pattern,filtering_regions,
      text_trace_offset,index_position,text,text_length,max_filtering_error);
  // Add to verified regions
  verified_region_t* verified_region;
  vector_alloc_new(verified_regions,verified_region_t,verified_region);
  verified_region->begin_position = index_position;
  verified_region->end_position = index_position + text_length;
  // Return number of filtering regions added (accepted)
  return num_matches_found;
}
