/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_verify.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region_verify.h"
#include "align.h"
#include "align_bpm_distance.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION  GEM_DEEP_DEBUG
//#define FILTERING_REGION_VERIFY_CHECK_KMER_FILTER

/*
 * Verify
 */
void filtering_region_verify_hamming(
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
int64_t filtering_region_verify_levenshtein(
    filtering_region_t* const candidate_region,
    const text_collection_t* const text_collection,
    const uint8_t* const key,const uint64_t key_length) {
  // Candidate
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,candidate_region->text_trace_offset);
  const uint8_t* const text = text_trace->text;
  // Check
  const uint64_t base_position = candidate_region->begin_position+candidate_region->base_position_offset;
  gem_slog("[Checking position %"PRIu64"]\n",base_position);
  gem_slog("\tRange [%"PRIu64",%"PRIu64"]\n",base_position,candidate_region->end_position);
  gem_slog("\tEffective Range [%"PRIu64",%"PRIu64"]\n",
      candidate_region->begin_position,candidate_region->end_position);
  const uint64_t eff_text_length =
      candidate_region->end_position - candidate_region->begin_position;
  uint64_t i, dp_position;
  const int64_t dp_distance = align_dp_compute_edit_distance(
      (const char * const)key,key_length,(const char * const)text,eff_text_length,true,&dp_position);
  gem_slog("\tDP-Alignment (distance=%"PRIu64",position=%"PRIu64")\n",dp_distance,dp_position);
  gem_slog("\tPattern: ");
  for (i=0;i<key_length;++i) gem_slog("%c",dna_decode(key[i]));
  gem_slog("\n\tText: ");
  for (i=0;i<eff_text_length;++i) gem_slog("%c",dna_decode(text[i]));
  gem_slog("\n");
  // Return distance
  return dp_distance;
}
bool filtering_region_verify(
    filtering_region_t* const filtering_region,const text_collection_t* const text_collection,
    search_parameters_t* const search_parameters,const pattern_t* const pattern) {
  PROF_START(GP_FC_VERIFY_CANDIDATE_REGION);
  // Parameters
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint8_t* const text = text_trace->text; // Candidate
  const uint64_t max_error = pattern->max_effective_filtering_error;
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
            text+filtering_region->base_position_offset,key,key_length,allowed_enc,max_error);
        if (filtering_region->align_distance != ALIGN_DISTANCE_INF) {
          // Adjust using text_offset
          filtering_region->align_match_begin_column += text_length;
          filtering_region->align_match_end_column += text_length;
          filtering_region->status = filtering_region_accepted;
          PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGION);
          return true;
        }
      }
      filtering_region->status = filtering_region_verified_discarded;
      PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGION);
      return false;
      break;
    }
    case alignment_model_levenshtein:
    case alignment_model_gap_affine: {
      const uint64_t eff_text_length = filtering_region->end_position - filtering_region->begin_position;
      // [Prefilter] Generalized Kmer-Counting filter
      const uint64_t test_positive = kmer_counting_filter(&pattern->kmer_counting,text,eff_text_length);
      if (test_positive==ALIGN_DISTANCE_INF) {
        filtering_region->status = filtering_region_verified_discarded;
        PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGION);
        PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
#ifdef FILTERING_REGION_VERIFY_CHECK_KMER_FILTER
        // DEBUG
        const bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
        bpm_compute_edit_distance(bpm_pattern,
            text,eff_text_length,&filtering_region->align_distance,
            &filtering_region->align_match_end_column,max_error,true);
        gem_cond_error_msg(filtering_region->align_distance != ALIGN_DISTANCE_INF,
            "Filtering.Region.Verify: K-mer filtering wrong discarding (edit-distance=%lu)",filtering_region->align_distance);
#endif
        return false;
      } else {
        PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
      }
      // [EditFilter] Myers's BPM algorithm
      const bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
      bpm_compute_edit_distance(bpm_pattern,
          text,eff_text_length,&filtering_region->align_distance,
          &filtering_region->align_match_end_column,max_error,true);
      if (filtering_region->align_distance != ALIGN_DISTANCE_INF) {
        filtering_region->status = filtering_region_accepted;
        PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGION);
        PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
        return true;
      } else {
        filtering_region->status = filtering_region_verified_discarded;
        PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGION);
        PROF_INC_COUNTER(GP_DISCARDED_REGIONS);
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
  PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGION);
  return false;
}
uint64_t filtering_region_verify_multiple_hits(
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
      num_matches_found = bpm_compute_edit_distance_all(bpm_pattern,filtering_regions,
          filtering_region->text_trace_offset,filtering_region->begin_position,
          text,text_length,max_filtering_error);
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
    vector_t* const filtering_regions,vector_t* const verified_regions,
    const text_collection_t* const text_collection,const uint64_t text_trace_offset,
    const uint64_t index_position,search_parameters_t* const search_parameters,
    const pattern_t* const pattern) {
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
  const uint64_t num_matches_found = bpm_compute_edit_distance_all(
      (bpm_pattern_t* const)&pattern->bpm_pattern,filtering_regions,
      text_trace_offset,index_position,text,text_length,max_filtering_error);
  PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS,num_matches_found);


  // Filter out already verified regions
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

  // Add to verified regions
  verified_region_t* verified_region;
  vector_alloc_new(verified_regions,verified_region_t,verified_region);
  verified_region->begin_position = index_position;
  verified_region->end_position = index_position + text_length;
  // Return number of filtering regions added (accepted)
  return num_matches_found;
}

