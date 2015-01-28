/*
 * PROJECT: GEMMapper
 * FILE: filtering_region.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region.h"

/*
 * Accessors
 */
GEM_INLINE text_trace_t* filtering_region_get_text_trace(
    const filtering_region_t* const filtering_region,const text_collection_t* const candidates_collection) {
  return text_collection_get_trace(candidates_collection,filtering_region->text_trace_offset);
}
/*
 * Sorting
 */
int region_matching_cmp_text_position(const region_matching_t* const a,const region_matching_t* const b) {
  return a->text_begin - b->text_begin;
}
GEM_INLINE void filtering_region_sort_regions_matching(const filtering_region_t* const filtering_region) {
  // Sort regions matching (region_matching_t) wrt their starting position in the text
  void* array = filtering_region->regions_matching;
  const size_t count = filtering_region->num_regions_matching;
  qsort(array,count,sizeof(region_matching_t),(int (*)(const void *,const void *))region_matching_cmp_text_position);
}
/*
 * Matching regions
 */
GEM_INLINE void filtering_region_exact_extend(
    filtering_region_t* const filtering_region,
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const bool* const allowed_enc) {
  // Extend all matching regions
  const uint64_t num_regions_matching = filtering_region->num_regions_matching;
  const uint64_t text_length = filtering_region->effective_end_position - filtering_region->effective_begin_position;
  const uint64_t last_region = num_regions_matching-1;
  uint64_t i, inc_coverage = 0;
  for (i=0;i<num_regions_matching;++i) {
    // Try to left extend
    region_matching_t* const region_matching = filtering_region->regions_matching + i;
    const int64_t left_key_max = (i==0) ? 0 : filtering_region->regions_matching[i-1].key_end;
    const int64_t left_text_max = (i==0) ? 0 : filtering_region->regions_matching[i-1].text_end;
    int64_t left_key = region_matching->key_begin-1;
    int64_t left_text = region_matching->text_begin-1;
    while (left_key_max<=left_key && left_text_max<=left_text) {
      // Check match
      const uint8_t candidate_enc = text[left_text];
      if (!allowed_enc[candidate_enc] || candidate_enc != key[left_key]) break;
      --left_key;
      --left_text;
      ++inc_coverage;
    }
    region_matching->key_begin = left_key+1;
    region_matching->text_begin = left_text+1;
    // Try to right extend
    const int64_t right_key_max = (i==last_region) ? key_length-1 : filtering_region->regions_matching[i+1].key_begin-1;
    const int64_t right_text_max = (i==last_region) ? text_length-1 : filtering_region->regions_matching[i+1].text_begin-1;
    int64_t right_key = region_matching->key_end;
    int64_t right_text = region_matching->text_end;
    while (right_key_max>=right_key && right_text_max>=right_text) {
      // Check match
      const uint8_t candidate_enc = text[right_text];
      if (!allowed_enc[candidate_enc] || candidate_enc != key[right_key]) break;
      ++right_key;
      ++right_text;
      ++inc_coverage;
    }
    region_matching->key_end = right_key;
    region_matching->text_end = right_text;
  }
  filtering_region->coverage += inc_coverage;
}
GEM_INLINE void filtering_region_chain_matching_regions(
    filtering_region_t* const filtering_region,
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const bool* const allowed_enc,
    const uint64_t max_error,mm_stack_t* const stack) {
  const uint64_t num_regions_matching = filtering_region->num_regions_matching;
  // Sort matching regions
  filtering_region_sort_regions_matching(filtering_region);
  // Check overlapping
  bool overlapping_text = false, unsorted_read = false, max_diff = false;
  region_matching_t* last_region_matching = NULL;
  uint64_t i, coverage = 0;
  for (i=0;i<num_regions_matching;++i) {
    region_matching_t* const region_matching = filtering_region->regions_matching + i;
    if (last_region_matching!=NULL) {
      if (last_region_matching->text_end > region_matching->text_begin) {
        overlapping_text = true; break; // TODO Can be fixed
      }
      if (last_region_matching->key_end > region_matching->key_begin) {
        unsorted_read = true; break;
      }
      const int64_t text_gap = region_matching->text_begin - last_region_matching->text_end;
      const int64_t key_gap = region_matching->key_begin - last_region_matching->key_end;
      const int64_t diff = text_gap-key_gap;
      if (ABS(diff) > max_error) {
        max_diff = true; break;
      }
    }
    coverage += region_matching->key_end - region_matching->key_begin;
    last_region_matching = region_matching;
  }
  if (overlapping_text || unsorted_read || max_diff) {
//    if (overlapping_text) fprintf(stderr,"> Overlapping-matches\n");
//    if (unsorted_read) fprintf(stderr,"> Unsorted-regions\n");
//    filtering_region_print_matching_regions(stderr,candidate_region,0);
    filtering_region->num_regions_matching = 0; // Disable region chaining
  } else {
    PROF_INC_COUNTER(GP_ACCEPTED_REGIONS_CHAINED);
    // Extend matching-regions
    filtering_region->coverage = coverage;
    PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS_COVERAGE,(100*filtering_region->coverage)/key_length);
    if (coverage < key_length) {
      filtering_region_exact_extend(filtering_region,key,key_length,text,allowed_enc);
    }
    PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS_EXT_COVERAGE,(100*filtering_region->coverage)/key_length);
  }
}
/*
 * (Re)Align
 */
GEM_INLINE void filtering_region_align(
    filtering_region_t* const filtering_region,const text_collection_t* const candidates_collection,
    const alignment_model_t alignment_model,const bool* const allowed_enc,
    const swg_penalties_t* const swg_penalties,const strand_t search_strand,
    const pattern_t* const pattern,const uint8_t* const key,const uint64_t key_length,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
  // Select Model
  if (filtering_region->align_distance==0 || alignment_model==alignment_model_none) {
    // Check column
    if (filtering_region->align_match_column == ALIGN_COLUMN_INF) {
      // Recalculate column
      region_matching_t* const region_matching = filtering_region->regions_matching;
      const uint64_t match_offset = region_matching->text_begin - region_matching->key_begin;
      filtering_region->align_match_column = key_length + match_offset - 1;
    }
    // Add exact match
    match_trace_t match_trace;
    matches_align_exact(matches,&match_trace,search_strand,swg_penalties,key_length,
        filtering_region->text_trace_offset,filtering_region->begin_position,
        filtering_region->align_distance,filtering_region->align_match_column+1);
    matches_add_match_trace_t(matches,&match_trace,true);
  } else {
    // Candidate
    const text_trace_t* const text_trace = filtering_region_get_text_trace(filtering_region,candidates_collection);
    uint8_t* const text = text_trace->text;
    // Select alignment model
    switch (alignment_model) {
      case alignment_model_hamming: {
        // Add hamming match
        match_trace_t match_trace;
        matches_align_hamming(matches,&match_trace,search_strand,allowed_enc,key,key_length,
            filtering_region->text_trace_offset,filtering_region->begin_position,
            text+(filtering_region->begin_position-filtering_region->effective_begin_position));
        matches_add_match_trace_t(matches,&match_trace,true);
        break;
      }
      case alignment_model_levenshtein: {
        // Add levenshtein match
        match_trace_t match_trace;
        matches_align_levenshtein(matches,&match_trace,search_strand,key,&pattern->bpm_pattern,
            filtering_region->text_trace_offset,filtering_region->begin_position,
            filtering_region->align_distance,text,filtering_region->align_match_column+1,
            filtering_region->regions_matching,filtering_region->num_regions_matching,mm_stack);
        if (match_trace.distance!=ALIGN_DISTANCE_INF) { // Double check
          matches_add_match_trace_t(matches,&match_trace,true);
        } else {
          filtering_region->align_distance = ALIGN_DISTANCE_INF;
        }
        break;
      }
      case alignment_model_gap_affine: {
        // Check matching-regions arrangement & find a region-chain
        //filtering_region_print_matching_regions(stderr,accepted_region,0);
        filtering_region_chain_matching_regions(filtering_region,key,key_length,
            text,allowed_enc,pattern->max_effective_filtering_error,mm_stack);
        //filtering_region_print_matching_regions(stderr,accepted_region,0);
        // Add Smith-Waterman-Gotoh match (affine-gap)
        const uint64_t text_length = filtering_region->effective_end_position-filtering_region->begin_position;
        match_trace_t match_trace;
        matches_align_smith_waterman_gotoh(matches,&match_trace,search_strand,allowed_enc,
            swg_penalties,key,key_length,filtering_region->text_trace_offset,
            filtering_region->begin_position,pattern->max_effective_filtering_error,
            text,text_length,filtering_region->regions_matching,filtering_region->num_regions_matching,mm_stack);
        //const uint64_t offset = match_trace.position-accepted_region->begin_position;
        //output_map_alignment_pretty(stderr,&match_trace,matches,
        //    (uint8_t* const)key,key_length,(uint8_t* const)text+offset,match_trace.effective_length,mm_stack);
        if (match_trace.distance!=ALIGN_DISTANCE_INF) { // Double check
          matches_add_match_trace_t(matches,&match_trace,true);
        } else {
          filtering_region->align_distance = ALIGN_DISTANCE_INF;
        }
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Set status aligned
  filtering_region->status = filtering_region_aligned;
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
        filtering_region->align_match_column = ALIGN_COLUMN_INF;
        return;
      }
    }
  }
  filtering_region->align_distance = mismatches;
  filtering_region->align_match_column = key_length-1;
}
GEM_INLINE int64_t filtering_region_verify_levenshtein(
    filtering_region_t* const candidate_region,
    const text_collection_t* const candidates_collection,
    const uint8_t* const key,const uint64_t key_length) {
  // Candidate
  const text_trace_t* const text_trace = filtering_region_get_text_trace(candidate_region,candidates_collection);
  const uint8_t* const text = text_trace->text;
  // Check
  gem_slog("[Checking position %lu]\n",candidate_region->begin_position);
  gem_slog("\tRange [%lu,%lu]\n",
      candidate_region->begin_position,candidate_region->effective_end_position);
  gem_slog("\tEffective Range [%lu,%lu]\n",
      candidate_region->effective_begin_position,candidate_region->effective_end_position);
  const uint64_t eff_text_length =
      candidate_region->effective_end_position - candidate_region->effective_begin_position;
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
    filtering_region_t* const candidate_region,const text_collection_t* const candidates_collection,
    const alignment_model_t alignment_model,const bool* const allowed_enc,const pattern_t* const pattern,
    const uint8_t* const key,const uint64_t key_length,const uint64_t max_filtering_error) {
  // DEBUG
  gem_cond_debug_block(DEBUG_ALIGN_CANDIDATES) {
    filtering_region_verify_levenshtein(candidate_region,candidates_collection,key,key_length);
  }
  /*
   * 1. Hamming switch
   */
  const text_trace_t* const text_trace = filtering_region_get_text_trace(candidate_region,candidates_collection);
  const uint8_t* const text = text_trace->text; // Candidate
  if (alignment_model==alignment_model_hamming) {
    const uint64_t text_length = candidate_region->effective_end_position - candidate_region->begin_position;
    if (text_length >= key_length) {
      const uint64_t text_offset = candidate_region->begin_position - candidate_region->effective_begin_position;
      filtering_region_verify_hamming(candidate_region,text+text_offset,
          key,key_length,allowed_enc,max_filtering_error);
      if (candidate_region->align_distance != ALIGN_DISTANCE_INF) {
        candidate_region->status = filtering_region_accepted;
        return true;
      }
    }
    candidate_region->status = filtering_region_discarded;
    return false;
  }
  /*
   * 2. Generalized Counting filter
   */
  const uint64_t eff_text_length = candidate_region->effective_end_position - candidate_region->effective_begin_position;
  // FIXME check that it might get to add 1 extra char
//    PROF_START(GP_FC_KMER_COUNTER_FILTER);
//    const uint64_t test_positive =
//        kmer_counting_filter(&pattern->kmer_counting,text,eff_text_length,max_effective_filtering_error);
//    PROF_STOP(GP_FC_KMER_COUNTER_FILTER);
  /*
   * 3. Myers's BPM algorithm
   */
  bpm_bound_distance_tiled((bpm_pattern_t* const)&pattern->bpm_pattern,
      text,eff_text_length,&candidate_region->align_distance,
      &candidate_region->align_match_column,max_filtering_error);
  if (candidate_region->align_distance != ALIGN_DISTANCE_INF) {
    PROF_INC_COUNTER(GP_LEVENSHTEIN_ACCEPTED);
    if (candidate_region->align_distance > max_filtering_error) {
      candidate_region->align_distance = max_filtering_error; // Adjust bound overestimations
    }
    candidate_region->begin_position = candidate_region->effective_begin_position; // Adjust position
    candidate_region->status = filtering_region_accepted;
    return true;
  }
  /*
   * 4. Local match (TODO Store as chunk matching ... try to extend borders)
   */
  candidate_region->status = filtering_region_discarded;
  return false;
}
/*
 * Display
 */
GEM_INLINE void filtering_region_print_matching_regions(
    FILE* const stream,filtering_region_t* const filtering_region,const uint64_t filtering_region_idx) {
  const uint64_t begin_position = filtering_region->effective_begin_position;
  fprintf(stream,"    #%lu -> [%lu,+%lu) \n",
      filtering_region_idx,begin_position,filtering_region->effective_end_position-begin_position);
  uint64_t j;
  for (j=0;j<filtering_region->num_regions_matching;++j) {
    region_matching_t* const region_matching = filtering_region->regions_matching + j;
    fprintf(stream,"      %lu.%lu) -> [%lu,%lu) ~> [+%lu,+%lu) \n",
        filtering_region_idx,j,region_matching->key_begin,region_matching->key_end,
        region_matching->text_begin,region_matching->text_end);
  }
}
