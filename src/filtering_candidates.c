/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "matches_align.h"
#include "archive_text_retrieve.h"

/*
 * Debug
 */
#define DEBUG_ALIGN_CANDIDATES  false
#define DEBUG_ALIGN_LEVENSHTEIN false
#define DEBUG_REGIONS_MATCHING  false

/*
 * Constants
 */
#define FC_INIT_REGIONS_BUFFER               100
#define FC_INIT_CANDIDATE_POSITIONS         1000

#define FC_DECODE_NUM_POSITIONS_PREFETCHED          10
#define FC_RETRIEVE_SAMPLE_NUM_POSITIONS_PREFETCHED 10

#define FC_DISTANCE_EXCEED    UINT64_MAX
#define FC_POSITION_DISCARDED UINT64_MAX

/*
 * Candidate Position
 */
// (1) Alias input fields
#define candidate_region_index_position          overloaded_field_0
// (2) Alias decoding positions
#define candidate_decode_distance                overloaded_field_0
#define candidate_decode_sampled_pos             overloaded_field_1
#define candidate_region_text_position           overloaded_field_3
// (3) Alias re-aligning candidates
#define candidate_effective_end_position         overloaded_field_0
#define candidate_effective_begin_position       overloaded_field_1
#define candidate_begin_position                 overloaded_field_2
#define candidate_text_trace_offset              overloaded_field_3
#define candidate_align_distance                 overloaded_field_4
#define candidate_align_match_column             overloaded_field_5
// Candidate Position
typedef struct {
  // Region Info
  region_t* candidate_region;
  locator_interval_t* locator_interval;
  // Internals (Overloaded fields)
  uint64_t overloaded_field_0;
  uint64_t overloaded_field_1;
  uint64_t overloaded_field_2;
  uint64_t overloaded_field_3;
} candidate_position_t;
/*
 * Candidate Region
 */
typedef struct {
  // Regions Matching
  uint64_t num_regions_matching;
  region_matching_t* regions_matching;
  uint64_t coverage;
  // Internals (Overloaded fields)
  uint64_t overloaded_field_0;
  uint64_t overloaded_field_1;
  uint64_t overloaded_field_2;
  uint64_t overloaded_field_3;
  uint64_t overloaded_field_4;
  uint64_t overloaded_field_5;
} candidate_region_t;

/*
 * Setup
 */
GEM_INLINE void filtering_candidates_init(filtering_candidates_t* const filtering_candidates) {
  // Candidates
  filtering_candidates->candidate_positions = vector_new(FC_INIT_CANDIDATE_POSITIONS,candidate_position_t);
  filtering_candidates->candidate_regions = vector_new(FC_INIT_CANDIDATE_POSITIONS,candidate_region_t);
  filtering_candidates->verified_regions = vector_new(FC_INIT_CANDIDATE_POSITIONS,verified_region_t);
  filtering_candidates->num_candidates_accepted = 0;
  // Region Buffer
  filtering_candidates->regions_buffer = vector_new(FC_INIT_REGIONS_BUFFER,region_t);
}
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
  vector_clear(filtering_candidates->candidate_positions);
  vector_clear(filtering_candidates->candidate_regions);
  vector_clear(filtering_candidates->verified_regions);
  filtering_candidates->num_candidates_accepted = 0;
  vector_clear(filtering_candidates->regions_buffer);
}
GEM_INLINE void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates) {
  vector_delete(filtering_candidates->candidate_positions);
  vector_delete(filtering_candidates->candidate_regions);
  vector_delete(filtering_candidates->verified_regions);
  vector_delete(filtering_candidates->regions_buffer);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t filtering_candidates_get_pending_candidate_regions(filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->candidate_regions);
}
/*
 * Adding candidate positions
 */
GEM_INLINE void filtering_candidates_add_interval(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t interval_lo,const uint64_t interval_hi,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t region_errors) {
  const uint64_t num_candidates = interval_hi-interval_lo;
  if (gem_expect_false(num_candidates==0)) return;
  // Store region
  region_t* region;
  vector_alloc_new(filtering_candidates->regions_buffer,region_t,region);
  region->start = region_start_pos;
  region->end = region_end_pos;
  region->degree = region_errors;
  // Store candidate positions (index-space)
  vector_reserve_additional(filtering_candidates->candidate_positions,num_candidates);
  candidate_position_t* candidate_position_index =
      vector_get_free_elm(filtering_candidates->candidate_positions,candidate_position_t);
  uint64_t index_position;
  for (index_position=interval_lo;index_position<interval_hi;++index_position) {
    candidate_position_index->candidate_region = region;
    candidate_position_index->candidate_region_index_position = index_position;
    ++candidate_position_index;
  }
  vector_add_used(filtering_candidates->candidate_positions,num_candidates);
}
GEM_INLINE void filtering_candidates_add_interval_set(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
GEM_INLINE void filtering_candidates_add_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t max_error) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
//  uint64_t it;
//  interval_t* result_interval = (interval_t*)vector_get_mem(result_vector) + init_int;
//  for (it=init_int; it<end_int; ++it, ++result_interval) {
//    if (result_interval->misms <= num_misms) {
//      ADD_INTERVAL_TO_FILTER_QUERIES(buffer_queries, result_interval->lo,
//        result_interval->hi, start_pos, end_pos, result_interval->misms);
//    }
//  }
// }
GEM_INLINE uint64_t filtering_candidates_get_pending_candidates(filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->candidate_positions);
}
/*
 * Candidate Accessors
 */
GEM_INLINE text_trace_t* filtering_candidate_get_text_trace(
    const text_collection_t* const candidates_collection,const candidate_region_t* const candidate_region) {
  return text_collection_get_trace(candidates_collection,candidate_region->candidate_text_trace_offset);
}
/*
 * Matching regions display
 */
GEM_INLINE void candidate_regions_matching_regions_print(
    FILE* const stream,candidate_region_t* const candidate_regions,const uint64_t candidate_position) {
  const uint64_t begin_position = candidate_regions->candidate_effective_begin_position;
  fprintf(stream,"    #%lu -> [%lu,+%lu) \n",
      candidate_position,begin_position,candidate_regions->candidate_effective_end_position-begin_position);
  uint64_t j;
  for (j=0;j<candidate_regions->num_regions_matching;++j) {
    region_matching_t* const region_matching = candidate_regions->regions_matching + j;
    fprintf(stream,"      %lu.%lu) -> [%lu,%lu) ~> [+%lu,+%lu) \n",candidate_position,j,
        region_matching->read_begin,region_matching->read_end,
        region_matching->text_begin,region_matching->text_end);
  }
}
GEM_INLINE void filtering_candidates_matching_regions_print(
    FILE* const stream,filtering_candidates_t* const filtering_candidates) {
  int64_t i;
  fprintf(stream,"[GEM]>Matching.Regions\n");
  fprintf(stream,"  => Initial.Regions\n");
  const uint64_t num_regions = vector_get_used(filtering_candidates->regions_buffer);
  region_t* const regions = vector_get_mem(filtering_candidates->regions_buffer,region_t);
  for (i=num_regions-1;i>=0;--i) {
    fprintf(stream,"    #%lu -> [%lu,%lu) \n",num_regions-i-1,regions[i].end,regions[i].start);
  }
  fprintf(stream,"  => Matching.Regions\n");
  const uint64_t num_candidate_regions = vector_get_used(filtering_candidates->candidate_regions);
  candidate_region_t* const candidate_regions = vector_get_mem(filtering_candidates->candidate_regions,candidate_region_t);
  for (i=0;i<num_candidate_regions;++i) {
    candidate_regions_matching_regions_print(stream,candidate_regions+i,i);
  }
}
/*
 * Candidate Alignment
 */
GEM_INLINE void filtering_accepted_regions_align_region(
    const text_collection_t* const candidates_collection,candidate_region_t* const accepted_region,
    const alignment_model_t alignment_model,const bool* const allowed_enc,
    const swg_penalties_t* const swg_penalties,const strand_t search_strand,
    const pattern_t* const pattern,const uint8_t* const key,const uint64_t key_length,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Select Model
  if (accepted_region->candidate_align_distance==0 || alignment_model==alignment_model_none) {
    // Add exact match
    match_trace_t match_trace;
    matches_align_exact(matches,&match_trace,search_strand,key_length,
        accepted_region->candidate_text_trace_offset,accepted_region->candidate_begin_position,
        accepted_region->candidate_align_distance,accepted_region->candidate_align_match_column+1);
    matches_add_match_trace_t(matches,&match_trace,true);
  } else {
    // Candidate
    const text_trace_t* const text_trace = filtering_candidate_get_text_trace(candidates_collection,accepted_region);
    const uint8_t* const text = text_trace->text;
    switch (alignment_model) {
      case alignment_model_hamming: {
        // Add hamming match
        match_trace_t match_trace;
        matches_align_hamming(matches,&match_trace,search_strand,allowed_enc,key,key_length,
            accepted_region->candidate_text_trace_offset,accepted_region->candidate_begin_position,
            text+(accepted_region->candidate_begin_position-accepted_region->candidate_effective_begin_position));
        matches_add_match_trace_t(matches,&match_trace,true);
        break;
      }
      case alignment_model_levenshtein: {
        // Add levenshtein match
        match_trace_t match_trace;
        matches_align_levenshtein(matches,&match_trace,search_strand,key,&pattern->bpm_pattern,
            accepted_region->candidate_text_trace_offset,accepted_region->candidate_begin_position,
            accepted_region->candidate_align_distance,text,accepted_region->candidate_align_match_column+1,
            accepted_region->regions_matching,accepted_region->num_regions_matching,mm_stack);
        matches_add_match_trace_t(matches,&match_trace,true);
        break;
      }
      case alignment_model_gap_affine: {
        // Add Smith-Waterman-Gotoh match (affine-gap)
        match_trace_t match_trace;
        matches_align_smith_waterman_gotoh(matches,&match_trace,
            search_strand,swg_penalties,key,key_length,accepted_region->candidate_text_trace_offset,
            accepted_region->candidate_begin_position,accepted_region->candidate_align_distance,
            text,accepted_region->candidate_effective_end_position-accepted_region->candidate_begin_position,
            accepted_region->regions_matching,accepted_region->num_regions_matching,mm_stack);
        matches_add_match_trace_t(matches,&match_trace,true);
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
GEM_INLINE void filtering_accepted_regions_align(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Model
  const alignment_model_t alignment_model = search_actual_parameters->search_parameters->alignment_model;
  // Pattern
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  search_parameters_t* const search_parameters = search_actual_parameters->search_parameters;
  const bool* const allowed_enc = search_parameters->allowed_enc;
  // Traverse all accepted candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->candidate_regions,candidate_region,candidate_pos,candidate_region_t) {
    if (candidate_region->candidate_align_distance == FC_DISTANCE_EXCEED) continue;
    // Align region
    filtering_accepted_regions_align_region(candidates_collection,candidate_region,
        alignment_model,allowed_enc,&search_parameters->swg_penalties,
        search_strand,pattern,key,key_length,matches,mm_stack);
  }
}
/*
 * Candidate Verify
 */
GEM_INLINE void filtering_candidate_region_verify_hamming(
    candidate_region_t* const candidate_region,const uint8_t* const text,
    const uint8_t* const key,const uint64_t key_length,
    const bool* const allowed_enc,uint64_t* const hamming_distance,const uint64_t max_mismatches) {
  // Check candidate
  uint64_t i, mismatches;
  for (i=0,mismatches=0;i<key_length;++i) {
    const uint8_t candidate_enc = text[i];
    // Check Mismatch
    if (!allowed_enc[candidate_enc] || candidate_enc != key[i]) {
      // Check Real Mismatch
      if (++mismatches > max_mismatches) {
        *hamming_distance = FC_DISTANCE_EXCEED;
        return;
      }
    }
  }
  *hamming_distance = mismatches;
  return;
}
GEM_INLINE void filtering_candidates_extend_matching_regions(
    const uint8_t* const key,const uint64_t key_length,
    candidate_region_t* const candidate_region,const uint8_t* const text,
    const bool* const allowed_enc) {
  // Check coverage
  if (candidate_region->coverage == key_length) return; // 100% coverage
  // Extend all matching regions
  const uint64_t num_regions_matching = candidate_region->num_regions_matching;
  const uint64_t effective_length = candidate_region->candidate_effective_end_position - candidate_region->candidate_effective_begin_position;
  const uint64_t last_region = num_regions_matching-1;
  uint64_t i, inc_coverage = 0;
  for (i=0;i<num_regions_matching;++i) {
    // Calculate limits
    region_matching_t* const region_matching = candidate_region->regions_matching + i;
    const int64_t left_read_max = (i==0) ? 0 : candidate_region->regions_matching[i-1].read_end+1;
    const int64_t right_read_max = (i==last_region) ? effective_length-1 : candidate_region->regions_matching[i+1].read_begin-1;
    const int64_t left_text_max = (i==0) ? 0 : candidate_region->regions_matching[i-1].text_end+1;
    const int64_t right_text_max = (i==last_region) ? effective_length-1 : candidate_region->regions_matching[i+1].text_begin-1;
    // Try to left extend
    int64_t left_read = region_matching->read_begin-1;
    int64_t left_text = region_matching->text_begin-1;
    while (left_read_max<=left_read && left_text_max<=left_text) {
      // Check match
      const uint8_t candidate_enc = text[left_text];
      if (!allowed_enc[candidate_enc] || candidate_enc != key[left_read]) break;
      --left_read;
      --left_text;
      ++inc_coverage;
    }
    region_matching->read_begin = left_read+1;
    region_matching->text_begin = left_text+1;
    // Try to right extend
    int64_t right_read = region_matching->read_end+1;
    int64_t right_text = region_matching->text_end+1;
    while (right_read_max>=right_read && right_text_max>=right_text) {
      // Check match
      const uint8_t candidate_enc = text[right_text];
      if (!allowed_enc[candidate_enc] || candidate_enc != key[right_read]) break;
      ++right_read;
      ++right_text;
      ++inc_coverage;
    }
    region_matching->read_end = right_read - 1;
    region_matching->text_end = right_text - 1;
  }
  candidate_region->coverage += inc_coverage;
  PROF_ADD_COUNTER(GP_FC_CANDIDATE_REGIONS_EXT_COVERAGE,(100*candidate_region->coverage)/key_length);
}
GEM_INLINE int64_t filtering_candidate_region_verify_levenshtein_debug(
    const uint8_t* const key,const uint64_t key_length,
    text_collection_t* const candidates_collection,
    const uint64_t candidate_pos,candidate_region_t* const candidate_region) {
  // Candidate
  const text_trace_t* const text_trace = filtering_candidate_get_text_trace(candidates_collection,candidate_region);
  const uint8_t* const text = text_trace->text;
  // Check
  gem_slog("(#%lu)[Checking position %lu]\n",candidate_pos,candidate_region->candidate_begin_position);
  gem_slog("\tRange [%lu,%lu]\n",
      candidate_region->candidate_begin_position,candidate_region->candidate_effective_end_position);
  gem_slog("\tEffective Range [%lu,%lu]\n",
      candidate_region->candidate_effective_begin_position,candidate_region->candidate_effective_end_position);
  const uint64_t eff_text_length =
      candidate_region->candidate_effective_end_position - candidate_region->candidate_effective_begin_position;
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
GEM_INLINE void filtering_candidate_region_verify_levenshtein(
    candidate_region_t* const candidate_region,const pattern_t* const pattern,
    const uint8_t* const candidate_text,const uint64_t candidate_length,
    uint64_t* const levenshtein_distance,uint64_t* const levenshtein_match_pos,const uint64_t max_error) {
  PROF_START(GP_FC_LEVENSHTEIN_BPM);
  // Cut-off + Quick-abandon
  bpm_get_distance__cutoff(
      &pattern->bpm_pattern,candidate_text,candidate_length,
      levenshtein_match_pos,levenshtein_distance,max_error);
  // DEBUG
  gem_cond_debug_block(DEBUG_ALIGN_LEVENSHTEIN) {
    uint64_t dp_position;
    const int64_t dp_distance = align_levenshtein_get_distance(
        (const char * const)pattern->key,pattern->key_length,
        (const char * const)candidate_text,candidate_length,true,&dp_position);
    if (dp_distance<=max_error && (dp_distance!=*levenshtein_distance || dp_position!=*levenshtein_match_pos)) {
      gem_slog("[BPM]Levenshtein-Alignment\n");
      // Print Pattern
      uint64_t i, j, p;
      gem_slog("\tPattern(%lu bases): ",pattern->key_length);
      for (i=0;i<pattern->key_length;++i) gem_slog("%c",dna_decode(pattern->key[i]));
      // Print BMP-Pattern
      const bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
      gem_slog("\n\tBPM-Pattern(%lu bases): ",bpm_pattern->pattern_length);
      for (i=0;i<bpm_pattern->pattern_length;++i) {
        const uint64_t block = i/8;
        const uint8_t mask = 1<<(i%8);
        for (j=0,p=0;j<DNA__N_RANGE;++j) {
          if (bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(j,block,
              bpm_pattern->pattern_num_words*bpm_pattern->pattern_word_size)] & mask) {
            gem_slog("%c",dna_decode(j));
            gem_cond_fatal_error_msg(p++>0,"More than 1 bit of PEQ are set on BPM-pattern");
          }
        }
        gem_cond_fatal_error_msg(p==0,"No bit of PEQ is set on BPM-pattern");
      }
      // Print Candidate
      gem_slog("\n\tCandidate: ");
      for (i=0;i<candidate_length;++i) gem_slog("%c",dna_decode(candidate_text[i]));
      gem_slog("\n\t\tDP-Alignment=(%lu,%lu)\tBPM-Alignment=(%lu,%lu)\n",
          dp_distance,dp_position,*levenshtein_distance,*levenshtein_match_pos);
      if (dp_distance<=max_error) {
        gem_cond_fatal_error_msg(dp_distance!=*levenshtein_distance,"Alignment distances don't match");
        gem_cond_fatal_error_msg(dp_position!=*levenshtein_match_pos,"Alignment positions don't match");
      }
    }
  }
  PROF_STOP(GP_FC_LEVENSHTEIN_BPM);
}
GEM_INLINE uint64_t filtering_candidate_region_verify(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,matches_t* const matches) {

  // TODO Remove discarded

  // Model
  const alignment_model_t alignment_model = search_actual_parameters->search_parameters->alignment_model;
  // Pattern
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const bool* const allowed_enc = search_actual_parameters->search_parameters->allowed_enc;
  // Matching Constraints
  const uint64_t max_search_matches = search_actual_parameters->search_parameters->max_search_matches;
  const uint64_t max_effective_filtering_error = pattern->max_effective_filtering_error;
  uint64_t total_accepted_regions = filtering_candidates->num_candidates_accepted;
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->candidate_regions,candidate_region,candidate_pos,candidate_region_t) {
    // Check matches accepted
    if (gem_expect_false(total_accepted_regions > max_search_matches)) break;
    // DEBUG
    gem_cond_debug_block(DEBUG_ALIGN_CANDIDATES) {
      filtering_candidate_region_verify_levenshtein_debug(key,key_length,candidates_collection,candidate_pos,candidate_region);
    }
    /*
     * 1. Hamming switch
     */
    const text_trace_t* const text_trace = filtering_candidate_get_text_trace(candidates_collection,candidate_region);
    const uint8_t* const text = text_trace->text; // Candidate
    if (alignment_model==alignment_model_hamming) {
      const uint64_t text_length = candidate_region->candidate_effective_end_position - candidate_region->candidate_begin_position;
      if (text_length >= key_length) {
        const uint64_t text_offset =
            candidate_region->candidate_begin_position - candidate_region->candidate_effective_begin_position;
        filtering_candidate_region_verify_hamming(candidate_region,text+text_offset,key,key_length,
            allowed_enc,&candidate_region->candidate_align_distance,max_effective_filtering_error);
        if (candidate_region->candidate_align_distance != FC_DISTANCE_EXCEED) ++total_accepted_regions;
      }
      continue; // Next
    }
    /*
     * 2. Extend regions
     */
//    filtering_candidates_extend_matching_regions(key,key_length,candidate_region,text,allowed_enc);
    /*
     * 3. Generalized Counting filter
     */
    const uint64_t eff_text_length =
        candidate_region->candidate_effective_end_position - candidate_region->candidate_effective_begin_position;
//    PROF_START(GP_FC_KMER_COUNTER_FILTER);
//    const uint64_t test_positive =
//        kmer_counting_filter(&pattern->kmer_counting,text,eff_text_length,max_effective_filtering_error);
//    PROF_STOP(GP_FC_KMER_COUNTER_FILTER);
//    if (test_positive==UINT64_MAX) {
//      PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
//    }
//    uint64_t n;
//    printf("+READ\n+");
//    for (n=0;n<pattern->key_length;++n) {
//      printf("%c",dna_decode(pattern->key[n]));
//    }
//    printf("\n+PATTERN\n+");
//    for (n=0;n<eff_text_length;++n) {
//      printf("%c",dna_decode(text[n]));
//    }
//    printf("\n");

    /*
     * 4. Myers's BPM algorithm
     */
    filtering_candidate_region_verify_levenshtein(candidate_region,pattern,text,eff_text_length,
        &candidate_region->candidate_align_distance,&candidate_region->candidate_align_match_column,
        max_effective_filtering_error);
    if (candidate_region->candidate_align_distance != FC_DISTANCE_EXCEED) {
//      if (test_positive==UINT64_MAX) fprintf(stderr,"Counting filter fails\n");
      PROF_INC_COUNTER(GP_FC_LEVENSHTEIN_ACCEPTED);
      candidate_region->candidate_begin_position = candidate_region->candidate_effective_begin_position;
      ++total_accepted_regions;
      continue; // Next
    }
    /*
     * 5. Local match (TODO Store as chunk matching ... try to extend borders)
     */
  }
  // Update
  const uint64_t accepted_regions = total_accepted_regions-filtering_candidates->num_candidates_accepted;
  filtering_candidates->num_candidates_accepted = total_accepted_regions;
  // DEBUG
  gem_cond_debug_block(DEBUG_REGIONS_MATCHING) {
    filtering_candidates_matching_regions_print(stderr,filtering_candidates);
  }
  // Return
  return accepted_regions;
}
/*
 * Retrieve all candidates(text) from the index
 */
GEM_INLINE void filtering_candidates_retrieve_candidate_regions(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const locator_t* const locator,const dna_text_t* const enc_text,mm_stack_t* const mm_stack) {
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->candidate_regions,text_candidate,candidate_pos,candidate_region_t) {
    // Skip discarded candidates
    if (text_candidate->candidate_begin_position == FC_POSITION_DISCARDED) continue;
    // Retrieve text(s)
    const uint64_t text_length = text_candidate->candidate_effective_end_position - text_candidate->candidate_effective_begin_position;
    archive_text_retrieve(locator,NULL,enc_text,candidates_collection,
        text_candidate->candidate_effective_begin_position,text_length,
        &text_candidate->candidate_text_trace_offset,mm_stack);
  }
}
/*
 * Filtering adjustment of the position wrt region/seed on which the candidate is based
 */
GEM_INLINE void filtering_candidates_adjust_position(
    candidate_position_t* const candidate,
    const uint64_t begin_offset,const uint64_t end_offset,const uint64_t boundary_error) {
  // Adjust Position
  const locator_interval_t* const locator_interval = candidate->locator_interval;
  uint64_t begin_position = (candidate->candidate_region_text_position > begin_offset) ?
      (candidate->candidate_region_text_position - begin_offset) : 0;
  uint64_t effective_begin_position;
  if (begin_position < locator_interval->begin_position) { // Adjust by locator-interval
    begin_position = locator_interval->begin_position; // Possible trim at the beginning
    effective_begin_position = locator_interval->begin_position;
  } else {
    effective_begin_position = (begin_position > boundary_error) ? begin_position-boundary_error : 0;
    if (effective_begin_position < locator_interval->begin_position) { // Adjust by locator-interval
      effective_begin_position = locator_interval->begin_position;
    }
  }
  uint64_t effective_end_position = candidate->candidate_region_text_position + end_offset + boundary_error;
  if (effective_end_position >= locator_interval->end_position) { // Adjust by locator-interval
    effective_end_position = locator_interval->end_position; // Possible trim at the end
  }
  candidate->candidate_begin_position = begin_position;
  candidate->candidate_effective_begin_position = effective_begin_position;
  candidate->candidate_effective_end_position = effective_end_position;
}
/*
 * Decode of all candidate positions (index-space -> text-space)
 */
GEM_INLINE void filtering_candidates_decode_candidates_positions(
    const locator_t* const locator,const fm_index_t* const fm_index,
    vector_t* const candidate_text_positions,const uint64_t key_length,const uint64_t boundary_error) {
  // Traverse all candidate positions in index-space
  VECTOR_ITERATE(candidate_text_positions,candidate,n,candidate_position_t) {
    // Lookup Position
    candidate->candidate_region_text_position = fm_index_lookup(fm_index,candidate->candidate_region_index_position);
    // Locate Position
    candidate->locator_interval = locator_lookup_interval(locator,candidate->candidate_region_text_position);
    // Adjust Position
    filtering_candidates_adjust_position(candidate,
        candidate->candidate_region->end,key_length-candidate->candidate_region->end,boundary_error);
  }
}
/*
 * Batch decode of all candidate positions (index-space -> text-space)
 *   (All the steps (CSA-lookup, rankQueries) are performed with prefetch-loops)
 */
//typedef struct {
//  uint64_t vector_rank;
//  uint64_t index_position;
//  uint64_t distance;
//  uint64_t used_slot;
//  bwt_block_locator_t bwt_block_locator;
//} fc_batch_decode_candidate;
//GEM_INLINE void filtering_candidates_decode_candidates_positions_batch_prefetched(
//    const locator_t* const locator,const fm_index_t* const fm_index,
//    vector_t* const candidate_text_positions,const uint64_t key_length,const uint64_t boundary_error) {
//  // Init
//  const uint64_t bwt_length = fm_index_get_length(fm_index);
//  const bwt_t* const bwt = fm_index->bwt;
//  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
//  candidate_position_t* const candidates = vector_get_mem(candidate_text_positions,candidate_position_t);
//  // Batch Decode
//  fc_batch_decode_candidate batch[FC_DECODE_NUM_POSITIONS_PREFETCHED];
//  const uint64_t num_candidate_text_positions = vector_get_used(candidate_text_positions);
//  // Initial fill batch
//  uint64_t current_position=0, i;
//  for (i=0;i<FC_DECODE_NUM_POSITIONS_PREFETCHED && current_position<num_candidate_text_positions;++current_position) {
//    if (!sampled_sa_is_sampled(sampled_sa,candidates[current_position].candidate_region_index_position)) {
//      batch[i].index_position = candidates[current_position].candidate_region_index_position;
//      batch[i].vector_rank = current_position;
//      batch[i].distance = 0;
//      batch[i].used_slot = true;
//      ++i;
//    }
//  }
//  const bool full_filled_batch = (i==FC_DECODE_NUM_POSITIONS_PREFETCHED);
//  for (;i<FC_DECODE_NUM_POSITIONS_PREFETCHED;++i) {
//    batch[i].used_slot = false;
//  }
//  // Full-prefetch loop for sampled-LF
//  if (full_filled_batch==FC_DECODE_NUM_POSITIONS_PREFETCHED) {
//    while (current_position<num_candidate_text_positions) {
//      for (i=0;i<FC_DECODE_NUM_POSITIONS_PREFETCHED;++i) {
//        bwt_prefetch(bwt,batch[i].index_position,&(batch[i].bwt_block_locator));
//      }
//      for (i=0;i<FC_DECODE_NUM_POSITIONS_PREFETCHED;++i) {
//        batch[i].index_position = bwt_prefetched_LF(bwt,batch[i].index_position,&(batch[i].bwt_block_locator));
//        ++(batch[i].distance);
//        if (sampled_sa_is_sampled(sampled_sa,batch[i].index_position)) {
//          candidates[batch[i].vector_rank].candidate_decode_sampled_pos = batch[i].index_position;
//          candidates[batch[i].vector_rank].candidate_decode_distance = batch[i].distance;
//          batch[i].used_slot = false;
//          // Select new candidate to decode
//          while (current_position < num_candidate_text_positions &&
//                 !sampled_sa_is_sampled(sampled_sa,candidates[current_position].candidate_region_index_position)) {
//            ++current_position;
//          }
//          if (current_position < num_candidate_text_positions) {
//            batch[i].index_position = candidates[current_position].candidate_region_index_position;
//            batch[i].vector_rank = current_position;
//            batch[i].distance = 0;
//            batch[i].used_slot = true;
//          }
//        }
//      }
//    }
//  }
//  // Solve remaining queries
//  for (i=0;i<FC_DECODE_NUM_POSITIONS_PREFETCHED;++i) {
//    if (batch[i].used_slot) {
//      do {
//        batch[i].index_position = bwt_LF(bwt,batch[i].index_position);
//        ++(batch[i].distance);
//      } while (!sampled_sa_is_sampled(sampled_sa,batch[i].index_position));
//      candidates[batch[i].vector_rank].candidate_decode_sampled_pos = batch[i].index_position;
//      candidates[batch[i].vector_rank].candidate_decode_distance = batch[i].distance;
//    }
//  }
//  // Prefetch SA-retrieve samples
//  uint64_t num_left_positions = num_candidate_text_positions;
//  current_position = 0;
//  while (num_left_positions < num_candidate_text_positions) {
//    const uint64_t batch_size = MIN(num_left_positions,FC_RETRIEVE_SAMPLE_NUM_POSITIONS_PREFETCHED);
//    const uint64_t batch_top = current_position+batch_size;
//    for (i=current_position;i<batch_top;++i) {
//      sampled_sa_prefetch_sample(sampled_sa,candidates[i].candidate_decode_sampled_pos);
//    }
//    for (i=current_position;i<batch_top;++i) {
//      candidates[i].candidate_region_text_position =
//          (sampled_sa_get_sample(sampled_sa,candidates[i].candidate_decode_sampled_pos) + candidates[i].candidate_decode_distance) % bwt_length;
//    }
//    current_position = batch_top;
//    num_left_positions -= batch_size;
//  }
//  // Adjust decoded position to the beginning of the read
//  candidate_position_t* candidate = vector_get_mem(candidate_text_positions,candidate_position_t);
//  for (current_position=0;current_position<num_candidate_text_positions;++current_position,++candidate) {
//    // Locate Position
//    candidates->locator_interval = locator_lookup_interval(locator,candidates->candidate_region_text_position);
//    // Adjust Position
//    filtering_candidates_adjust_position(candidates,
//        candidates->candidate_region->end,key_length-candidates->candidate_region->end,boundary_error);
//  }
//}
/*
 * Sorting candidates/positions/matching-regions
 */
int candidate_positions_cmp_position(const candidate_position_t* const a,const candidate_position_t* const b) {
  return a->candidate_begin_position - b->candidate_begin_position;
}
int regions_matching_cmp_position(const region_matching_t* const a,const region_matching_t* const b) {
  return a->text_begin - b->text_begin;
}
int verified_candidate_positions_cmp(const uint64_t* const a,const uint64_t* const b) {
  return *a - *b;
}
GEM_INLINE void filtering_candidates_sort_candidate_positions(const filtering_candidates_t* const filtering_candidates) {
  // Sort candidates positions (candidate_position_t) wrt their position in the text
  qsort(vector_get_mem(filtering_candidates->candidate_positions,candidate_position_t),
      vector_get_used(filtering_candidates->candidate_positions),sizeof(candidate_position_t),
      (int (*)(const void *,const void *))candidate_positions_cmp_position);
}
GEM_INLINE void filtering_candidates_sort_regions_matching(const candidate_region_t* const candidate_region) {
  // Sort regions matching (region_matching_t) wrt their starting position in the text
  qsort(candidate_region->regions_matching,candidate_region->num_regions_matching,
      sizeof(region_matching_t),(int (*)(const void *,const void *))regions_matching_cmp_position);
}

/*
 * Filtering discard duplicates & add positions to the list of verified positions
 */
//GEM_INLINE void filtering_candidates_add_verified_positions(
//    filtering_candidates_t* const filtering_candidates,const uint64_t num_added_positions);
//GEM_INLINE uint64_t filtering_candidates_discard_duplicates__add_to_verified(
//    filtering_candidates_t* const filtering_candidates,const uint64_t max_delta_difference);
//GEM_INLINE uint64_t filtering_candidates_discard_duplicates(
//    filtering_candidates_t* const filtering_candidates,const uint64_t max_delta_difference);




/*
 * Compose matching regions
 */
GEM_INLINE candidate_region_t* filtering_candidates_matching_regions_create(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t first_candidate_idx,const uint64_t last_candidate_idx,mm_stack_t* const mm_stack) {
  // Fetch candidate-positions
  candidate_position_t* const candidate_positions = vector_get_mem(filtering_candidates->candidate_positions,candidate_position_t);
  const uint64_t num_regions_matching = last_candidate_idx-first_candidate_idx+1;
  // Allow new matching candidate-region
  candidate_region_t* candidate_region;
  vector_alloc_new(filtering_candidates->candidate_regions,candidate_region_t,candidate_region);
  candidate_region->candidate_begin_position = candidate_positions[first_candidate_idx].candidate_begin_position;
  candidate_region->candidate_effective_begin_position = candidate_positions[first_candidate_idx].candidate_effective_begin_position;
  candidate_region->candidate_effective_end_position = candidate_positions[last_candidate_idx].candidate_effective_end_position;
  candidate_region->regions_matching = mm_stack_calloc(mm_stack,num_regions_matching,region_matching_t,false);
  candidate_region->num_regions_matching = num_regions_matching;
  uint64_t i, coverage = 0;
  for (i=0;i<num_regions_matching;++i) {
    region_matching_t* const region_matching = candidate_region->regions_matching + i;
    candidate_position_t* const candidate_position = candidate_positions + first_candidate_idx + i;
    // Region error
    region_matching->error = candidate_position->candidate_region->degree;
    // Read coordinates [Inclusive]
    region_matching->read_begin = candidate_position->candidate_region->end;
    region_matching->read_end = candidate_position->candidate_region->start - 1;
    // Text coordinates (relative to the effective begin position) [Inclusive]
    const uint64_t region_length = region_matching->read_end - region_matching->read_begin + 1;
    region_matching->text_begin = candidate_position->candidate_region_text_position - candidate_region->candidate_effective_begin_position;
    region_matching->text_end = region_matching->text_begin + region_length - 1;
    // Next
    coverage += region_length;
  }
  candidate_region->coverage = coverage;
  return candidate_region;
}
GEM_INLINE void filtering_candidates_matching_regions_check_arrangement(
    const dna_text_t* const enc_text,candidate_region_t* const candidate_region,
    mm_stack_t* const stack) {
  const uint64_t num_regions_matching = candidate_region->num_regions_matching;
  uint64_t i, overlapping=0;
  region_matching_t* last_region_matching = NULL;
  for (i=0;i<num_regions_matching;++i) {
    region_matching_t* const region_matching = candidate_region->regions_matching + i;
    if (last_region_matching!=NULL && last_region_matching->read_end >= region_matching->read_end) {
      overlapping=1;
    }
    last_region_matching = region_matching;
  }
  if (overlapping) {
    candidate_regions_matching_regions_print(stderr,candidate_region,0);
    // Retrieve text(s)
    const uint64_t text_length =
        candidate_region->candidate_effective_end_position - candidate_region->candidate_effective_begin_position;
    uint8_t* const text = dna_text_retrieve_sequence(
        enc_text,candidate_region->candidate_effective_begin_position,text_length,stack);
    fprintf(stderr,">>>READ\n");
    uint64_t n;
    for (n=0;n<text_length;++n) {
      fprintf(stderr,"%c",dna_decode(text[n]));
    }
    fprintf(stderr,"\n");
  }
}
GEM_INLINE uint64_t filtering_candidates_matching_regions_compose(
    filtering_candidates_t* const filtering_candidates,const dna_text_t* const enc_text,
    const uint64_t key_length,const uint64_t max_delta_difference,mm_stack_t* const mm_stack) {
  // Sort candidate positions (text-space)
  filtering_candidates_sort_candidate_positions(filtering_candidates);
  // Traverse positions and eliminate duplicates
  const uint64_t num_candidate_positions = vector_get_used(filtering_candidates->candidate_positions);
  candidate_position_t* const candidate_positions = vector_get_mem(filtering_candidates->candidate_positions,candidate_position_t);
  uint64_t candidate_idx = 0;
  while (candidate_idx < num_candidate_positions) {
    // Determine the positions belonging to the same region
    uint64_t last_position = candidate_positions[candidate_idx].candidate_begin_position;
    uint64_t group_idx = candidate_idx + 1;
    while (group_idx < num_candidate_positions) {
      const uint64_t position = candidate_positions[group_idx].candidate_begin_position;
      const uint64_t delta = position - last_position;
      if (delta > max_delta_difference) break; // Doesn't belong to the group. Stop!
      // Next
      last_position = position;
      ++group_idx;
    }
    // Create a region candidate with the positions from [candidate_idx] to [group_idx-1]
    filtering_candidates_matching_regions_create(filtering_candidates,candidate_idx,group_idx-1,mm_stack);
    // PROF_ADD_COUNTER(GP_FC_CANDIDATE_REGIONS_COVERAGE,(100*candidate_region->coverage)/key_length);  // TODO
    // Sort matching regions
    // filtering_candidates_sort_regions_matching(candidate_region); // TODO
    // Check matching-regions arrangement
    // filtering_candidates_matching_regions_check_arrangement(enc_text,candidate_region,mm_stack); // TODO
    // Next group
    const uint64_t num_regions_matching = group_idx-candidate_idx;
    candidate_idx += num_regions_matching;
  }
  // Clear candidate positions
  vector_clear(filtering_candidates->candidate_positions);
  // DEBUG
  gem_cond_debug_block(DEBUG_REGIONS_MATCHING) {
    filtering_candidates_matching_regions_print(stderr,filtering_candidates);
  }

//  // Add to verified positions // TODO
//  uint64_t candidate_idx = 0, verified_idx = 0;
//  const uint64_t num_verified_positions = vector_get_used(filtering_candidates->verified_positions);
//  const uint64_t* const verified_positions = vector_get_mem(filtering_candidates->verified_positions,uint64_t);
//  // Merge verified candidate positions with accepted positions
//  filtering_candidates_add_verified_positions(filtering_candidates,num_accepted_positions);

  // Return number of accepted positions
  const uint64_t num_accepted_positions = vector_get_used(filtering_candidates->candidate_regions);
  return num_accepted_positions;
}
/*
 * Verify filtering candidates
 */
GEM_INLINE uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,const pattern_t* const pattern,
    const search_actual_parameters_t* const search_actual_parameters,mm_stack_t* const mm_stack) {
  PROF_START(GP_FC_PROCESS_CANDIDATES);
  // Check non-empty pending candidates set
  uint64_t pending_candidates = filtering_candidates_get_pending_candidates(filtering_candidates);
  PROF_ADD_COUNTER(GP_FC_NUM_CANDIDATE_POSITIONS,pending_candidates);
  if (pending_candidates==0) {
    PROF_STOP(GP_FC_PROCESS_CANDIDATES); return 0; // Nothing to do
  }

  // Batch decode+adjust of all positions of the candidates (cip_begin_position = decoded(candidate_region_index_position))
  PROF_START(GP_FC_DECODE_POSITIONS);
  const uint64_t key_length = pattern->key_length;
  const uint64_t boundary_error = search_actual_parameters->max_filtering_error_nominal;
//  if (pending_candidates < FC_DECODE_NUM_POSITIONS_PREFETCHED) {
    filtering_candidates_decode_candidates_positions(
        locator,fm_index,filtering_candidates->candidate_positions,key_length,boundary_error);
//  } else {  // TODO Enable batch decode
//    filtering_candidates_decode_candidates_positions_batch_prefetched(
//        locator,fm_index,filtering_candidates->candidate_text_positions,key_length,boundary_error);
//  }
  PROF_STOP(GP_FC_DECODE_POSITIONS);

  // Compose matching regions into candidate regions (also filter out duplicated positions or already checked)
  PROF_START(GP_FC_COMPOSE_REGIONS);
  pending_candidates = filtering_candidates_matching_regions_compose(
      filtering_candidates,enc_text,key_length,boundary_error,mm_stack);
  PROF_STOP(GP_FC_COMPOSE_REGIONS);
  PROF_ADD_COUNTER(GP_FC_NUM_CANDIDATE_REGIONS,pending_candidates);

  // Return total candidate regions
  PROF_STOP(GP_FC_PROCESS_CANDIDATES);
  return pending_candidates;
}
GEM_INLINE uint64_t filtering_candidates_verify(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const text_collection,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  PROF_START(GP_FC_VERIFICATION);

  // Process candidates
  uint64_t pending_candidates = filtering_candidates_process_candidates(filtering_candidates,
      locator,fm_index,enc_text,pattern,search_actual_parameters,mm_stack);
  if (pending_candidates==0) { PROF_STOP(GP_FC_VERIFICATION); return 0; }

  // Retrieve text-candidates
  PROF_START(GP_FC_RETRIEVE_CANDIDATE_REGIONS);
  filtering_candidates_retrieve_candidate_regions(filtering_candidates,text_collection,locator,enc_text,mm_stack);
  PROF_STOP(GP_FC_RETRIEVE_CANDIDATE_REGIONS);

  // Verify candidates
  PROF_START(GP_FC_VERIFY_CANDIDATE_REGIONS);
  pending_candidates = filtering_candidate_region_verify(
      filtering_candidates,text_collection,pattern,search_strand,search_actual_parameters,matches);
  PROF_ADD_COUNTER(GP_FC_NUM_ACCEPTED_REGIONS,pending_candidates);
  PROF_STOP(GP_FC_VERIFY_CANDIDATE_REGIONS);
  if (pending_candidates==0) { PROF_STOP(GP_FC_VERIFICATION); return 0; }

  // Align accepted candidates
  PROF_START(GP_FC_REALIGN_CANDIDATE_REGIONS);
  matches_hint_add_match_trace(matches,pending_candidates); // Hint to matches
  filtering_accepted_regions_align(filtering_candidates,text_collection,
      pattern,search_strand,search_actual_parameters,matches,mm_stack);
  PROF_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS);

  PROF_STOP(GP_FC_VERIFICATION);
  return pending_candidates;
}
/*
 * BPM-Buffer API (Verification)
 */
GEM_INLINE uint64_t filtering_candidates_bpm_buffer_add(
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,mm_stack_t* const mm_stack) {
  // Check number of pending regions
  const uint64_t pending_candidates = filtering_candidates_get_pending_candidate_regions(filtering_candidates);
  if (pending_candidates==0) return 0;
  // Add the pattern to the buffer (add a new query)
  bpm_gpu_buffer_put_pattern(bpm_gpu_buffer,pattern);
  // Traverse all candidates (text-space) & add them to the buffer
  const candidate_region_t* candidate_region = vector_get_mem(filtering_candidates->candidate_regions,candidate_region_t);
  uint64_t candidate_pos;
  for (candidate_pos=0;candidate_pos<pending_candidates;++candidate_pos,++candidate_region) {
    const uint64_t eff_candidate_begin_position = candidate_region->candidate_effective_begin_position;
    const uint64_t eff_text_length = candidate_region->candidate_effective_end_position - eff_candidate_begin_position;
    bpm_gpu_buffer_put_candidate(bpm_gpu_buffer,eff_candidate_begin_position,eff_text_length);
  }
  // Return the final number of candidates added to the buffer
  return pending_candidates;
}
GEM_INLINE void filtering_candidates_bpm_buffer_align(
    const text_collection_t* const text_collection,const dna_text_t* const enc_text,pattern_t* const pattern,
    const strand_t search_strand,const search_actual_parameters_t* const search_actual_parameters,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end,
    matches_t* const matches,mm_stack_t* const mm_stack) {
  // Count total candidates
  const uint64_t total_candidates = candidate_offset_end-candidate_offset_begin;
  if (gem_expect_false(total_candidates==0)) return;
  PROF_START(GP_FC_REALIGN_CANDIDATE_REGIONS);
  // Hint to matches
  matches_hint_add_match_trace(matches,total_candidates);
  // Fetch Parameters
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_effective_filtering_error = pattern->max_effective_filtering_error;
  search_parameters_t* const search_parameters = search_actual_parameters->search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  const bool* const allowed_enc = search_parameters->allowed_enc;
  // Traverse all candidates
  uint64_t i;
  for (i=candidate_offset_begin;i<candidate_offset_end;++i) {
    uint32_t levenshtein_distance, levenshtein_match_pos;
    bpm_gpu_buffer_get_candidate_result(bpm_gpu_buffer,i,&levenshtein_distance,&levenshtein_match_pos);
//    // DEBUG
//    uint32_t candidate_text_position, candidate_length;
//    bpm_gpu_buffer_get_candidate(bpm_gpu_buffer,i,&candidate_text_position,&candidate_length);
//    fprintf(stderr,"F p=%lu e=%lu c=%lu\n",candidate_text_position,levenshtein_distance,levenshtein_match_pos);
    if (levenshtein_distance <= max_effective_filtering_error) {
      // Get the accepted candidate
      uint64_t candidate_text_position;
      uint32_t candidate_length;
      bpm_gpu_buffer_get_candidate(bpm_gpu_buffer,i,&candidate_text_position,&candidate_length);
      // Allocate text-trace
      const uint64_t text_trace_offset = text_collection_new_trace(text_collection);
      text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
      text_trace->text = dna_text_retrieve_sequence(enc_text,candidate_text_position,candidate_length,mm_stack);
      text_trace->length = candidate_length;
      // Configure accepted candidate (DTO)
      candidate_region_t accepted_region;
      accepted_region.candidate_text_trace_offset = text_trace_offset;
      accepted_region.candidate_begin_position = candidate_text_position;
      accepted_region.candidate_effective_begin_position = candidate_text_position;
      accepted_region.candidate_effective_end_position = candidate_text_position+candidate_length;
      accepted_region.candidate_align_distance = levenshtein_distance;
      accepted_region.candidate_align_match_column = levenshtein_match_pos;
      accepted_region.regions_matching = NULL;
      accepted_region.num_regions_matching = 0;
      filtering_accepted_regions_align_region(text_collection,&accepted_region,
          alignment_model,allowed_enc,&search_parameters->swg_penalties,
          search_strand,pattern,key,key_length,matches,mm_stack);
    }
  }
  PROF_STOP(GP_FC_REALIGN_CANDIDATE_REGIONS);
}


