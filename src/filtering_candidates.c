/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"

/*
 * Debug
 */
#define DEBUG_ALIGN_CANDIDATES  false
#define DEBUG_ALIGN_LEVENSHTEIN true

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
 * Candidates
 */
// (1) Alias input fields
#define cip_region_index_position          overloaded_field_0
// (2) Alias decoding positions
#define cip_distance                       overloaded_field_0
#define cip_sampled_pos                    overloaded_field_1
#define cip_region_text_position           overloaded_field_1
// (3) Alias re-aligning candidates
#define cip_eff_candidate_end_position     overloaded_field_0
#define cip_eff_candidate_begin_position   overloaded_field_1
#define cip_candidate_begin_position       overloaded_field_2
#define cip_text_trace_offset              overloaded_field_3
//// (4) Alias output filtering
//#define cip_alg_distance                   overloaded_field_2
//#define cip_alg_position                   overloaded_field_3
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
 * Setup
 */
GEM_INLINE void filtering_candidates_init(filtering_candidates_t* const filtering_candidates) {
  // Candidates Positions
  filtering_candidates->candidate_positions = vector_new(FC_INIT_CANDIDATE_POSITIONS,candidate_position_t);
  // Checked Positions
  filtering_candidates->verified_candidate_positions = vector_new(FC_INIT_CANDIDATE_POSITIONS,uint64_t);
  // Candidates accepted
  filtering_candidates->num_candidates_accepted = 0;
  filtering_candidates->max_candidates_accepted = ALL;
  // Internals
  filtering_candidates->regions_buffer = vector_new(FC_INIT_REGIONS_BUFFER,region_t);
}
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
  vector_clear(filtering_candidates->regions_buffer);
  vector_clear(filtering_candidates->candidate_positions);
  filtering_candidates->num_candidates_accepted = 0;
  filtering_candidates->max_candidates_accepted = ALL;
  vector_clear(filtering_candidates->verified_candidate_positions);
}
GEM_INLINE void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates) {
  vector_delete(filtering_candidates->regions_buffer);
  vector_delete(filtering_candidates->candidate_positions);
  vector_delete(filtering_candidates->verified_candidate_positions);
}
/*
 * Max candidates accepted
 */
GEM_INLINE void filtering_candidates_set_max_candidates_accepted(
    filtering_candidates_t* const filtering_candidates,const uint64_t max_candidates_accepted) {
  filtering_candidates->max_candidates_accepted = max_candidates_accepted;
}
GEM_INLINE bool filtering_candidates_is_max_candidates_reached(filtering_candidates_t* const filtering_candidates) {
  return filtering_candidates->num_candidates_accepted >= filtering_candidates->max_candidates_accepted;
}
/*
 * Add Candidates
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
    candidate_position_index->cip_region_index_position = index_position;
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
    const text_collection_t* const candidates_collection,const candidate_position_t* const text_candidate) {
  return text_collection_get_trace(candidates_collection,text_candidate->cip_text_trace_offset);
}
/*
 *
 */
GEM_INLINE void filtering_candidates_align_hamming(
    candidate_position_t* const text_candidate,const uint8_t* const text,
    const uint8_t* const key,const uint64_t key_length,
    const bool* const allowed_enc,uint64_t* const hamming_distance,const uint64_t max_mismatches) {
  // Check candidate
  uint64_t i, mismatches;
  for (i=0,mismatches=0;i<key_length;++i) {
    const uint64_t candidate_enc = text[i];
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
GEM_INLINE void filtering_candidates_align_levenshtein(
    candidate_position_t* const text_candidate,const pattern_t* const pattern,
    const uint8_t* const candidate_text,const uint64_t candidate_length,
    uint64_t* const levenshtein_distance,uint64_t* const levenshtein_match_pos,const uint64_t max_error) {
  // TODO: Cut-off + Quick-abandon (if diagonal cannot reach max_error, then quit())
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
}
GEM_INLINE void filtering_candidates_align(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,matches_t* const matches) {
  // Pattern
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  // Matching Constraints
  const uint64_t max_search_matches = search_actual_parameters->search_parameters->max_search_matches;
//  const bool* const allowed_enc = search_actual_parameters->search_parameters->allowed_enc;
  const uint64_t max_effective_filtering_error = pattern->max_effective_filtering_error;
  uint64_t num_candidates_accepted = filtering_candidates->num_candidates_accepted;
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->candidate_positions,text_candidate,candidate_pos,candidate_position_t) {
    // Skip discarded candidates
    if (text_candidate->cip_candidate_begin_position == FC_POSITION_DISCARDED) continue;
    // Check matches accepted
    if (gem_expect_false(num_candidates_accepted > max_search_matches)) break;
    // Candidate
    const text_trace_t* const text_trace = filtering_candidate_get_text_trace(candidates_collection,text_candidate);
    const uint8_t* const text = text_trace->text;
    /*
     * DEBUG
     */
    gem_cond_debug_block(DEBUG_ALIGN_CANDIDATES) {
      gem_slog("(#%lu)[Checking position %lu]\n",candidate_pos,text_candidate->cip_candidate_begin_position);
      gem_slog("\tRange [%lu,%lu]\n",
          text_candidate->cip_candidate_begin_position,text_candidate->cip_eff_candidate_end_position);
      gem_slog("\tEffective Range [%lu,%lu]\n",
          text_candidate->cip_eff_candidate_begin_position,text_candidate->cip_eff_candidate_end_position);
      const uint64_t eff_text_length =
          text_candidate->cip_eff_candidate_end_position - text_candidate->cip_eff_candidate_begin_position;
      uint64_t i, dp_position;
      const int64_t dp_distance = align_levenshtein_get_distance(
          (const char * const)key,key_length,(const char * const)text,eff_text_length,true,&dp_position);
      gem_slog("\tDP-Alignment (distance=%lu,position=%lu)\n",dp_distance,dp_position);
      gem_slog("\tPattern: ");
      for (i=0;i<key_length;++i) gem_slog("%c",dna_decode(key[i]));
      gem_slog("\n\tText: ");
      for (i=0;i<eff_text_length;++i) gem_slog("%c",dna_decode(text[i]));
      gem_slog("\n");
    }
//    // 1. Check Hamming distance
//    const uint64_t text_length =
//        text_candidate->cip_eff_candidate_end_position - text_candidate->cip_candidate_begin_position;
//    if (text_length >= key_length) {
//      const uint64_t text_offset =
//          text_candidate->cip_candidate_begin_position - text_candidate->cip_eff_candidate_begin_position;
//      // TODO Check this mathematically embedding distances => 1d
//      uint64_t hamming_distance;
//      filtering_candidates_align_hamming(text_candidate,text+text_offset,
//          key,key_length,allowed_enc,&hamming_distance,max_effective_filtering_error);
//      if (hamming_distance != FC_DISTANCE_EXCEED) {
//        // Store match
//        matches_add_match_trace_mark(
//            matches,text_candidate->cip_text_trace_offset,
//            text_candidate->cip_candidate_begin_position,hamming_distance,
//            text_offset,key_length,search_strand,true);
//        ++num_candidates_accepted;
//        continue; // Next
//      }
//    }
    // 2. Check Levenshtein distance
    // 2.1 Generalized Counting filter // TODO
    // 2.2 Myers's BPM algorithm
    const uint64_t eff_text_length =
        text_candidate->cip_eff_candidate_end_position - text_candidate->cip_eff_candidate_begin_position;
    uint64_t levenshtein_distance, levenshtein_match_pos;
    filtering_candidates_align_levenshtein(
        text_candidate,pattern,text,eff_text_length,
        &levenshtein_distance,&levenshtein_match_pos,max_effective_filtering_error);
    if (levenshtein_distance != FC_DISTANCE_EXCEED) {
      // Store match
      matches_add_match_trace_mark(
          matches,text_candidate->cip_text_trace_offset,
          text_candidate->cip_eff_candidate_begin_position,levenshtein_distance,
          0,levenshtein_match_pos,search_strand,true);
      ++num_candidates_accepted;
      continue; // Next
    }

    // 3. Local match (TODO Store as chunk matching ... try to extend borders)
  }
  // Update
  filtering_candidates->num_candidates_accepted = num_candidates_accepted;


//          if (delta_position > position) { // Initial trim
//            initial_trim = true;
//            initial_trim_len = delta_position-position;
//            offset_key += initial_trim_len;
//            offset_key_len -= initial_trim_len;
//          }
//          if (delta_position+delta_decoded_len < position+key_len) { //  Final trim
//            final_trim = true;
//            final_trim_len = (position+key_len) - (delta_position+delta_decoded_len);
//            offset_key_len -= final_trim_len;
//          }

//          if (offset_key_len<(int64_t)min_anchor_size) continue; // No enough anchor

//
//    //
//    // Check single-indel distance
//    //
//    START_TIMER(TSC_CHECK_BIG_INDEL);
//    const idx_t indel_delta = queries->f_gap_length;
//    const uint64_t indel_key_len = key_len+2*indel_delta;
//    if (indel_delta==0 || indel_delta<=max_distance || indel_delta>max_indel_len ||
//        key_len < indel_key_len || position < indel_delta ||
//        position+key_len+indel_delta >= a->text_length) continue;
//    INC_COUNTER(GSC_CHECK_BIG_INDEL);
//    if (fmi_decode(a,position-indel_delta,
//        indel_key_len,decoded)!=indel_key_len) continue; // Gap found !
//    correct = fmi_check_single_indel_candidate(a,key,position,decoded+indel_delta,
//        key_len,indel_delta,allowed_chars,mismatch_mask,max_mismatches,
//        &match_position,&match_end_position,&match_mismatches,misms);
//    STOP_TIMER(TSC_CHECK_BIG_INDEL);
//    if (correct) {  // CORRECT SBI, cool!!
//      INC_COUNTER(GSC_BIG_INDEL);
//      FMI_MATCHES_CHECK_DUPLICATES__STORE_MATCH(match_position,
//          (match_end_position-match_position),match_mismatches);
//    }
//  }
}
/*
 * Retrieve all candidates(text) from the index
 */
GEM_INLINE void filtering_candidates_retrieve_candidates(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const candidates_collection,
    const locator_t* const locator,const dna_text_t* const enc_text) {
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->candidate_positions,text_candidate,candidate_pos,candidate_position_t) {
    // Skip discarded candidates
    if (text_candidate->cip_candidate_begin_position == FC_POSITION_DISCARDED) continue;
    // Allocate text-trace
    const uint64_t text_trace_offset = text_collection_new_trace(candidates_collection);
    text_candidate->cip_text_trace_offset = text_trace_offset; // Link it with the candidate
    text_trace_t* const text_trace = text_collection_get_trace(candidates_collection,text_trace_offset);
    text_trace->text = dna_text_get_buffer(enc_text) + text_candidate->cip_eff_candidate_begin_position;
    PREFETCH(text_trace->text); // Prefetch text // TODO Hint later on (LLC)
    text_trace->length = text_candidate->cip_eff_candidate_end_position - text_candidate->cip_eff_candidate_begin_position;

//    // Allocate trace-block [[ TODO GRAPH]]
//    const uint64_t trace_block_offset = text_collection_allocate_trace_blocks(candidates_collection,1);
//    text_trace->trace_blocks_offset = trace_block_offset;
//    text_trace->trace_length = 1;
//    trace_block_t* const trace_block = text_collection_get_trace_block(candidates_collection,trace_block_offset);
//    trace_block->position = text_candidate->position;
//    trace_block->length = text_candidate->effective_text_length;

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
  uint64_t begin_position = (candidate->cip_region_text_position > begin_offset) ?
      (candidate->cip_region_text_position - begin_offset) : 0;
  uint64_t effective_begin_position;
  if (begin_position < locator_interval->begin_position) { // Adjust by locator-interval
    begin_position = locator_interval->begin_position; // Possible trim at the beginning
    effective_begin_position = locator_interval->begin_position;
  } else {
    effective_begin_position = (begin_position > boundary_error) ? begin_position-boundary_error : 0;
    if (begin_position < locator_interval->begin_position) { // Adjust by locator-interval
      effective_begin_position = locator_interval->begin_position;
    }
  }
  uint64_t effective_end_position = candidate->cip_region_text_position + end_offset + boundary_error;
  if (effective_end_position >= locator_interval->end_position) { // Adjust by locator-interval
    effective_end_position = locator_interval->begin_position; // Possible trim at the end
  }
  candidate->cip_candidate_begin_position = begin_position;
  candidate->cip_eff_candidate_begin_position = effective_begin_position;
  candidate->cip_eff_candidate_end_position = effective_end_position;
}
/*
 * Filters all the regions candidates in @positions against the index.
 * Check whether this regions belong to an alignment with errors (mismatches/indels)
 */
GEM_INLINE void filtering_candidates_decode_candidates_positions(
    const locator_t* const locator,const fm_index_t* const fm_index,
    vector_t* const candidate_text_positions,const uint64_t key_length,const uint64_t boundary_error) {
  // Traverse all candidate positions in index-space
  VECTOR_ITERATE(candidate_text_positions,candidate,n,candidate_position_t) {
    // Lookup Position
    candidate->cip_region_text_position = fm_index_lookup(fm_index,candidate->cip_region_index_position);
    // Locate Position
    candidate->locator_interval = locator_lookup_interval(locator,candidate->cip_region_text_position);
    // Adjust Position
    filtering_candidates_adjust_position(candidate,
        candidate->candidate_region->end,key_length-candidate->candidate_region->end,boundary_error);
  }
}
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
//    if (!sampled_sa_is_sampled(sampled_sa,candidates[current_position].cip_region_index_position)) {
//      batch[i].index_position = candidates[current_position].cip_region_index_position;
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
//          candidates[batch[i].vector_rank].cip_sampled_pos = batch[i].index_position;
//          candidates[batch[i].vector_rank].cip_distance = batch[i].distance;
//          batch[i].used_slot = false;
//          // Select new candidate to decode
//          while (current_position < num_candidate_text_positions &&
//                 !sampled_sa_is_sampled(sampled_sa,candidates[current_position].cip_region_index_position)) {
//            ++current_position;
//          }
//          if (current_position < num_candidate_text_positions) {
//            batch[i].index_position = candidates[current_position].cip_region_index_position;
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
//      candidates[batch[i].vector_rank].cip_sampled_pos = batch[i].index_position;
//      candidates[batch[i].vector_rank].cip_distance = batch[i].distance;
//    }
//  }
//  // Prefetch SA-retrieve samples
//  uint64_t num_left_positions = num_candidate_text_positions;
//  current_position = 0;
//  while (num_left_positions < num_candidate_text_positions) {
//    const uint64_t batch_size = MIN(num_left_positions,FC_RETRIEVE_SAMPLE_NUM_POSITIONS_PREFETCHED);
//    const uint64_t batch_top = current_position+batch_size;
//    for (i=current_position;i<batch_top;++i) {
//      sampled_sa_prefetch_sample(sampled_sa,candidates[i].cip_sampled_pos);
//    }
//    for (i=current_position;i<batch_top;++i) {
//      candidates[i].cip_region_text_position =
//          (sampled_sa_get_sample(sampled_sa,candidates[i].cip_sampled_pos) + candidates[i].cip_distance) % bwt_length;
//    }
//    current_position = batch_top;
//    num_left_positions -= batch_size;
//  }
//  // Adjust decoded position to the beginning of the read
//  candidate_position_t* candidate = vector_get_mem(candidate_text_positions,candidate_position_t);
//  for (current_position=0;current_position<num_candidate_text_positions;++current_position,++candidate) {
//    // Locate Position
//    candidates->locator_interval = locator_lookup_interval(locator,candidates->cip_region_text_position);
//    // Adjust Position
//    filtering_candidates_adjust_position(candidates,
//        candidates->candidate_region->end,key_length-candidates->candidate_region->end,boundary_error);
//  }
//}
/*
 * Sorting candidates/positions
 */
int candidate_positions_cmp_position(const candidate_position_t* const a,const candidate_position_t* const b) {
  return a->cip_candidate_begin_position - b->cip_candidate_begin_position;
}
int verified_candidate_positions_cmp(const uint64_t* const a,const uint64_t* const b) {
  return *a - *b;
}
GEM_INLINE void filtering_candidates_sort_candidate_positions(const filtering_candidates_t* const filtering_candidates) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(filtering_candidates->candidate_positions,candidate_position_t),
      vector_get_used(filtering_candidates->candidate_positions),sizeof(candidate_position_t),
      (int (*)(const void *,const void *))candidate_positions_cmp_position);
}
/*
 * Filtering discard duplicates & add positions to the list of verified positions
 */
GEM_INLINE uint64_t filtering_candidates_get_previous_valid_candidate_position(
    const candidate_position_t* const candidates_vector,int64_t* const candidates_vector_idx) {
  if (*candidates_vector_idx < 0) return -1;
  // Load current position
  uint64_t candidates_vector_position = candidates_vector[*candidates_vector_idx].cip_candidate_begin_position;
  // Skip discarded positions of candidates_vector
  while (candidates_vector_position == FC_POSITION_DISCARDED) {
    --(*candidates_vector_idx); // Load next
    if (*candidates_vector_idx < 0) return -1;
    candidates_vector_position = candidates_vector[*candidates_vector_idx].cip_candidate_begin_position;
  }
  // Return
  return candidates_vector_position;
}
GEM_INLINE void filtering_candidates_add_verified_positions(
    filtering_candidates_t* const filtering_candidates,const uint64_t num_added_positions) {
  // Check zero added
  if (num_added_positions == 0) return;
  // Vectors
  vector_t* const candidate_positions = filtering_candidates->candidate_positions;
  vector_t* const verified_candidate_positions = filtering_candidates->verified_candidate_positions;
  // Reserve added ones
  vector_reserve_additional(verified_candidate_positions,num_added_positions);
  // Prepare vector Sentinels
  int64_t candidates_vector_idx = vector_get_used(candidate_positions)-1;
  int64_t verified_vector_idx = vector_get_used(verified_candidate_positions)-1;
  int64_t merged_vector_idx = candidates_vector_idx+num_added_positions;
  candidate_position_t* const candidates_vector = vector_get_mem(candidate_positions,candidate_position_t);
  uint64_t* const verified_vector = vector_get_mem(verified_candidate_positions,uint64_t);
  // Positions Sentinels (Load current positions)
  uint64_t candidates_vector_position =
      filtering_candidates_get_previous_valid_candidate_position(candidates_vector,&candidates_vector_idx);
  uint64_t verified_vector_position = (verified_vector_idx>=0) ? verified_vector[verified_vector_idx] : -1;
  while (candidates_vector_idx >= 0 && verified_vector_idx >= 0) {
    // Add to merged_vector
    if (candidates_vector_position > verified_vector_position) {
      verified_vector[merged_vector_idx--] = candidates_vector_position;
      // Load previous position
      if (--candidates_vector_idx < 0) break;
      candidates_vector_position =
          filtering_candidates_get_previous_valid_candidate_position(candidates_vector,&candidates_vector_idx);
    } else { // candidates_vector_position < verified_vector_position
      GEM_INTERNAL_CHECK(candidates_vector_position!=verified_vector_position,"Merging equal positions (invalid)");
      verified_vector[merged_vector_idx--] = verified_vector_position;
      // Load previous position
      if (--verified_vector_idx < 0) break;
      verified_vector_position = (verified_vector_idx>=0) ? verified_vector[verified_vector_idx] : -1;
    }
  }
  // Add remaining
  while (verified_vector_idx >= 0) {
    verified_vector[merged_vector_idx--] = verified_vector[verified_vector_idx--];
  }
  while (candidates_vector_idx >= 0) {
    verified_vector[merged_vector_idx--] = candidates_vector_position;
    --candidates_vector_idx;
    candidates_vector_position =
        filtering_candidates_get_previous_valid_candidate_position(candidates_vector,&candidates_vector_idx);
  }
  // Set used
  vector_add_used(verified_candidate_positions,num_added_positions);
}
GEM_INLINE uint64_t filtering_candidates_discard_duplicates__add_to_verified(
    filtering_candidates_t* const filtering_candidates,const uint64_t max_delta_difference) {
  // Sort candidate positions (text-space)
  filtering_candidates_sort_candidate_positions(filtering_candidates);
  // Traverse positions and eliminate duplicates
  vector_t* const candidate_positions = filtering_candidates->candidate_positions;
  vector_t* const verified_candidate_positions = filtering_candidates->verified_candidate_positions;
  const uint64_t num_verified_positions = vector_get_used(verified_candidate_positions);
  const uint64_t* const verified_positions = vector_get_mem(verified_candidate_positions,uint64_t);
  uint64_t last_position = UINT64_MAX, verified_idx = 0, num_accepted_positions = 0;
  VECTOR_ITERATE(candidate_positions,text_candidate,n,candidate_position_t) {
    // Check the position for duplicates (against previous position)
    const uint64_t position = text_candidate->cip_candidate_begin_position;
    const uint64_t delta = last_position<=position ? position-last_position : UINT64_MAX;
    if (delta <= max_delta_difference) { // TODO delta=0 => Discard & delta>0 => Extend check
      text_candidate->cip_candidate_begin_position = FC_POSITION_DISCARDED;
      continue; // Repeated position
    }
    // Check the position for duplicates (against previous verified positions)
    while (verified_idx < num_verified_positions) {
      const int64_t delta = ((int64_t)position - (int64_t)verified_positions[verified_idx]);
      if (ABS(delta) <= max_delta_difference) {
        text_candidate->cip_candidate_begin_position = FC_POSITION_DISCARDED;
        continue; // Already checked position
      }
      if (delta > 0) ++verified_idx; // Next
    }
    // Accept the candidate
    ++num_accepted_positions;
    last_position = position;
  }
  // Merge verified candidate positions with accepted positions
  filtering_candidates_add_verified_positions(filtering_candidates,num_accepted_positions);
  // Return number of accepted positions
  return num_accepted_positions;
}
GEM_INLINE uint64_t filtering_candidates_discard_duplicates(
    filtering_candidates_t* const filtering_candidates,const uint64_t max_delta_difference) {
  // Sort candidate positions (text-space)
  filtering_candidates_sort_candidate_positions(filtering_candidates);
  // Traverse positions and eliminate duplicates
  vector_t* const candidate_positions = filtering_candidates->candidate_positions;
  uint64_t last_position = UINT64_MAX, num_accepted_positions = 0;
  VECTOR_ITERATE(candidate_positions,text_candidate,n,candidate_position_t) {
    // Check the position for duplicates (against previous position)
    const uint64_t position = text_candidate->cip_candidate_begin_position;
    const uint64_t delta = last_position<=position ? position-last_position : UINT64_MAX;
    if (delta <= max_delta_difference) { // TODO delta=0 => Discard & delta>0 => Extend check
      text_candidate->cip_candidate_begin_position = FC_POSITION_DISCARDED;
      continue; // Repeated position
    }
    // Accept the candidate
    ++num_accepted_positions;
    last_position = position;
  }
  // Return number of accepted positions
  return num_accepted_positions;
}
/*
 * Batch decode of all candidate positions (index-space -> text-space)
 *   (All the steps (CSA-lookup, rankQueries) are performed with prefetch-loops)
 */
GEM_INLINE void filtering_candidates_verify_pending(
    filtering_candidates_t* const filtering_candidates,text_collection_t* const text_collection,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,matches_t* const matches) {

  // Check non-empty pending candidates set
  uint64_t pending_candidates = filtering_candidates_get_pending_candidates(filtering_candidates);
  if (pending_candidates==0) return; // Nothing to do

  // Batch decode+adjust of all positions of the candidates (cip_begin_position = decoded(cip_region_index_position))
  const uint64_t key_length = pattern->key_length;
  const uint64_t boundary_error = search_actual_parameters->max_search_error_nominal;
//  if (pending_candidates < FC_DECODE_NUM_POSITIONS_PREFETCHED) {
    filtering_candidates_decode_candidates_positions(
        locator,fm_index,filtering_candidates->candidate_positions,key_length,boundary_error);
//  } else {  // TODO Enable batch decode
//    filtering_candidates_decode_candidates_positions_batch_prefetched(
//        locator,fm_index,filtering_candidates->candidate_text_positions,key_length,boundary_error);
//  }

  // Filter out duplicated positions (or already checked)
  pending_candidates = filtering_candidates_discard_duplicates__add_to_verified(filtering_candidates,boundary_error);
  if (pending_candidates==0) return;

  // Retrieve text-candidates
  filtering_candidates_retrieve_candidates(filtering_candidates,text_collection,locator,enc_text);

  // Verify candidates
  matches_hint_add_match_trace(matches,pending_candidates); // Hint to matches
  filtering_candidates_align(filtering_candidates,text_collection,pattern,search_strand,search_actual_parameters,matches);
}
GEM_INLINE uint64_t filtering_candidates_add_to_bpm_buffer(
    filtering_candidates_t* const filtering_candidates,
    const locator_t* const locator,const fm_index_t* const fm_index,
    const dna_text_t* const enc_text,const pattern_t* const pattern,const strand_t search_strand,
    const search_actual_parameters_t* const search_actual_parameters,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  // Batch decode+adjust of all positions of the candidates (cip_begin_position = decoded(cip_region_index_position))
  const uint64_t key_length = pattern->key_length;
  const uint64_t boundary_error = search_actual_parameters->max_search_error_nominal;
//  if (pending_candidates < FC_DECODE_NUM_POSITIONS_PREFETCHED) {
    filtering_candidates_decode_candidates_positions(
        locator,fm_index,filtering_candidates->candidate_positions,key_length,boundary_error);
//  } else {  // TODO Enable batch decode
//    filtering_candidates_decode_candidates_positions_batch_prefetched(
//        locator,fm_index,filtering_candidates->candidate_text_positions,key_length,boundary_error);
//  }

  // Filter out duplicated positions
  const uint64_t pending_candidates =
      filtering_candidates_discard_duplicates(filtering_candidates,boundary_error);
  if (pending_candidates==0) return 0;

  // Add the pattern to the buffer (add a new query)
  bpm_gpu_buffer_put_pattern(bpm_gpu_buffer,(bpm_pattern_t* const)&pattern->bpm_pattern);
  // Traverse all candidates (text-space) & add them to the buffer
  const uint64_t num_candidates = vector_get_used(filtering_candidates->candidate_positions);
  candidate_position_t* text_candidate = vector_get_mem(filtering_candidates->candidate_positions,candidate_position_t);
  uint64_t candidate_pos;
  for (candidate_pos=0;candidate_pos<num_candidates;++candidate_pos,++text_candidate) {
    const uint64_t eff_candidate_begin_position = text_candidate->cip_eff_candidate_begin_position;
    const uint64_t eff_text_length = text_candidate->cip_eff_candidate_end_position - eff_candidate_begin_position;
    bpm_gpu_buffer_put_candidate(bpm_gpu_buffer,eff_candidate_begin_position,eff_text_length);
  }

  // Clear candidates (makes no sense to keep storing anything)
  filtering_candidates_clear(filtering_candidates);

  // Return the final number of candidates added to the buffer
  return num_candidates;
}



