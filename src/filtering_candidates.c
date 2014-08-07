/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "locator.h"

/*
 * Constants
 */
#define FC_INIT_REGIONS_BUFFER               100
#define FC_INIT_CANDIDATE_POSITIONS         1000

#define FC_DECODE_NUM_POSITIONS_PREFETCHED          10
#define FC_RETRIEVE_SAMPLE_NUM_POSITIONS_PREFETCHED 10

#define FC_DISTANCE_EXCEED UINT64_MAX

/*
 * Filtering Candidates
 */
struct _filtering_candidates_t {
  /* Pending candidates */
  vector_t* candidate_positions;              // Candidates positions (candidate_position_t)
  /* Checked Positions */
  vector_t* verified_candidate_positions;     // Verified positions (uint64_t)
  /* Internals */
  vector_t* regions_buffer;                   // Regions Buffer (region_t)
  text_collection_t* candidates_collection;   // Candidates Text-Collection (Stores candidates Texts)
};

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
// (4) Alias output filtering
#define cip_alg_distance                   overloaded_field_2
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
GEM_INLINE void filtering_candidates_new(filtering_candidates_t* const filtering_candidates) {
  // Candidates Positions
  filtering_candidates->candidate_positions = vector_new(FC_INIT_CANDIDATE_POSITIONS,candidate_position_t);
  // Checked Positions
  filtering_candidates->verified_candidate_positions = vector_new(FC_INIT_CANDIDATE_POSITIONS,uint64_t);
  // Internals
  filtering_candidates->regions_buffer = vector_new(FC_INIT_REGIONS_BUFFER,region_t);
  filtering_candidates->candidates_collection = text_collection_new();
}
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
  vector_clear(filtering_candidates->regions_buffer);
  vector_clear(filtering_candidates->candidate_positions);
  text_collection_clear(filtering_candidates->candidates_collection);
  vector_clear(filtering_candidates->verified_candidate_positions);
}
GEM_INLINE void filtering_candidates_delete(filtering_candidates_t* const filtering_candidates) {
  vector_delete(filtering_candidates->regions_buffer);
  vector_delete(filtering_candidates->candidate_positions);
  text_collection_delete(filtering_candidates->candidates_collection);
  vector_delete(filtering_candidates->verified_candidate_positions);
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
    const uint8_t* const allowed_enc,const uint64_t max_mismatches) {
  // Check candidate
  uint64_t i, mismatches;
  for (i=0,mismatches=0;i<key_length;++i) {
    const uint64_t candidate_enc = text[i];
    // Check Mismatch
    if (!allowed_enc[candidate_enc] || candidate_enc != key[i]) {
      // Check Real Mismatch
      if (++mismatches > max_mismatches) {
        text_candidate->cip_alg_distance = FC_DISTANCE_EXCEED;
        return;
      }
    }
  }
  text_candidate->cip_alg_distance = mismatches;
  return;
}

GEM_INLINE void filtering_candidates_align_levenshtein(
    candidate_position_t* const text_candidate,pattern_t* const pattern,
    const uint8_t* const candidate_text,const uint64_t candidate_length,const uint64_t max_mismatches) {
  uint64_t position, distance;
  bpm_pattern_t* const bpm_pattern = &pattern->bpm_pattern;
  bpm_get_distance__cutoff(
      bpm_pattern,candidate_text,candidate_length,
      &position,&distance,max_mismatches);
  // TODO
}
GEM_INLINE void filtering_candidates_align(
    approximate_search_t* const approximate_search,matches_t* const matches,
    filtering_candidates_t* const filtering_candidates) {
  const approximate_search_parameters_t* const search_parameters = approximate_search->search_parameters;
  const text_collection_t* const candidates_collection = filtering_candidates->candidates_collection;
  // Pattern
  const pattern_t* const pattern = &approximate_search->pattern;
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  // Matching Constraints
  const uint8_t* const allowed_enc = search_parameters->allowed_enc;
  const uint64_t max_effective_filtering_error = pattern->max_effective_filtering_error;
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->candidate_positions,text_candidate,candidate_pos,candidate_position_t) {
    // Candidate
    const text_trace_t* const text_trace = filtering_candidate_get_text_trace(candidates_collection,text_candidate);
    const uint8_t* const text = text_trace->text;
    // 1. Check Hamming distance
    const uint64_t text_length =
        text_candidate->cip_eff_candidate_end_position - text_candidate->cip_candidate_begin_position;
    if (text_length >= key_length) {
      const uint64_t boundary_offset =
          text_candidate->cip_candidate_begin_position - text_candidate->cip_eff_candidate_begin_position;
      filtering_candidates_align_hamming(text_candidate,text+boundary_offset,key,key_length,allowed_enc,1); // TODO Check this mathematically
//      if (distance != FC_DISTANCE_EXCEED) { // TODO
//        // Store match
//        continue; // Next
//      }
    }
    // 2. Check Levenshtein distance
    // 2.1 Generalized Counting filter // TODO
    // 2.2 Myers's BPM algorithm
    const uint64_t eff_text_length =
        text_candidate->cip_eff_candidate_end_position - text_candidate->cip_eff_candidate_begin_position;
    filtering_candidates_align_levenshtein(pattern,text,eff_text_length,max_effective_filtering_error);
//    if (distance != FC_DISTANCE_EXCEED) { // TODO
//      // Store match
//      continue; // Next
//    }

    // 3. Local match (TODO Store as chunk matching ... try to extend borders)
  }


//  ch_t *key = search_params->key;
//  const uint64_t* const peq = search_params->peq;
//  const uint64_t key_len = search_params->key_len;
//  ch_t* mismatch_mask = search_params->mismatch_mask;
//  const uint64_t max_indel_len = search_params->max_indel_len;
//  const bool* const allowed_chars = search_params->allowed_chars;
//  const uint64_t min_anchor_size = search_params->min_anchor_size;
//  // const bool unique_mapping = search_params->unique_mapping; // TODO
//
//  START_STATS(SC_CHECK_CAND);
//
//  uint64_t* cntrs = vector_get_mem(matches->rbuf_counts);
//  pos_match_t* out = (pos_match_t*)vector_get_mem(matches->rbuf_pos)+vector_get_used(matches->rbuf_pos);
//  mismatch* misms = (mismatch*)vector_get_mem(matches->rbuf_mismatches)+vector_get_used(matches->rbuf_mismatches);
//
//  const uint64_t max_pattern_ext = max_indel_len+UNSAFE_MAX(max_distance,max_mismatches);
//  vector_prepare(mpool->fbuf2, ch_t, 2*max_pattern_ext+key_len); // We use fbuf2 to store the decoded read
//  ch_t* decoded = vector_get_mem(mpool->fbuf2);
//
//  uint64_t num_matches, i;
//  uint64_t match_mismatches;
//  idx_t match_position, match_end_position;
//  for (i = 0, num_matches=0; i < num_queries; ++i, ++queries) {
//    INC_COUNTER(GSC_CHECKED_MATCHES);
//
//    //
//    // Check Hamming distance
//    //
//    const idx_t position = queries->cip_begin_position;
//    bool correct = true;
//    if (position+key_len < a->text_length) {
//      // Decode
//      START_STATS(SC_CHECK_DECODE);
//      idx_t delta_position = position-UNSAFE_MIN(max_distance,queries->cip_begin_position);
//      idx_t delta = position-delta_position;
//      idx_t delta_end_key = UNSAFE_MIN(position+key_len+max_distance,a->text_length);
//      idx_t delta_key_len = delta_end_key-delta_position;
//      idx_t delta_decoded_len = fmi_decode(a,delta_position,delta_key_len,decoded);
//      bool exists_gap = false;
//      if (__builtin_expect(delta_decoded_len < key_len+delta,false)) { // There is a gap (Separator)
//        exists_gap = true;
//        const idx_t gap_pos = delta_position+delta_decoded_len;
//        if (delta_decoded_len < delta) { // Delta gap, so still can try hamming
//          delta_position = gap_pos+1; // Correct the position
//          delta_end_key = UNSAFE_MIN(delta_position+key_len+max_distance,a->text_length);
//          delta_key_len = delta_end_key-delta_position;
//          delta = position-delta_position;
//          assert(((int64_t)position-(int64_t)delta_position)>=0);
//          delta_decoded_len = fmi_decode(a,delta_position,delta_key_len,decoded);
//          STOP_STATS(SC_CHECK_DECODE);
//          if (delta_decoded_len >= key_len+delta) {
//            START_TIMER(TSC_CHECK_HAMMING);
//            correct = fmi_matches_check_hamming_match(
//                key,key_len,mismatch_mask,allowed_chars,decoded+delta,
//                max_mismatches,misms,&match_mismatches,mpool);
//            STOP_TIMER(TSC_CHECK_HAMMING);
//          } else {
//            correct = false;
//          }
//        } else {
//          correct = false; // We go straight to levenshtein
//          const idx_t chunk_pos = position+queries->end;
//          if (gap_pos < chunk_pos) {
//            delta_position = gap_pos+1;
//            delta_end_key = UNSAFE_MIN(delta_position+key_len+max_distance,a->text_length);
//          } else { // gap_pos > chunk_pos
//            delta_end_key = gap_pos+1;
//          }
//          delta_key_len = delta_end_key-delta_position;
//          delta_decoded_len = fmi_decode(a,delta_position,delta_key_len,decoded);
//          STOP_STATS(SC_CHECK_DECODE);
//        }
//      } else {
//        STOP_STATS(SC_CHECK_DECODE); START_TIMER(TSC_CHECK_HAMMING);
//        correct = fmi_matches_check_hamming_match(
//            key,key_len,mismatch_mask,allowed_chars,decoded+delta,
//            max_mismatches,misms,&match_mismatches,mpool);
//        STOP_TIMER(TSC_CHECK_HAMMING);
//      }
//      if (correct) { // CORRECT Hamming, cool!!
//        INC_COUNTER(GSC_HAMMING_HIT);
//        if (store_hamming_matches) {
//          // FMI_MATCHES_STORE_NEW_MATCH(out,position,key_len,match_mismatches);
//          FMI_MATCHES_CHECK_DUPLICATES__STORE_MATCH(position,key_len,match_mismatches);
//        } else {
//          FMI_MATCHES_ADJUST_MAX_MISMS_USING_DELTA(match_mismatches);
//        }
//      } else if (max_distance > 0 && min_anchor_size > 0) {
//        //
//        // Check Levenshtein distance
//        //
//        if (delta_decoded_len>=key_len || (key_len-delta_decoded_len)<=max_distance) {
//          START_TIMER(TSC_CHECK_LEVENSHTEIN);
//          correct = fmi_matches_check_levenshtein_match(
//              key,peq,key_len,decoded,delta_decoded_len,max_distance,true,
//              misms,&match_position,&match_end_position,&match_mismatches,mpool);
//          STOP_TIMER(TSC_CHECK_LEVENSHTEIN);
//        } else if (exists_gap) {
//          bool initial_trim = false, final_trim = false;
//          idx_t initial_trim_len = 0, final_trim_len = 0;
//          ch_t* offset_key = key;
//          int64_t offset_key_len = key_len;
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
//          START_TIMER(TSC_CHECK_LEVENSHTEIN);
//          // Recompute the peq vector
//          uint64_t* const offset_peq = (uint64_t*)vector_get_mem(mpool->vbuf1)+(key_len*5);
//          fmi_prepare_pattern_eq(a,offset_key,offset_key_len,offset_peq,allowed_chars,mismatch_mask);
//          // Adjust the maximum differences allowed
//          uint64_t max_anchor_distance = (offset_key_len*search_params->max_differences)/key_len;
//          max_anchor_distance += matches->bad_quality_misms; // TODO: Tight
//          if (max_anchor_distance > offset_key_len) max_anchor_distance=offset_key_len;
//          correct = fmi_matches_check_levenshtein_match(
//              offset_key,offset_peq,offset_key_len,decoded,delta_decoded_len,max_anchor_distance,true,
//              misms,&match_position,&match_end_position,&match_mismatches,mpool);
//          STOP_TIMER(TSC_CHECK_LEVENSHTEIN);
//          if (correct) {
//            if (initial_trim) {
//              FMI_MATCHES_ADD_INITIAL_TRIM(misms,initial_trim_len,match_position,match_mismatches);
//            }
//            if (final_trim) {
//              const uint64_t total_len = offset_key_len+(initial_trim?initial_trim_len:0);
//              FMI_MATCHES_ADD_FINAL_TRIM(misms,final_trim_len,total_len,match_end_position,match_mismatches);
//            }
//          }
//        } else {
//          correct = false;
//        }
//        if (correct) { // CORRECT Levenshtein, cool!!
//          INC_COUNTER(GSC_LEVENSHTEIN);
//          match_position+=delta_position; match_end_position+=delta_position;
//          FMI_MATCHES_CHECK_DUPLICATES__STORE_MATCH(match_position,
//              (match_end_position-match_position),match_mismatches);
//        }
//      }
//    }
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
//
//  // Update mismatches and positions counters
//  vector_update_used(matches->rbuf_mismatches,misms);
//  vector_update_used(matches->rbuf_pos,out);
//
//  STOP_STATS(SC_CHECK_CAND);
}
/*
 *
 */
GEM_INLINE void filtering_candidates_retrieve_candidates(
    approximate_search_t* const approximate_search,filtering_candidates_t* const filtering_candidates) {
  const locator_t* const locator = approximate_search->locator;
  const dna_text_t* const enc_text = approximate_search->enc_text;
  text_collection_t* const candidates_collection = filtering_candidates->candidates_collection;
  // Traverse all candidates (text-space)
  VECTOR_ITERATE(filtering_candidates->candidate_positions,candidate,candidate_pos,candidate_position_t) {
    // Allocate text-trace
    const uint64_t text_trace_offset = text_collection_new_trace(candidates_collection);
    candidate->cip_text_trace_offset = text_trace_offset; // Link it with the candidate
    text_trace_t* const text_trace = text_collection_get_trace(candidates_collection,text_trace_offset);
    text_trace->text = dna_text_get_buffer(enc_text) + candidate->cip_eff_candidate_begin_position;
    PREFETCH(text_trace->text); // Prefetch text // TODO Hint later on (LLC)
    text_trace->length = candidate->cip_eff_candidate_end_position - candidate->cip_eff_candidate_begin_position;

//    // Allocate trace-block [[GRAPH]]
//    const uint64_t trace_block_offset = text_collection_allocate_trace_blocks(candidates_collection,1);
//    text_trace->trace_blocks_offset = trace_block_offset;
//    text_trace->trace_length = 1;
//    trace_block_t* const trace_block = text_collection_get_trace_block(candidates_collection,trace_block_offset);
//    trace_block->position = text_candidate->position;
//    trace_block->length = text_candidate->effective_text_length;

  }
}
/*
 *
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
  uint64_t end_position = begin_position + end_offset + boundary_error;
  if (end_position >= locator_interval->end_position) { // Adjust by locator-interval
    end_position = locator_interval->begin_position; // Possible trim at the end
  }
  candidate->cip_candidate_begin_position = begin_position;
  candidate->cip_eff_candidate_begin_position = effective_begin_position;
  candidate->cip_eff_candidate_end_position = end_position;
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
 *
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
GEM_INLINE void filtering_candidates_sort_verified_candidate_positions(filtering_candidates_t* const filtering_candidates) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(filtering_candidates->verified_candidate_positions,uint64_t),
      vector_get_used(filtering_candidates->verified_candidate_positions),sizeof(uint64_t),
      (int (*)(const void *,const void *))verified_candidate_positions_cmp);
}
/*
 *
 */
GEM_INLINE uint64_t filtering_candidates_discard_duplicates(
    approximate_search_t* const approximate_search,filtering_candidates_t* const filtering_candidates) {

  // Sort candidate positions (text-space)
  filtering_candidates_sort_candidate_positions(filtering_candidates);

  // Traverse positions and eliminate duplicates
  const uint64_t pending_candidates = vector_get_used(filtering_candidates->candidate_positions);
  candidate_position_t* const candidates =
      vector_get_mem(filtering_candidates->candidate_positions,candidate_position_t);
  uint64_t* verified_candidate_position =
      vector_get_mem(filtering_candidates->verified_candidate_positions,uint64_t);
  uint64_t num_accepted_positions = 0;
  uint64_t i, checked_idx;
  uint64_t last_in_pos = UINT64_MAX, last_checked_pos = UINT64_MAX;
  for (i=0;i<pending_candidates;++i) {


    // Check the position for duplicates (against previous position)
    const uint64_t position = candidates[i].cip_candidate_begin_position;
    const uint64_t delta = last_in_pos<=position ? position-last_in_pos : UINT64_MAX;
    if (delta==0) continue; // Repeated position // FIXME: Delta threshold
    // Check the position for duplicates (against previous verified positions)
    if (ihash_is_contained(filtering_candidates->verified_candidate_text_positions,position)) continue; // Repeated position
    text_candidates[num_queries].position = position;
    text_candidates[num_queries].ctp_effective_text_length =
        index_candidates[i].cip_eff_candidate_end_position - index_candidates[i].cip_eff_candidate_begin_position;
    text_candidates[num_queries].candidate_region = index_candidates[i].candidate_region;

    ++num_accepted_positions;
    last_in_pos = position;


  }
  vector_set_used(filtering_candidates->candidate_positions,num_accepted_positions);

  // Sort verified candidate positions
  filtering_candidates_sort_verified_candidate_positions(filtering_candidates);

  return num_accepted_positions;
}
/*
 * Batch decode of all candidate positions (index-space -> text-space)
 *   (All the steps (CSA-lookup, rankQueries) are performed with prefetch-loops)
 */
GEM_INLINE void filtering_candidates_verify_pending(
    approximate_search_t* const approximate_search,matches_t* const matches,
    filtering_candidates_t* const filtering_candidates) {

  // Check non-empty pending candidates set
  uint64_t pending_candidates = filtering_candidates_get_pending_candidates(filtering_candidates);
  if (pending_candidates==0) return; // Nothing to do

  // Batch decode+adjust of all positions of the candidates (cip_begin_position = decoded(cip_region_index_position))
  const fm_index_t* const fm_index = approximate_search->fm_index;
  const locator_t* const locator = approximate_search->locator;
  const uint64_t key_length = approximate_search->pattern.key_length;
  const uint64_t boundary_error = approximate_search->search_parameters->max_search_error_nominal;
//  if (pending_candidates < FC_DECODE_NUM_POSITIONS_PREFETCHED) {
    filtering_candidates_decode_candidates_positions(
        locator,fm_index,filtering_candidates->candidate_positions,key_length,boundary_error);
//  } else {  // TODO Enable batch decode
//    filtering_candidates_decode_candidates_positions_batch_prefetched(
//        locator,fm_index,filtering_candidates->candidate_text_positions,key_length,boundary_error);
//  }

  // Filter out duplicated positions (or already checked)
  if ((pending_candidates=filtering_candidates_discard_duplicates(approximate_search,filtering_candidates))==0) return;

  // Retrieve text-candidates
  filtering_candidates_retrieve_candidates(approximate_search,filtering_candidates);

  // Verify candidates
  matches_hint_add_match_trace(matches,pending_candidates); // Hint to matches
  filtering_candidates_align(approximate_search,matches,filtering_candidates);
}



