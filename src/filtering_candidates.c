/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"

/*
 * Candidate Internals
 */
// General fields
#define end_key_pos    overload_field_1
#define end_reg_idx    overload_field_2
// Used in decode positions
#define f_distance     overload_field_0
#define f_sampled_pos  overload_field_1
#define f_csa_pos      overload_field_1
// Used in indel match check
#define f_gap_length   overload_field_0

/*
 * Setup
 */
GEM_INLINE void filtering_candidates_new(filtering_candidates_t* const filtering_candidates,mm_slab_t* const mm_slab) {
//  // Allocate vectors
//  filtering_candidates->pending_candidates = svector_new(mm_slab,candidate_t);
//  filtering_candidates->checked_positions = svector_new(mm_slab,uint64_t);
//  // Init appending iterator
//  svector_iterator_new(&filtering_candidates->pending_candidates_iterator,
//      filtering_candidates->pending_candidates,SVECTOR_WRITE_ITERATOR,0);
}
GEM_INLINE void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
//  // Clear vectors
//  svector_clear(filtering_candidates->pending_candidates);
//  svector_clear(filtering_candidates->checked_positions);
//  // Reset appending iterator
//  svector_iterator_new(&filtering_candidates->pending_candidates_iterator,
//      filtering_candidates->pending_candidates,SVECTOR_WRITE_ITERATOR,0);
}
GEM_INLINE void filtering_candidates_delete(filtering_candidates_t* const filtering_candidates) {
//  svector_delete(filtering_candidates->pending_candidates);
//  svector_delete(filtering_candidates->checked_positions);
}
/*
 * Add Candidates
 */
GEM_INLINE void filtering_candidates_add_text_position(
    filtering_candidates_t* const filtering_candidates,const uint64_t index_position) {
  // TODO
}
GEM_INLINE void filtering_candidates_add_interval(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t interval_lo,const uint64_t interval_hi,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t region_errors) {
//  svector_iterator_t* const pending_candidates_iterator = &filtering_candidates->pending_candidates_iterator;
//  uint64_t index_position;
//  for (index_position=interval_lo;index_position<interval_hi;++index_position) {
//    candidate_t* const candidate = svector_iterator_get_element(pending_candidates_iterator,candidate_t);
//    candidate->start = region_start_pos;
//    candidate->end = region_end_pos;
//    candidate->errors = region_errors;
//    candidate->end_reg_idx = index_position;
//    svector_write_iterator_next(pending_candidates_iterator);
//  }
}
GEM_INLINE void filtering_candidates_add_interval_set(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos) {
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
}
GEM_INLINE void filtering_candidates_add_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,interval_set_t* const interval_set,
    const uint64_t region_start_pos,const uint64_t region_end_pos,const uint64_t max_error) {
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
}
//  register uint64_t it;
//  register interval_t* result_interval = (interval_t*)vector_get_mem(result_vector) + init_int;
//  for (it=init_int; it<end_int; ++it, ++result_interval) {
//    if (result_interval->misms <= num_misms) {
//      ADD_INTERVAL_TO_FILTER_QUERIES(buffer_queries, result_interval->lo,
//        result_interval->hi, start_pos, end_pos, result_interval->misms);
//    }
//  }
// }
GEM_INLINE uint64_t filtering_candidates_get_pending_candidates(filtering_candidates_t* const filtering_candidates) {
  // TODO
  return 0;
}
/*
 * Filters all the regions candidates in @positions against the index.
 * Check whether this regions belong to an alignment with mismatches and, possibly,
 * a single indel (of a limited length)
 */
GEM_INLINE void filtering_candidates_verify_pending(
    filtering_candidates_t* const filtering_candidates,matches_t* const matches,
    const bool use_levenshtein_distance,const bool store_hamming_matches) {
//  register const uint64_t num_checked_positions = vector_get_used(checked_positions);
//
//  // Check non-empty positions set
//  if (vector_get_used(in_queries)==0) {
//    if (num_checked_positions > 0) {
//      vector_prepare(out_positions,idx_t,num_checked_positions);
//      memcpy(vector_get_mem(out_positions),
//          vector_get_mem(checked_positions),sizeof(idx_t)*num_checked_positions);
//      vector_set_used(out_positions,num_checked_positions);
//    }
//    return;
//  }
//
//  // Batch decode of all positions of the candidates => end_key_pos = decoded(end_reg_idx);
//  register const uint64_t key_len = search_params->key_len;
//  if (vector_get_used(in_queries) < NUMBER_DECODE_QUERIES_PREFETCHED) {
//    fmi_matches_decode_positions_candidates(a,key_len,matches,in_queries);
//  } else {
//    fmi_matches_prefetched_decode_positions_candidates(a,key_len,matches,in_queries);
//  }
//
//  // Sorting the positions
//  filter_query_t_sort(vector_get_mem(in_queries),vector_get_used(in_queries));
//
//  // Allocate space for the output positions vector
//  register const idx_t *checked_pos = vector_get_mem(checked_positions);
//  vector_prepare(out_positions,idx_t,vector_get_used(in_queries)+num_checked_positions);
//  register idx_t *out_pos = vector_get_mem(out_positions);
//
//  // Traverse positions and eliminate duplicates and store indel positions
//  register const uint64_t num_in_queries = vector_get_used(in_queries);
//  vector_prepare(mpool->buf3,filter_query,num_in_queries);
//  register filter_query* queries = vector_get_mem(in_queries);
//  register filter_query* pending_queries = vector_get_mem(mpool->buf3);
//  register uint64_t num_queries, i, checked_idx;
//  register idx_t last_in_pos = IDX_T_MAX, last_checked_pos = IDX_T_MAX;
//  for (i=0,num_queries=0,checked_idx=0; i<num_in_queries && queries[i].end_key_pos!=IDX_T_MAX; ++i) {
//    register const uint64_t pos = queries[i].end_key_pos;
//    register const uint64_t in_delta = last_in_pos<=pos ? pos-last_in_pos : IDX_T_MAX;
//    if (in_delta==0) continue; // Repeated position
//    // Synch with current query position
//    while (checked_idx<num_checked_positions && checked_pos[checked_idx]<=pos) {
//      *out_pos = checked_pos[checked_idx]; ++out_pos;
//      last_checked_pos = checked_pos[checked_idx++];
//    }
//    register const uint64_t checked_delta = last_checked_pos<=pos ? pos-last_checked_pos : IDX_T_MAX;
//    if (checked_delta==0) continue; // Repeated position
//    pending_queries[num_queries] = queries[i];
//    pending_queries[num_queries++].f_gap_length = UNSAFE_MIN(in_delta,checked_delta);
//    *out_pos = pos; ++out_pos; last_in_pos = pos;
//  }
//  while (checked_idx<num_checked_positions) { // Copy remaining ones
//    *out_pos = checked_pos[checked_idx++]; ++out_pos;
//  }
//  vector_update_used(out_positions,out_pos);
//  if (num_queries==0) return; // Nothing to do
//
//  // Calculate the effective number of mismatches/differences
//  register uint64_t eff_max_mismatches, eff_max_distance;
//  eff_max_mismatches = UNSAFE_MAX(search_params->max_mismatches,search_params->incomplete_strata);
//  eff_max_distance = UNSAFE_MAX(search_params->max_differences,search_params->incomplete_strata);
//  register const uint64_t max_align_diff = key_len-search_params->min_anchor_size;
//  if (eff_max_mismatches > max_align_diff) eff_max_mismatches=max_align_diff;
//  // eff_max_distance += matches->bad_quality_misms;
//  if (eff_max_distance > max_align_diff) eff_max_distance=max_align_diff;
//  if (search_params->max_differences == 0 || !use_levenshtein_distance) eff_max_distance=0;
//
//  // Allocate space for the positions
//  register const uint64_t max_misms = UNSAFE_MAX(eff_max_mismatches+matches->bad_quality_misms,eff_max_distance)+2;
//  matches_prepare(matches,max_misms+1);
//  vector_reserve(matches->rbuf_pos,vector_get_used(matches->rbuf_pos)+num_queries);
//  vector_reserve(matches->rbuf_mismatches,vector_get_used(matches->rbuf_mismatches)+max_misms*num_queries);
//
//  // Prepare peq vector
//  if (eff_max_distance > 0) {
//    vector_prepare(mpool->vbuf1,uint64_t,key_len*5*2);
//    search_params->peq = vector_get_mem(mpool->vbuf1);
//    fmi_prepare_pattern_eq(a,search_params->key,key_len,search_params->peq,
//        search_params->allowed_chars,search_params->mismatch_mask);
//  }
//
//  // Check candidates
//  fmi_matches_check_candidates(a,matches,search_params,store_hamming_matches,
//      eff_max_mismatches,eff_max_distance,pending_queries,num_queries,mpool);
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////







//
//
///*
// *
// */
//GEM_INLINE void fmi_matches_check_candidates(
//    const _FMI_* const a,matches* const matches,
//    fmi_search_parameters* const search_params,const bool store_hamming_matches,
//    const uint64_t max_mismatches,const uint64_t max_distance,
//    filter_query* queries,const uint64_t num_queries,vector_pool* const mpool) {
//  register ch_t *key = search_params->key;
//  register const uint64_t* const peq = search_params->peq;
//  register const uint64_t key_len = search_params->key_len;
//  register ch_t* mismatch_mask = search_params->mismatch_mask;
//  register const uint64_t max_indel_len = search_params->max_indel_len;
//  register const bool* const allowed_chars = search_params->allowed_chars;
//  register const uint64_t min_anchor_size = search_params->min_anchor_size;
//  // register const bool unique_mapping = search_params->unique_mapping; // TODO
//
//  START_STATS(SC_CHECK_CAND);
//
//  register uint64_t* cntrs = vector_get_mem(matches->rbuf_counts);
//  register pos_match_t* out = (pos_match_t*)vector_get_mem(matches->rbuf_pos)+vector_get_used(matches->rbuf_pos);
//  register mismatch* misms = (mismatch*)vector_get_mem(matches->rbuf_mismatches)+vector_get_used(matches->rbuf_mismatches);
//
//  register const uint64_t max_pattern_ext = max_indel_len+UNSAFE_MAX(max_distance,max_mismatches);
//  vector_prepare(mpool->fbuf2, ch_t, 2*max_pattern_ext+key_len); // We use fbuf2 to store the decoded read
//  register ch_t* decoded = vector_get_mem(mpool->fbuf2);
//
//  register uint64_t num_matches, i;
//  uint64_t match_mismatches;
//  idx_t match_position, match_end_position;
//  for (i = 0, num_matches=0; i < num_queries; ++i, ++queries) {
//    INC_COUNTER(GSC_CHECKED_MATCHES);
//
//    //
//    // Check Hamming distance
//    //
//    register const idx_t position = queries->end_key_pos;
//    register bool correct = true;
//    if (position+key_len < a->text_length) {
//      // Decode
//      START_STATS(SC_CHECK_DECODE);
//      register idx_t delta_position = position-UNSAFE_MIN(max_distance,queries->end_key_pos);
//      register idx_t delta = position-delta_position;
//      register idx_t delta_end_key = UNSAFE_MIN(position+key_len+max_distance,a->text_length);
//      register idx_t delta_key_len = delta_end_key-delta_position;
//      register idx_t delta_decoded_len = fmi_decode(a,delta_position,delta_key_len,decoded);
//      register bool exists_gap = false;
//      if (__builtin_expect(delta_decoded_len < key_len+delta,false)) { // There is a gap (Separator)
//        exists_gap = true;
//        register const idx_t gap_pos = delta_position+delta_decoded_len;
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
//          register const idx_t chunk_pos = position+queries->end;
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
//          register bool initial_trim = false, final_trim = false;
//          register idx_t initial_trim_len = 0, final_trim_len = 0;
//          register ch_t* offset_key = key;
//          register int64_t offset_key_len = key_len;
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
//          register uint64_t* const offset_peq = (uint64_t*)vector_get_mem(mpool->vbuf1)+(key_len*5);
//          fmi_prepare_pattern_eq(a,offset_key,offset_key_len,offset_peq,allowed_chars,mismatch_mask);
//          // Adjust the maximum differences allowed
//          register uint64_t max_anchor_distance = (offset_key_len*search_params->max_differences)/key_len;
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
//              register const uint64_t total_len = offset_key_len+(initial_trim?initial_trim_len:0);
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
//    register const idx_t indel_delta = queries->f_gap_length;
//    register const uint64_t indel_key_len = key_len+2*indel_delta;
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
//}
//
//
//
//
///*
// *
// */
//GEM_TOPLEVEL_INLINE bool fmi_matches_check_levenshtein_match(
//    const ch_t* const key,const uint64_t* const peq,const uint64_t key_length,
//    const ch_t* const decoded_text,const uint64_t text_length,
//    const uint64_t distance,const bool allow_indels,
//    mismatch* const misms,idx_t* match_position,idx_t* match_end_position,
//    uint64_t* match_distance,vector_pool* const mpool);
