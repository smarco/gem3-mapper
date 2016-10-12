/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_levenshtein_query.h"

/*
 * Setup
 */
void nsearch_query_init(
    nsearch_query_t* const nsearch_query,
    fm_index_t* const fm_index) {
  // Init interval
#ifdef NSEARCH_ENUMERATE
  nsearch_query->fm_2interval.backward_lo = 0;
  nsearch_query->fm_2interval.backward_hi = 1;
  nsearch_query->fm_2interval.forward_lo = 0;
  nsearch_query->fm_2interval.forward_hi = 1;
#else
  fm_index_2query_init(fm_index,&nsearch_query->fm_2interval);
#endif
  // Init query cut-offs
  nsearch_query->num_optimization_steps = 0;
  nsearch_query->num_eq_candidates_steps = 0;
  nsearch_query->prev_num_candidates = UINT64_MAX;
}

/*
 * Standard search query
 */
void nsearch_levenshtein_query(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t current_position,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out) {
#ifdef NSEARCH_ENUMERATE
  nsearch_schedule->pending_searches->text[current_position] = char_enc;
  *lo_out = 0; *hi_out = 1;
#else
  fm_index_t* const fm_index = nsearch_schedule->archive->fm_index;
  *lo_out = bwt_erank(fm_index->bwt,char_enc,lo_in);
  *hi_out = bwt_erank(fm_index->bwt,char_enc,hi_in);
#endif
}
/*
 * Scheduled search precomputed-query
 */
void nsearch_levenshtein_scheduled_precompute_query(
    nsearch_schedule_t* const nsearch_schedule,
    const bool forward_search,
    nsearch_query_t* const nsearch_query) {
#ifndef NSEARCH_ENUMERATE
  fm_index_t* const fm_index = nsearch_schedule->archive->fm_index;
  if (forward_search) {
    fm_index_2query_forward_precompute(
        fm_index,&nsearch_query->fm_2interval,
        &nsearch_query->lo_2erank_elms,&nsearch_query->hi_2erank_elms);
  } else {
    fm_index_2query_backward_precompute(
        fm_index,&nsearch_query->fm_2interval,
        &nsearch_query->lo_2erank_elms,&nsearch_query->hi_2erank_elms);
  }
#endif
}
/*
 * Scheduled search query
 */
uint64_t nsearch_levenshtein_scheduled_query_exact(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const int64_t current_position,
    const uint8_t char_enc,
    nsearch_query_t* const nsearch_query) {
  // Store character
  nsearch_operation->text[current_position] = char_enc;
  nsearch_operation->text_position = current_position+1;
  // Query character
  fm_2interval_t* const fm_2interval = &nsearch_query->fm_2interval;
#ifdef NSEARCH_ENUMERATE
  fm_2interval->backward_lo = 0;
  fm_2interval->backward_hi = 1;
  fm_2interval->forward_lo = 0;
  fm_2interval->forward_hi = 1;
  // Return
  return 1;
#else
  fm_index_t* const fm_index = nsearch_schedule->archive->fm_index;
  if (forward_search) {
    fm_index_2query_forward_query(fm_index,fm_2interval,fm_2interval,char_enc);
  } else {
    fm_index_2query_backward_query(fm_index,fm_2interval,fm_2interval,char_enc);
  }
  // Return
  return fm_2interval->backward_hi - fm_2interval->backward_lo;
#endif
}
uint64_t nsearch_levenshtein_scheduled_query(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const int64_t current_position,
    const uint8_t char_enc,
    nsearch_query_t* const nsearch_query,
    nsearch_query_t* const next_nsearch_query) {
  // Store character
  nsearch_operation->text[current_position] = char_enc;
  nsearch_operation->text_position = current_position+1;
  // Query character
  fm_2interval_t* const fm_2interval = &next_nsearch_query->fm_2interval;
#ifdef NSEARCH_ENUMERATE
  next_nsearch_query->fm_2interval_out.backward_lo = 0;
  next_nsearch_query->fm_2interval_out.backward_hi = 1;
  next_nsearch_query->fm_2interval_out.forward_lo = 0;
  next_nsearch_query->fm_2interval_out.forward_hi = 1;
  return 1;
#else
  if (forward_search) {
    fm_index_2query_precomputed_forward_query(
        &nsearch_query->lo_2erank_elms,&nsearch_query->hi_2erank_elms,
        &nsearch_query->fm_2interval,fm_2interval,char_enc);
  } else {
    fm_index_2query_precomputed_backward_query(
        &nsearch_query->lo_2erank_elms,&nsearch_query->hi_2erank_elms,
        &nsearch_query->fm_2interval,fm_2interval,char_enc);
  }
  // Return
  return fm_2interval->backward_hi - fm_2interval->backward_lo;
#endif
}
