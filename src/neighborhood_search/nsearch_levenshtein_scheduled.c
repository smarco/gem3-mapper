/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_levenshtein_scheduled.h"
#include "neighborhood_search/nsearch_levenshtein.h"
#include "neighborhood_search/nsearch_levenshtein_state.h"
#include "neighborhood_search/nsearch_partition.h"

/*
 * Levenshtein bidirectional query
 */
void nsearch_levenshtein_scheduled_directional_query(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const int64_t current_position,
    fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out,
    const uint8_t char_enc) {
  NSEARCH_PROF_ADD_NODE(nsearch_schedule); // PROFILE
  // Store character
  nsearch_operation->text[current_position] = char_enc;
  nsearch_operation->text_position = current_position+1;
  // Query character
#ifdef NSEARCH_ENUMERATE
  fm_2interval_out->backward_lo = 0; fm_2interval_out->backward_hi = 1;
  fm_2interval_out->forward_lo = 0; fm_2interval_out->forward_hi = 1;
#else
  if (forward_search) {
    fm_index_2query_precomputed_forward_query(
        lo_2erank_elms,hi_2erank_elms,fm_2interval_in,fm_2interval_out,char_enc);
  } else {
    fm_index_2query_precomputed_backward_query(
        lo_2erank_elms,hi_2erank_elms,fm_2interval_in,fm_2interval_out,char_enc);
  }
#endif
}
void nsearch_levenshtein_scheduled_directional_query_exact(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const int64_t current_position,
    fm_2interval_t* const fm_2interval,
    const uint8_t char_enc) {
  NSEARCH_PROF_ADD_NODE(nsearch_schedule); // PROFILE
  // Store character
  nsearch_operation->text[current_position] = char_enc;
  nsearch_operation->text_position = current_position+1;
  // Query character
#ifdef NSEARCH_ENUMERATE
  fm_2interval->backward_lo = 0; fm_2interval->backward_hi = 1;
  fm_2interval->forward_lo = 0; fm_2interval->forward_hi = 1;
#else
  fm_index_t* const fm_index = nsearch_schedule->search->archive->fm_index;
  if (forward_search) {
    fm_index_2query_forward_query(fm_index,fm_2interval,fm_2interval,char_enc);
  } else {
    fm_index_2query_backward_query(fm_index,fm_2interval,fm_2interval,char_enc);
  }
#endif
}
uint64_t nsearch_levenshtein_scheduled_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    fm_2interval_t* const fm_2interval,
    const uint64_t align_distance) {
  NSEARCH_PROF_ADD_SOLUTION(nsearch_schedule); // PROFILE
  // Print Search Text
  nsearch_operation->text_position = text_length;
#ifdef NSEARCH_ENUMERATE
  nsearch_operation_state_print_global_text(stdout,nsearch_operation); // PRINT TEXT
  return 1;
#else
  filtering_candidates_t* const filtering_candidates = nsearch_schedule->search->filtering_candidates;
  search_parameters_t* const search_parameters = nsearch_schedule->search->search_parameters;
  pattern_t* const pattern = &nsearch_schedule->search->pattern;
  const uint64_t region_begin = nsearch_operation->global_key_begin;
  const uint64_t region_end = nsearch_operation->global_key_end;
  filtering_candidates_add_region_interval(filtering_candidates,
      search_parameters,pattern,fm_2interval->backward_lo,
      fm_2interval->backward_hi,region_begin,region_end,align_distance);
  return fm_2interval->backward_hi-fm_2interval->backward_lo;
#endif
}
/*
 * Scheduled-Operation check distance
 */
bool nsearch_levenshtein_scheduled_search_distance_pass(
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    const uint64_t global_align_distance) {
  // Check global error
  const uint64_t max_error = nsearch_operation->max_global_error;
  const uint64_t min_global_error = nsearch_operation->min_global_error;
  if (min_global_error > global_align_distance || global_align_distance > max_error) return false;
  // Check local error
  const uint64_t min_local_error = nsearch_operation->min_local_error;
  const uint64_t global_key_length = nsearch_operation->global_key_end - nsearch_operation->global_key_begin;
  const uint64_t local_key_length = nsearch_operation->local_key_end - nsearch_operation->local_key_begin;
  const uint64_t local_align_distance =
      nsearch_levenshtein_state_get_local_align_distance(
          &nsearch_operation->nsearch_state,local_key_length,
          global_key_length,text_length,max_error);
  return (local_align_distance >= min_local_error);
}
/*
 * Scheduled-Operation Search (Performs each operational search)
 */
uint64_t nsearch_levenshtein_scheduled_search_operation_next(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    fm_2interval_t* const fm_2interval,
    const uint64_t align_distance) {
  //  // Debug
  //  if (nsearch_operation_state_text_eq(nsearch_operation,"ACGGTAC",nsearch_schedule->mm_stack)) {
  //    nsearch_operation_state_print(stderr,nsearch_operation,nsearch_schedule->key);
  //    printf("Here");
  //  }
  if (pending_searches==0) {
    if (nsearch_levenshtein_scheduled_search_distance_pass(nsearch_operation,text_length,align_distance)) {
      // EOS
      return nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,text_length,fm_2interval,align_distance);
    } else {
      // Keep doing queries in within this operation
      return nsearch_levenshtein_scheduled_search_operation_query(
          nsearch_schedule,pending_searches,nsearch_operation,fm_2interval);
    }
  } else {
    if (nsearch_levenshtein_scheduled_search_distance_pass(nsearch_operation,text_length,align_distance)) {
      // Next search-operation
      return nsearch_levenshtein_scheduled_search(
          nsearch_schedule,pending_searches,nsearch_operation,fm_2interval);
    } else {
      // Keep doing queries in within this operation
      return nsearch_levenshtein_scheduled_search_operation_query(
          nsearch_schedule,pending_searches,nsearch_operation,fm_2interval);
    }
  }
}
uint64_t nsearch_levenshtein_scheduled_search_operation_query(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    fm_2interval_t* const fm_2interval) {
  // Parameters Current Operation
  nsearch_levenshtein_state_t* const nsearch_state = &nsearch_operation->nsearch_state;
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  const uint8_t* const key_chunk = nsearch_schedule->key + nsearch_operation->global_key_begin;
  const uint64_t key_chunk_length = nsearch_operation->global_key_end - nsearch_operation->global_key_begin;
  // Parameters Current Text
  const uint64_t text_position = nsearch_operation->text_position;
  uint64_t min_val, align_distance;
  uint64_t total_matches = 0;
  // Parameters error
  search_parameters_t* const search_parameters = nsearch_schedule->search->search_parameters;
  const uint64_t ns_filtering_threshold = search_parameters->region_profile_model.ns_filtering_threshold;
  const uint64_t max_error = nsearch_operation->max_global_error;
  // Consider adding epsilon-char
  if (pending_searches > 0) {
    align_distance =
        nsearch_levenshtein_state_get_global_align_distance(
            nsearch_state,key_chunk_length,text_position,max_error);
    if (nsearch_levenshtein_scheduled_search_distance_pass(nsearch_operation,text_position,align_distance)) {
      total_matches += nsearch_levenshtein_scheduled_search(
          nsearch_schedule,pending_searches,nsearch_operation,fm_2interval);
      return total_matches;
    }
  }
  // Precompute all ranks
  fm_2erank_elms_t lo_2erank_elms, hi_2erank_elms;
#ifndef NSEARCH_ENUMERATE
  fm_index_t* const fm_index = nsearch_schedule->search->archive->fm_index;
  if (forward_search) {
    fm_index_2query_forward_precompute(fm_index,fm_2interval,&lo_2erank_elms,&hi_2erank_elms);
  } else {
    fm_index_2query_backward_precompute(fm_index,fm_2interval,&lo_2erank_elms,&hi_2erank_elms);
  }
#endif
  // Search all characters (expand node)
  fm_2interval_t next_fm_2interval;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    nsearch_levenshtein_state_compute_chararacter_banded(
        nsearch_state,forward_search,key_chunk,key_chunk_length,
        text_position,char_enc,max_error,&min_val,&align_distance);
    if (min_val > max_error) continue; // Check max-error constraints
    // Query
    nsearch_levenshtein_scheduled_directional_query(
        nsearch_schedule,nsearch_operation,forward_search,
        text_position,&lo_2erank_elms,&hi_2erank_elms,
        fm_2interval,&next_fm_2interval,char_enc);
    const uint64_t num_candidates = next_fm_2interval.backward_hi - next_fm_2interval.backward_lo;
    if (num_candidates <= ns_filtering_threshold) {
      if (num_candidates==0) continue;
      total_matches += nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,
          text_position+1,&next_fm_2interval,align_distance);
    } else {
      // Next
      total_matches += nsearch_levenshtein_scheduled_search_operation_next(
          nsearch_schedule,pending_searches,nsearch_operation,
          text_position+1,&next_fm_2interval,align_distance);
    }
  }
  return total_matches;
}
uint64_t nsearch_levenshtein_scheduled_search_operation_query_exact(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    fm_2interval_t* const fm_2interval) {
  // Parameters
  search_parameters_t* const search_parameters = nsearch_schedule->search->search_parameters;
  const uint64_t ns_filtering_threshold = search_parameters->region_profile_model.ns_filtering_threshold;
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  const uint8_t* const key_chunk = nsearch_schedule->key + nsearch_operation->global_key_begin;
  const uint64_t key_chunk_length = nsearch_operation->global_key_end - nsearch_operation->global_key_begin;
  // Search all characters (expand node)
  fm_2interval_t next_fm_2interval = *fm_2interval;
  uint64_t i;
  for (i=0;i<key_chunk_length;++i) {
    const uint64_t key_idx = forward_search ? i : (key_chunk_length-i-1);
    const uint8_t char_enc = key_chunk[key_idx];
    // Query
    nsearch_levenshtein_scheduled_directional_query_exact(
        nsearch_schedule,nsearch_operation,forward_search,
        nsearch_operation->text_position,&next_fm_2interval,char_enc);
    const uint64_t num_candidates = next_fm_2interval.backward_hi - next_fm_2interval.backward_lo;
    if (num_candidates <= ns_filtering_threshold) {
      if (num_candidates==0) return 0;
      return nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,
          nsearch_operation->text_position+1,&next_fm_2interval,0);
    }
  }
  // Keep searching
  if (pending_searches==0) {
    return nsearch_levenshtein_scheduled_terminate(
        nsearch_schedule,nsearch_operation,
        nsearch_operation->text_position+1,&next_fm_2interval,0);
  } else {
    return nsearch_levenshtein_scheduled_search(
        nsearch_schedule,pending_searches,nsearch_operation,&next_fm_2interval);
  }
}
/*
 * Scheduled Search (Chains scheduled operations)
 */
uint64_t nsearch_levenshtein_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const current_nsearch_operation,
    fm_2interval_t* const fm_2interval) {
  // Check current operation
  if (current_nsearch_operation==NULL) {
    // Prepare a new search
    // nsearch_schedule_print_pretty(stderr,nsearch_schedule); // DEBUG
    const uint64_t next_pending_searches = pending_searches-1;
    nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
    next_nsearch_operation->text_position = 0;
    // Search
    fm_2interval_t init_fm_2interval;
#ifdef NSEARCH_ENUMERATE
    init_fm_2interval.backward_lo = 0; init_fm_2interval.backward_hi = 1;
    init_fm_2interval.forward_lo = 0;  init_fm_2interval.forward_hi = 1;
#else
    fm_index_2query_init(nsearch_schedule->search->archive->fm_index,&init_fm_2interval);
#endif
    if (next_nsearch_operation->max_global_error == 0) {
      return nsearch_levenshtein_scheduled_search_operation_query_exact(nsearch_schedule,
          next_pending_searches,next_nsearch_operation,&init_fm_2interval);
    } else {
      return nsearch_levenshtein_scheduled_search_operation_query(nsearch_schedule,
          next_pending_searches,next_nsearch_operation,&init_fm_2interval);
    }
  } else {
    // Continue searching the next chunk (Prepare DP forward/backward)
    const uint64_t next_pending_searches = pending_searches-1;
    nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
    // Prepare Search-Operation
    const bool reverse_sequence =
        (current_nsearch_operation->search_direction != next_nsearch_operation->search_direction);
    const bool supercondensed = nsearch_operation_chained_prepare(
        current_nsearch_operation,next_nsearch_operation,
        nsearch_schedule->key,nsearch_schedule->key_length,reverse_sequence);
    if (!supercondensed) return 0;
    // Keep searching
    return nsearch_levenshtein_scheduled_search_operation_query(nsearch_schedule,
        next_pending_searches,next_nsearch_operation,fm_2interval);
  }
}
