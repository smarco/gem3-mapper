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
 * NSearch Query Helper
 */
typedef struct {
  fm_2interval_t fm_2interval;
  fm_2erank_elms_t lo_2erank_elms;
  fm_2erank_elms_t hi_2erank_elms;
  uint64_t max_steps;
} nsearch_query_t;

/*
 * Levenshtein bidirectional query
 */
void nsearch_levenshtein_scheduled_precompute_query(
    nsearch_schedule_t* const nsearch_schedule,
    const bool forward_search,
    nsearch_query_t* const nsearch_query) {
#ifndef NSEARCH_ENUMERATE
  fm_index_t* const fm_index = nsearch_schedule->search->archive->fm_index;
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
  fm_index_t* const fm_index = nsearch_schedule->search->archive->fm_index;
  if (forward_search) {
    fm_index_2query_forward_query(fm_index,fm_2interval,fm_2interval,char_enc);
  } else {
    fm_index_2query_backward_query(fm_index,fm_2interval,fm_2interval,char_enc);
  }
  // Return
  return fm_2interval->backward_hi - fm_2interval->backward_lo;
#endif
}
bool nsearch_levenshtein_scheduled_cutoff(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t num_candidates,
    nsearch_query_t* const nsearch_query,
    nsearch_query_t* const next_nsearch_query) {
  // Parameters
  search_parameters_t* const search_parameters = nsearch_schedule->search->search_parameters;
  const uint64_t ns_filtering_threshold = search_parameters->region_profile_model.ns_filtering_threshold;
  // Check number of candidates
  if (num_candidates <= ns_filtering_threshold) {
    if (nsearch_query->max_steps==0) {
      next_nsearch_query->max_steps = search_parameters->region_profile_model.max_steps;
      return false;
    } else {
      if (nsearch_query->max_steps == 1) return true;
      next_nsearch_query->max_steps = nsearch_query->max_steps - 1;
      return false;
    }
  } else {
    next_nsearch_query->max_steps = 0;
    return false;
  }
}
uint64_t nsearch_levenshtein_scheduled_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    nsearch_query_t* const nsearch_query,
    const uint64_t align_distance) {
  // Parameters
  fm_2interval_t* const fm_2interval = &nsearch_query->fm_2interval;
#ifdef NSEARCH_ENUMERATE
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,text_length);
  PROF_ADD_COUNTER(GP_NS_CANDIDATES_GENERATED,1);
  // Search Text
  nsearch_operation->text_position = text_length;
  nsearch_operation_state_print_global_text(stdout,nsearch_operation); // Print
  return 1;
#else
  // Parameters
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  filtering_candidates_t* const filtering_candidates = nsearch_schedule->search->filtering_candidates;
  search_parameters_t* const search_parameters = nsearch_schedule->search->search_parameters;
  pattern_t* const pattern = &nsearch_schedule->search->pattern;
  nsearch_operation->text_position = text_length;
  // Compute region boundaries
  const uint64_t max_error = nsearch_schedule->max_error;
  const uint64_t global_key_begin = nsearch_operation->global_key_begin;
  const uint64_t global_key_end = nsearch_operation->global_key_end;
  const uint64_t max_region_scope = text_length + max_error;
  uint64_t region_begin, region_end;
  if (forward_search) {
    region_begin = BOUNDED_SUBTRACTION(global_key_begin,max_error,0);
    region_end = global_key_begin + max_region_scope;
  } else {
    region_begin = BOUNDED_SUBTRACTION(global_key_end,max_region_scope,0);
    region_end = global_key_end + max_error;
  }
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,text_length);
  PROF_ADD_COUNTER(GP_NS_BRANCH_CANDIDATES_GENERATED,
      (fm_2interval->backward_hi-fm_2interval->backward_lo));
  // Add interval
  bool limited;
  filtering_candidates_add_region_interval(
      filtering_candidates,search_parameters,pattern,fm_2interval->backward_lo,
      fm_2interval->backward_hi,region_begin,region_end,align_distance,&limited);
  return fm_2interval->backward_hi-fm_2interval->backward_lo;
#endif
}
/*
 * Scheduled-Operation check distance
 */
bool nsearch_levenshtein_scheduled_align_distance_pass(
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
//
uint64_t nsearch_levenshtein_scheduled_search_next(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const current_nsearch_operation,
    nsearch_query_t* const nsearch_query);
uint64_t nsearch_levenshtein_scheduled_operation_step_query(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    nsearch_query_t* const nsearch_query);
//
uint64_t nsearch_levenshtein_scheduled_operation_step_next(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    nsearch_query_t* const nsearch_query,
    const uint64_t align_distance) {
  //  // Debug
  //  if (nsearch_operation_state_text_eq(nsearch_operation,"ACGGTAC",nsearch_schedule->mm_stack)) {
  //    nsearch_operation_state_print(stderr,nsearch_operation,nsearch_schedule->key);
  //  }
  if (pending_searches==0) {
    if (nsearch_levenshtein_scheduled_align_distance_pass(nsearch_operation,text_length,align_distance)) {
      // EOS
      return nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,text_length,nsearch_query,align_distance);
    } else {
      // Keep doing queries in within this operation
      return nsearch_levenshtein_scheduled_operation_step_query(
          nsearch_schedule,pending_searches,nsearch_operation,nsearch_query);
    }
  } else {
    if (nsearch_levenshtein_scheduled_align_distance_pass(nsearch_operation,text_length,align_distance)) {
      // Next search-operation
      return nsearch_levenshtein_scheduled_search_next(
          nsearch_schedule,pending_searches,nsearch_operation,nsearch_query);
    } else {
      // Keep doing queries in within this operation
      return nsearch_levenshtein_scheduled_operation_step_query(
          nsearch_schedule,pending_searches,nsearch_operation,nsearch_query);
    }
  }
}
uint64_t nsearch_levenshtein_scheduled_operation_step_query(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    nsearch_query_t* const nsearch_query) {
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
  const uint64_t max_error = nsearch_operation->max_global_error;
  // Consider adding epsilon-char
  if (pending_searches > 0) {
    align_distance =
        nsearch_levenshtein_state_get_global_align_distance(
            nsearch_state,key_chunk_length,text_position,max_error);
    if (nsearch_levenshtein_scheduled_align_distance_pass(nsearch_operation,text_position,align_distance)) {
      total_matches += nsearch_levenshtein_scheduled_search_next(
          nsearch_schedule,pending_searches,nsearch_operation,nsearch_query);
      NSEARCH_PROF_NODE(nsearch_schedule,total_matches); // PROFILE
      return total_matches;
    }
  }
  // Search all characters (expand node)
  nsearch_query_t next_nsearch_query;
  uint8_t char_enc;
  nsearch_levenshtein_scheduled_precompute_query(
      nsearch_schedule,forward_search,nsearch_query); // Precompute all ranks
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    nsearch_levenshtein_state_compute_chararacter_banded(
        nsearch_state,forward_search,key_chunk,key_chunk_length,
        text_position,char_enc,max_error,&min_val,&align_distance);
    if (min_val > max_error) { // Check max-error constraints
      PROF_INC_COUNTER(GP_NS_NODES_DP_PREVENT);
      continue;
    }
    // Query
    const uint64_t num_candidates = nsearch_levenshtein_scheduled_query(
        nsearch_schedule,nsearch_operation,forward_search,
        text_position,char_enc,nsearch_query,&next_nsearch_query);
    if (num_candidates==0) continue;
    const bool finish_search = nsearch_levenshtein_scheduled_cutoff(
        nsearch_schedule,num_candidates,nsearch_query,&next_nsearch_query);
    if (finish_search) {
      total_matches += nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,
          text_position+1,&next_nsearch_query,align_distance); // Terminate
    } else {
      total_matches += nsearch_levenshtein_scheduled_operation_step_next(
          nsearch_schedule,pending_searches,nsearch_operation,
          text_position+1,&next_nsearch_query,align_distance); // Next
    }
  }
  NSEARCH_PROF_NODE(nsearch_schedule,total_matches); // PROFILE
  return total_matches;
}
uint64_t nsearch_levenshtein_scheduled_operation_step_query_exact(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    nsearch_query_t* const nsearch_query) {
  // Parameters
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  const uint8_t* const key_chunk = nsearch_schedule->key + nsearch_operation->global_key_begin;
  const uint64_t key_chunk_length = nsearch_operation->global_key_end - nsearch_operation->global_key_begin;
  // Search all characters (expand node)
  uint64_t i;
  for (i=0;i<key_chunk_length;++i) {
    const uint64_t key_idx = forward_search ? i : (key_chunk_length-i-1);
    const uint8_t char_enc = key_chunk[key_idx];
    // Query
    const uint64_t num_candidates = nsearch_levenshtein_scheduled_query_exact(
        nsearch_schedule,nsearch_operation,forward_search,
        nsearch_operation->text_position,char_enc,nsearch_query);
    NSEARCH_PROF_NODE(nsearch_schedule,num_candidates); // PROFILE
    if (num_candidates==0) return 0;
    const bool finish_search = nsearch_levenshtein_scheduled_cutoff(
        nsearch_schedule,num_candidates,nsearch_query,nsearch_query);
    if (finish_search) {
      return nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,
          nsearch_operation->text_position,nsearch_query,0);
    }
  }
  // Keep searching
  if (pending_searches==0) {
    return nsearch_levenshtein_scheduled_terminate(
        nsearch_schedule,nsearch_operation,nsearch_operation->text_position+1,nsearch_query,0);
  } else {
    return nsearch_levenshtein_scheduled_search_next(
        nsearch_schedule,pending_searches,nsearch_operation,nsearch_query);
  }
}
/*
 * Scheduled Search (Chains scheduled operations)
 */
uint64_t nsearch_levenshtein_scheduled_search_next(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const current_nsearch_operation,
    nsearch_query_t* const nsearch_query) {
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
  return nsearch_levenshtein_scheduled_operation_step_query(nsearch_schedule,
      next_pending_searches,next_nsearch_operation,nsearch_query);
}
uint64_t nsearch_levenshtein_scheduled_search(nsearch_schedule_t* const nsearch_schedule) {
  // DEBUG
  // nsearch_schedule_print_pretty(stderr,nsearch_schedule);
  // Prepare a new search
  const uint64_t next_pending_searches = nsearch_schedule->num_pending_searches-1;
  nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
  next_nsearch_operation->text_position = 0;
  // Prepare search-query
  nsearch_query_t nsearch_query;
  nsearch_query.max_steps = 0;
#ifdef NSEARCH_ENUMERATE
  nsearch_query.fm_2interval.backward_lo = 0;
  nsearch_query.fm_2interval.backward_hi = 1;
  nsearch_query.fm_2interval.forward_lo = 0;
  nsearch_query.fm_2interval.forward_hi = 1;
#else
  fm_index_2query_init(nsearch_schedule->search->archive->fm_index,&nsearch_query.fm_2interval);
#endif
  // Search
  if (next_nsearch_operation->max_global_error == 0) {
    return nsearch_levenshtein_scheduled_operation_step_query_exact(nsearch_schedule,
        next_pending_searches,next_nsearch_operation,&nsearch_query);
  } else {
    return nsearch_levenshtein_scheduled_operation_step_query(nsearch_schedule,
        next_pending_searches,next_nsearch_operation,&nsearch_query);
  }
}
