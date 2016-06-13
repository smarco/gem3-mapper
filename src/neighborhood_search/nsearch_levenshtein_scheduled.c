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
 * Levenshtein Scheduled Search
 */
uint64_t nsearch_levenshtein_scheduled_search_operation(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation);
uint64_t nsearch_levenshtein_scheduled_search_operation_next(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    const uint64_t align_distance) {
  // Parameters
  const uint64_t max_error = nsearch_operation->max_local_error;
  const uint64_t min_error = nsearch_operation->min_local_error;
  // Check Local-Search finished
  uint64_t num_matches = 0;
  if (pending_searches==0) {
    if (min_error <= align_distance && align_distance <= max_error) {
      // Print Search Text
      nsearch_operation->text_position = text_length;
      nsearch_operation_state_print_global_text(stdout,nsearch_operation); // PRINT TEXT
//      nsearch_levenshtein_print_trace(stderr,nsearch_schedule); // PRINT TRACE (to find duplicates)
//      if (nsearch_operation_state_global_text_cmp(
//          nsearch_operation,"AGGTTAC",nsearch_schedule->mm_stack)==0) {
//        nsearch_levenshtein_print_status(stderr,nsearch_schedule,pending_searches);// DEBUG
//        printf("----\n");
//      }
      ++nsearch_schedule->ns_nodes_success; // PROFILE
      ++num_matches; // Account the current match
    } else {
      // Keep searching
      num_matches += nsearch_levenshtein_scheduled_search_operation(
              nsearch_schedule,pending_searches,nsearch_operation);
    }
    return num_matches;
  } else {
    // Step into next search-operation
    if (min_error <= align_distance && align_distance <= max_error) {
      // Propagate last row active cells & perform next operation
      num_matches += nsearch_levenshtein_scheduled_search(
          nsearch_schedule,pending_searches,nsearch_operation);
    } else {
      // Keep searching
      num_matches += nsearch_levenshtein_scheduled_search_operation(
            nsearch_schedule,pending_searches,nsearch_operation);
    }
  }
  return num_matches;
}
uint64_t nsearch_levenshtein_scheduled_search_operation(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation) {
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
  const uint64_t max_error = nsearch_operation->max_local_error;
  const uint64_t min_error = nsearch_operation->min_local_error;
  // Consider adding epsilon-char
  align_distance = nsearch_levenshtein_state_get_align_distance(nsearch_state,key_chunk_length,0);
  if (pending_searches>0 && min_error <= align_distance && align_distance <= max_error) {
    total_matches += nsearch_levenshtein_scheduled_search(
        nsearch_schedule,pending_searches,nsearch_operation);
    return total_matches;
  }
  // Search all characters (expand node)
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    nsearch_levenshtein_state_compute_chararacter(
        nsearch_state,forward_search,key_chunk,key_chunk_length,
        text_position,char_enc,max_error,&min_val,&align_distance);
    // Check max-error constraints
    if (min_val > max_error) continue;
    // Search character
    nsearch_operation->text[text_position] = char_enc;
    nsearch_operation->text_position = text_position+1;
    ++(nsearch_schedule->ns_nodes); // PROFILING
    // DEBUG
    // nsearch_operation_state_print(stderr,nsearch_operation,nsearch_schedule->key); // DEBUG
    // nsearch_levenshtein_print_status(stderr,nsearch_schedule,pending_searches); // DEBUG
    // nsearch_levenshtein_state_print_text(stderr,nsearch_state,forward_search); // DEBUG
    total_matches += nsearch_levenshtein_scheduled_search_operation_next(
        nsearch_schedule,pending_searches,nsearch_operation,text_position+1,align_distance);
  }
  return total_matches;
}
uint64_t nsearch_levenshtein_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const current_nsearch_operation) {
  // Check current operation
  if (current_nsearch_operation==NULL) {
    // nsearch_schedule_print_pretty(stderr,nsearch_schedule); // DEBUG
    const uint64_t next_pending_searches = pending_searches-1;
    nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
    nsearch_levenshtein_state_t* const next_nsearch_state = &next_nsearch_operation->nsearch_state;
    const bool next_operation_towards_end =
        (next_nsearch_operation->search_direction==direction_forward) ?
            (next_nsearch_operation->global_key_begin == 0) :
            (next_nsearch_operation->global_key_end == nsearch_schedule->key_length);
    // Prepare a new search
    next_nsearch_operation->global_text_length = 0;
    next_nsearch_operation->text_position = 0;
    if (next_operation_towards_end) {
      nsearch_levenshtein_state_prepare_supercondensed_neighbourhood(next_nsearch_state);
    } else {
      nsearch_levenshtein_state_prepare_full_neighbourhood(next_nsearch_state);
    }
    return nsearch_levenshtein_scheduled_search_operation(
        nsearch_schedule,next_pending_searches,next_nsearch_operation);
  } else {
    // Continue searching the next chunk (Prepare DP forward/backward)
    const uint64_t next_pending_searches = pending_searches-1;
    nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
    // Prepare Search
    if (current_nsearch_operation->search_direction == next_nsearch_operation->search_direction) {
      nsearch_operation_chained_prepare_forward(
          current_nsearch_operation,next_nsearch_operation,
          nsearch_schedule->key,nsearch_schedule->key_length);
    } else {
      // Compute reverse
      nsearch_operation_compute_reverse(
          current_nsearch_operation,nsearch_schedule->nsearch_operation_aux,
          nsearch_schedule->key,nsearch_schedule->key_length);
      // Chain forward
      nsearch_operation_chained_prepare_forward(
          nsearch_schedule->nsearch_operation_aux,next_nsearch_operation,
          nsearch_schedule->key,nsearch_schedule->key_length);
    }
    // Keep searching
    return nsearch_levenshtein_scheduled_search_operation(
        nsearch_schedule,next_pending_searches,next_nsearch_operation);
  }
}
