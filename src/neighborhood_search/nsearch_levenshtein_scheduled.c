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
    nsearch_operation_t* const nsearch_operation) {
  // Parameters
  nsearch_levenshtein_state_t* const nsearch_state = &nsearch_operation->nsearch_state;
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  const uint64_t max_error = nsearch_operation->max_global_error;
  const uint8_t* const key_chunk = nsearch_schedule->key + nsearch_operation->local_key_begin;
  const uint64_t key_chunk_length = nsearch_operation->local_key_end - nsearch_operation->local_key_begin;
  const uint64_t text_position = nsearch_operation->text_position;
  // Search all characters (expand node)
  uint64_t total_matches = 0, num_matches;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    uint64_t min_val, align_distance;
    const bool active_column =
        nsearch_levenshtein_state_compute_chararacter(
            nsearch_state,forward_search,key_chunk,key_chunk_length,
            text_position,char_enc,max_error,&min_val,&align_distance);
    if (active_column && nsearch_state->first_active_column==0) {
      nsearch_state->first_active_column = text_position+1;
    }
    // Check max-error constraints
    if (min_val > max_error) continue; // || text_position < nsearch_schedule->max_text_length
    // Search character
    nsearch_operation->text[text_position] = char_enc;
    nsearch_operation->text_position = text_position+1;
    ++(nsearch_schedule->ns_nodes); // PROFILING
    // DEBUG
    // nsearch_levenshtein_print_search_trace(stderr,nsearch_schedule,pending_searches); // DEBUG
    // nsearch_levenshtein_state_print_text(stderr,nsearch_state,forward_search); // DEBUG
    // nsearch_levenshtein_print_pair_key_text(stderr,nsearch_schedule,nsearch_operation); // DEBUG
    // Check Local-Search finished
    if (pending_searches==0) {
      if (align_distance <= max_error) {
        // Print Search Text
        nsearch_operation_state_print_global_text(stdout,nsearch_operation,forward_search);
        ++nsearch_schedule->ns_nodes_success; // PROFILE
        ++total_matches; // Account the current match
      }
      // Keep searching
      total_matches += nsearch_levenshtein_scheduled_search_operation(
          nsearch_schedule,pending_searches,nsearch_operation);
    } else {
      // Keep searching // TODO Cleverer ways of finding out that num_matches==0
      num_matches = nsearch_levenshtein_scheduled_search_operation(
          nsearch_schedule,pending_searches,nsearch_operation);
      if (num_matches==0 && nsearch_state->first_active_column!=0) {
        // Propagate last row active cells & perform next operation
        total_matches += nsearch_levenshtein_scheduled_search(
            nsearch_schedule,pending_searches,nsearch_operation);
      } else {
        total_matches += num_matches;
      }
    }
  }
  return total_matches;
}
uint64_t nsearch_levenshtein_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const current_nsearch_operation) {
  /*
   * Initialize a new search
   */
  if (current_nsearch_operation==NULL) {
    const uint64_t next_pending_searches = pending_searches-1;
    nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
    nsearch_levenshtein_state_t* const next_nsearch_state = &next_nsearch_operation->nsearch_state;
    // Prepare DP
    next_nsearch_operation->global_text_length = 0;
    next_nsearch_operation->text_position = 0;
    nsearch_levenshtein_state_prepare_full_neighbourhood(next_nsearch_state);
    return nsearch_levenshtein_scheduled_search_operation(
        nsearch_schedule,next_pending_searches,next_nsearch_operation);
  }
//  /*
//   * Continue searching the next chunk (Prepare DP forward/backward)
//   */
//  const uint64_t next_pending_searches = pending_searches-1;
//  nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
//  nsearch_levenshtein_state_t* const current_nsearch_state = &current_nsearch_operation->nsearch_state;
//  nsearch_levenshtein_state_t* const next_nsearch_state = &next_nsearch_operation->nsearch_state;
//  // Prepare DP
//  const bool forward_search = (next_nsearch_operation->search_direction == direction_forward);
//  if (current_nsearch_operation->search_direction == next_nsearch_operation->search_direction) {
//    // Prepare context (Isolate active columns)
//    nsearch_operation_chained_prepare_sequence(current_nsearch_state,next_nsearch_state);
//    nsearch_operation_chained_prepare_context(current_nsearch_state,next_nsearch_state);
//    // Chain search (Same direction)
//    const uint64_t current_min_error = current_nsearch_operation->min_global_error;
//    nsearch_levenshtein_state_prepare_chained(current_nsearch_state,next_nsearch_state,current_min_error);
//    // Search context (between one tile and the next)
//    const uint64_t next_max_error = next_nsearch_operation->max_global_error;
//    const uint8_t* const key_chunk = nsearch_schedule->key + next_nsearch_state->key_begin;
//    const uint64_t key_chunk_length = next_nsearch_state->key_end - next_nsearch_state->key_begin;
//    nsearch_levenshtein_state_compute_text(
//        next_nsearch_state,forward_search,key_chunk,key_chunk_length,
//        next_nsearch_state->text,next_nsearch_state->context_length,next_max_error);
//    // nsearch_levenshtein_state_print(stderr,next_nsearch_state,forward_search,nsearch_schedule->key); // DEBUG
//  } else {
//    nsearch_levenshtein_state_print(stderr,current_nsearch_state,
//        current_nsearch_operation->search_direction == direction_forward,nsearch_schedule->key); // DEBUG
//    // Reverse DP (Reverse sequence, recompute DP & isolate active columns)
//    const uint64_t current_max_error = current_nsearch_operation->max_global_error;
//    nsearch_operation_chained_prepare_reverse_sequence(current_nsearch_state,next_nsearch_state);
//    nsearch_levenshtein_state_prepare_full_neighbourhood(next_nsearch_state);
//    const uint64_t global_text_length = next_nsearch_state->global_text_length;
//    const uint64_t global_key_begin = current_nsearch_state->key_end - global_text_length;
//    const uint8_t* const key_chunk = nsearch_schedule->key + next_nsearch_state->key_begin;
//    const uint64_t current_key_chunk_length = next_nsearch_state->key_end - next_nsearch_state->key_begin;
//    nsearch_levenshtein_state_compute_text(
//        next_nsearch_state,forward_search,nsearch_schedule->key,
//        next_nsearch_state->global_text,next_nsearch_state->global_text_length,current_max_error);
//    nsearch_levenshtein_state_print(stderr,next_nsearch_state,forward_search,nsearch_schedule->key); // DEBUG
//    // Chain search (Reversed direction)
//    const uint64_t current_min_error = current_nsearch_operation->min_global_error;
//    nsearch_levenshtein_state_prepare_chained(next_nsearch_state,next_nsearch_state,current_min_error);
//    // Search context (between one tile and the next)
//    const uint64_t next_max_error = next_nsearch_operation->max_global_error;
//    nsearch_levenshtein_state_compute_text(
//        next_nsearch_state,forward_search,nsearch_schedule->key,
//        next_nsearch_state->text,next_nsearch_state->context_length,next_max_error);
//    nsearch_levenshtein_state_print(stderr,next_nsearch_state,forward_search,nsearch_schedule->key); // DEBUG
//  }
//  /*
//   * Keep searching
//   */
//  return nsearch_levenshtein_scheduled_search_operation(
//      nsearch_schedule,next_pending_searches,next_nsearch_operation);
  return 0; // TODO
}
