/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "neighborhood_search/nsearch_levenshtein_scheduled.h"
#include "neighborhood_search/nsearch_levenshtein.h"
#include "neighborhood_search/nsearch_levenshtein_state.h"
#include "neighborhood_search/nsearch_levenshtein_query.h"
#include "neighborhood_search/nsearch_levenshtein_control.h"
#include "neighborhood_search/nsearch_partition.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/candidates/filtering_candidates_align.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_verify.h"

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
uint64_t nsearch_levenshtein_scheduled_operation_step_next(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t text_length,
    nsearch_query_t* const nsearch_query,
    const uint64_t align_distance) {
  //  // Debug
  //  if (nsearch_operation_state_text_eq(nsearch_operation,"ACGGTAC",nsearch_schedule->mm_allocator)) {
  //    nsearch_operation_state_print(stderr,nsearch_operation,nsearch_schedule->key);
  //  }
  if (pending_searches==0) {
    if (nsearch_levenshtein_scheduled_align_distance_pass(nsearch_operation,text_length,align_distance)) {
      // EOS
      return nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,
          text_length,nsearch_query,align_distance,true);
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
  const uint8_t* const key_chunk = nsearch_schedule->pattern->key + nsearch_operation->global_key_begin;
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
        nsearch_state,forward_search,key_chunk,
        key_chunk_length,text_position,char_enc,max_error,
        &min_val,&align_distance,nsearch_schedule->mm_allocator);
    if (min_val > max_error) { // Check max-error constraints
      PROF_INC_COUNTER(GP_NS_NODES_DP_PREVENT);
      continue;
    }
    // Query
    const uint64_t num_candidates = nsearch_levenshtein_scheduled_query(
        nsearch_schedule,nsearch_operation,forward_search,
        text_position,char_enc,nsearch_query,&next_nsearch_query);
    if (num_candidates==0) continue;
    const bool finish_search = nsearch_levenshtein_candidates_cutoff(
        nsearch_schedule,num_candidates,nsearch_query,&next_nsearch_query);
    if (finish_search) {
      total_matches += nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,
          text_position+1,&next_nsearch_query,align_distance,false); // Terminate
    } else {
      total_matches += nsearch_levenshtein_scheduled_operation_step_next(
          nsearch_schedule,pending_searches,nsearch_operation,
          text_position+1,&next_nsearch_query,align_distance); // Next
    }
    if (nsearch_schedule->quick_abandon) break; // Quick abandon
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
  const uint8_t* const key_chunk = nsearch_schedule->pattern->key + nsearch_operation->global_key_begin;
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
    const bool finish_search = nsearch_levenshtein_candidates_cutoff(
        nsearch_schedule,num_candidates,nsearch_query,nsearch_query);
    if (finish_search) {
      return nsearch_levenshtein_scheduled_terminate(
          nsearch_schedule,nsearch_operation,
          nsearch_operation->text_position,nsearch_query,0,false);
    }
  }
  // Keep searching
  if (pending_searches==0) {
    return nsearch_levenshtein_scheduled_terminate(
        nsearch_schedule,nsearch_operation,
        nsearch_operation->text_position+1,nsearch_query,0,true);
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
      nsearch_schedule->pattern->key,nsearch_schedule->pattern->key_length,
      reverse_sequence,nsearch_schedule->mm_allocator);
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
  nsearch_query_init(&nsearch_query,nsearch_schedule->archive->fm_index);
  // Search
  if (next_nsearch_operation->max_global_error == 0) {
    return nsearch_levenshtein_scheduled_operation_step_query_exact(nsearch_schedule,
        next_pending_searches,next_nsearch_operation,&nsearch_query);
  } else {
    return nsearch_levenshtein_scheduled_operation_step_query(nsearch_schedule,
        next_pending_searches,next_nsearch_operation,&nsearch_query);
  }
}
