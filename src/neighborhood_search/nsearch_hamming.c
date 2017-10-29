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

#include "align/alignment.h"
#include "neighborhood_search/nsearch_hamming.h"
#include "neighborhood_search/nsearch_partition.h"
#include "neighborhood_search/nsearch_schedule.h"
#include "neighborhood_search/nsearch_filtering.h"
#include "fm_index/fm_index_query.h"
#include "filtering/candidates/filtering_candidates_accessors.h"

/*
 * Query
 */
void nsearch_hamming_query(
    nsearch_schedule_t* const nsearch_schedule,
    const int64_t current_position,
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
void nsearch_hamming_directional_query(
    nsearch_schedule_t* const nsearch_schedule,
    const search_direction_t search_direction,
    const int64_t current_position,
    const uint8_t char_enc,
    fm_2interval_t* const fm_2interval_in,
    fm_2interval_t* const fm_2interval_out) {
#ifdef NSEARCH_ENUMERATE
  nsearch_schedule->pending_searches->text[current_position] = char_enc;
  fm_2interval_out->backward_lo = 0;
  fm_2interval_out->backward_hi = 1;
  fm_2interval_out->forward_lo = 0;
  fm_2interval_out->forward_hi = 1;
#else
  fm_index_t* const fm_index = nsearch_schedule->archive->fm_index;
  if (search_direction==direction_forward) {
    fm_index_2query_forward_query(fm_index,fm_2interval_in,fm_2interval_out,char_enc);
  } else {
    fm_index_2query_backward_query(fm_index,fm_2interval_in,fm_2interval_out,char_enc);
  }
#endif
}
void nsearch_hamming_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t lo,
    const uint64_t hi,
    const uint64_t align_distance) {
  // PROFILE
  PROF_ADD_COUNTER(GP_NS_SEARCH_DEPTH,nsearch_schedule->pattern->key_length);
  PROF_ADD_COUNTER(GP_NS_BRANCH_CANDIDATES_GENERATED,(hi-lo));
#ifdef NSEARCH_ENUMERATE
  const uint8_t* const text = nsearch_schedule->pending_searches->text;
  dna_buffer_print(stdout,text,nsearch_schedule->pattern->key_length,false);
  fprintf(stdout,"\n");
#else
  filtering_candidates_t* const filtering_candidates = nsearch_schedule->filtering_candidates;
  search_parameters_t* const search_parameters = nsearch_schedule->search_parameters;
  pattern_t* const pattern = nsearch_schedule->pattern;
  filtering_candidates_add_positions_from_interval(
      filtering_candidates,search_parameters,pattern,
      lo,hi,0,pattern->key_length,align_distance);
#endif
}
/*
 * Brute Force
 */
void nsearch_hamming_brute_force_step(
    nsearch_schedule_t* const nsearch_schedule,
    const int64_t current_position,
    const uint64_t current_error,
    const uint64_t lo,
    const uint64_t hi) {
  uint64_t next_lo, next_hi;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    const uint64_t next_error = current_error +
        (char_enc!=(nsearch_schedule->pattern->key[current_position]) ? 1 : 0);
    if (next_error <= nsearch_schedule->max_error) {
      // Query
      nsearch_hamming_query(nsearch_schedule,current_position,char_enc,lo,hi,&next_lo,&next_hi);
      if (next_lo >= next_hi) continue;
      // Next step
      const int64_t next_position = current_position - 1;
      if (next_position >= 0) {
        nsearch_hamming_brute_force_step(nsearch_schedule,next_position,next_error,next_lo,next_hi);
      } else {
        nsearch_hamming_terminate(nsearch_schedule,next_lo,next_hi,next_error);
      }
    }
  }
  NSEARCH_PROF_NODE(nsearch_schedule,1); // PROFILE
}
void nsearch_hamming_brute_force(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROF_START(GP_NS_GENERATION);
  // Init
  nsearch_schedule_init(search->nsearch_schedule,nsearch_model_hamming,
      search->max_search_error,search->archive,
      &search->pattern,&search->region_profile,
      search->search_parameters,search->filtering_candidates,
      matches);
#ifdef NSEARCH_ENUMERATE
  const uint64_t init_lo = 0;
  const uint64_t init_hi = 1;
#else
  const uint64_t init_lo = 0;
  const uint64_t init_hi = fm_index_get_length(search->nsearch_schedule->archive->fm_index);
#endif
  // Search
  nsearch_hamming_brute_force_step(search->nsearch_schedule,
      search->nsearch_schedule->pattern->key_length-1,0,init_lo,init_hi);
  // PROFILE
  // nsearch_schedule_print_profile(stderr,search->nsearch_schedule);
  PROF_ADD_COUNTER(GP_NS_NODES,search->nsearch_schedule->profile.ns_nodes);
  PROF_ADD_COUNTER(GP_NS_NODES_SUCCESS,search->nsearch_schedule->profile.ns_nodes_success);
  PROF_ADD_COUNTER(GP_NS_NODES_FAIL,search->nsearch_schedule->profile.ns_nodes_fail);
  PROF_STOP(GP_NS_GENERATION);
}
/*
 * Scheduled search
 */
void nsearch_hamming_scheduled_search_operation_next(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    fm_2interval_t* const fm_2interval,
    const uint64_t local_error,
    const uint64_t global_error,
    const int64_t current_position) {
  NSEARCH_PROF_NODE(nsearch_schedule,1); // PROFILE
  // Compute limits
  const search_direction_t search_direction = nsearch_operation->search_direction;
  int64_t next_position;
  bool limit_reached;
  if (search_direction==direction_forward) {
    next_position = current_position + 1;
    limit_reached = (next_position >= (int64_t)nsearch_operation->local_key_end);
  } else {
    next_position = current_position - 1;
    limit_reached = (next_position < (int64_t)nsearch_operation->local_key_begin);
  }
  // Search all characters
  fm_2interval_t next_fm_2interval;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute distance function & constraints
    const uint64_t next_local_error = local_error +
        ((char_enc!=nsearch_schedule->pattern->key[current_position]) ? 1 : 0);
    if (next_local_error+global_error > nsearch_operation->max_global_error) continue;
    if (limit_reached) {
      if (next_local_error < nsearch_operation->min_local_error) continue;
      if (next_local_error+global_error < nsearch_operation->min_global_error) continue;
    }
    // Query character
    nsearch_hamming_directional_query(nsearch_schedule,
        search_direction,current_position,char_enc,fm_2interval,&next_fm_2interval);
    if (next_fm_2interval.backward_lo >= next_fm_2interval.backward_hi) continue;
    // Keep searching
    if (!limit_reached) {
      nsearch_hamming_scheduled_search_operation_next(
          nsearch_schedule,pending_searches,nsearch_operation,
          &next_fm_2interval,next_local_error,global_error,next_position);
    } else {
      // Check pending searches
      if (pending_searches==0) {
        nsearch_hamming_terminate(
            nsearch_schedule,next_fm_2interval.backward_lo,
            next_fm_2interval.backward_hi,next_local_error+global_error);
      } else {
        // Search next chunk
        const uint64_t next_pending_searches = pending_searches-1;
        nsearch_operation_t* const next_nsearch_operation =
            nsearch_schedule->pending_searches + next_pending_searches;
        const int64_t next_begin_position =
            (next_nsearch_operation->search_direction==direction_forward) ?
                next_nsearch_operation->local_key_begin :
                next_nsearch_operation->local_key_end - 1;
        nsearch_hamming_scheduled_search_operation_next(
            nsearch_schedule,next_pending_searches,next_nsearch_operation,
            &next_fm_2interval,0,next_local_error+global_error,next_begin_position);
      }
    }
  }
}
void nsearch_hamming_scheduled_search(nsearch_schedule_t* const nsearch_schedule) {
  // Parameters
  const uint64_t pending_searches = nsearch_schedule->num_pending_searches - 1;
  nsearch_operation_t* const nsearch_operation = nsearch_schedule->pending_searches + pending_searches;
  // Init query
  fm_2interval_t fm_2interval;
#ifndef NSEARCH_ENUMERATE
  fm_index_2query_init(nsearch_schedule->archive->fm_index,&fm_2interval);
#endif
  // Launch search
  // nsearch_schedule_print_pretty(stderr,nsearch_schedule); // DEBUG
  const int64_t begin_position =
      (nsearch_operation->search_direction==direction_forward) ?
          nsearch_operation->local_key_begin :
          nsearch_operation->local_key_end - 1;
  nsearch_hamming_scheduled_search_operation_next(nsearch_schedule,
      pending_searches,nsearch_operation,&fm_2interval,0,0,begin_position);
}
/*
 * Hamming Neighborhood Search
 */
void nsearch_hamming(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Search
  nsearch_schedule_init(
      search->nsearch_schedule,nsearch_model_hamming,
      search->max_search_error,search->archive,
      &search->pattern,&search->region_profile,
      search->search_parameters,search->filtering_candidates,
      matches);
  nsearch_schedule_search(search->nsearch_schedule);
  nsearch_filtering(search->nsearch_schedule);
  // PROFILE
#ifdef GEM_PROFILE
  // nsearch_schedule_print_profile(stderr,search->nsearch_schedule);
  if (search->filtering_candidates != NULL) {
    PROF_ADD_COUNTER(GP_NS_SEARCH_CANDIDATES_GENERATED,
        filtering_candidates_get_num_positions(search->filtering_candidates));
  }
#endif
}
/*
 * Neighborhood Search (Preconditioned by region profile)
 */
void nsearch_hamming_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Search
  nsearch_schedule_init(
      search->nsearch_schedule,nsearch_model_hamming,
      search->max_search_error,search->archive,
      &search->pattern,&search->region_profile,
      search->search_parameters,search->filtering_candidates,
      matches);
  nsearch_schedule_search_preconditioned(search->nsearch_schedule);
  nsearch_filtering(search->nsearch_schedule);
  // PROFILE
#ifdef GEM_PROFILE
  // nsearch_schedule_print_profile(stderr,&nsearch_schedule);
  if (search->filtering_candidates != NULL) {
    PROF_ADD_COUNTER(GP_NS_SEARCH_CANDIDATES_GENERATED,
        filtering_candidates_get_num_positions(search->filtering_candidates));
  }
#endif
}
