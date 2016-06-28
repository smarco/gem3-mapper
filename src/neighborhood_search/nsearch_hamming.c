/*
 * PROJECT: GEMMapper
 * FILE: nsearch_hamming.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_hamming.h"
#include "neighborhood_search/nsearch_partition.h"
#include "neighborhood_search/nsearch_schedule.h"
#include "align/align.h"
#include "fm_index/fm_index_query.h"

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
  NSEARCH_PROF_ADD_NODE(nsearch_schedule);
#ifdef NSEARCH_ENUMERATE
  nsearch_schedule->nsearch_operation_aux->global_text[current_position] = char_enc;
  *lo_out = 0; *hi_out = 1;
#else
  fm_index_t* const fm_index = nsearch_schedule->search->archive->fm_index;
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
  NSEARCH_PROF_ADD_NODE(nsearch_schedule);
#ifdef NSEARCH_ENUMERATE
  nsearch_schedule->nsearch_operation_aux->global_text[current_position] = char_enc;
#else
  fm_index_t* const fm_index = nsearch_schedule->search->archive->fm_index;
  if (search_direction==direction_forward) {
    fm_index_2query_forward(fm_index,char_enc,fm_2interval_in,fm_2interval_out);
  } else {
    fm_index_2query_backward(fm_index,char_enc,fm_2interval_in,fm_2interval_out);
  }
#endif
}
void nsearch_hamming_terminate(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t lo,
    const uint64_t hi,
    const uint64_t align_distance) {
  NSEARCH_PROF_ADD_SOLUTION(nsearch_schedule);
#ifdef NSEARCH_ENUMERATE
  const uint8_t* const text = nsearch_schedule->nsearch_operation_aux->global_text;
  dna_buffer_print(stdout,text,nsearch_schedule->key_length,false);
  fprintf(stdout,"\n");
#else
  filtering_candidates_t* const filtering_candidates = nsearch_schedule->search->filtering_candidates;
  search_parameters_t* const search_parameters = nsearch_schedule->search->search_parameters;
  pattern_t* const pattern = &nsearch_schedule->search->pattern;
  filtering_candidates_add_region_interval(filtering_candidates,
      search_parameters,pattern,lo,hi,0,pattern->key_length,align_distance);
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
    const uint64_t next_error = current_error + (char_enc!=(nsearch_schedule->key[current_position]) ? 1 : 0);
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
}
void nsearch_hamming_brute_force(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Init
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_hamming,search,matches);
  // Search
  TIMER_START(&nsearch_schedule.profile.ns_timer);
#ifdef NSEARCH_ENUMERATE
  const uint64_t init_lo = 0;
  const uint64_t init_hi = 1;
#else
  const uint64_t init_lo = 0;
  const uint64_t init_hi = fm_index_get_length(nsearch_schedule.search->archive->fm_index);
#endif
  nsearch_hamming_brute_force_step(&nsearch_schedule,
      nsearch_schedule.key_length-1,0,init_lo,init_hi);
  TIMER_STOP(&nsearch_schedule.profile.ns_timer);
  // nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
/*
 * Scheduled search
 */
void nsearch_hamming_scheduled_search_operation(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,
    fm_2interval_t* const fm_2interval,
    const uint64_t local_error,
    const uint64_t global_error,
    const int64_t current_position) {
  // COmpute limits
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
    const uint64_t next_local_error = local_error + ((char_enc!=nsearch_schedule->key[current_position]) ? 1 : 0);
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
      nsearch_hamming_scheduled_search_operation(
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
        nsearch_hamming_scheduled_search_operation(
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
  fm_index_2query_init(nsearch_schedule->search->archive->fm_index,&fm_2interval);
#endif
  // Launch search
  // nsearch_schedule_print_pretty(stderr,nsearch_schedule); // DEBUG
  const int64_t begin_position =
      (nsearch_operation->search_direction==direction_forward) ?
          nsearch_operation->local_key_begin :
          nsearch_operation->local_key_end - 1;
  nsearch_hamming_scheduled_search_operation(nsearch_schedule,
      pending_searches,nsearch_operation,&fm_2interval,0,0,begin_position);
}
/*
 * Hamming Neighborhood Search
 */
void nsearch_hamming(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Search
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_hamming,search,matches);
  nsearch_schedule_search(&nsearch_schedule);
  // nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
/*
 * Neighborhood Search (Preconditioned by region profile)
 */
void nsearch_hamming_preconditioned(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Search
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_hamming,search,matches);
  nsearch_schedule_search_preconditioned(&nsearch_schedule);
  // nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
