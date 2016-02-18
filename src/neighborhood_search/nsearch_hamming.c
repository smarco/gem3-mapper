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
#include "data_structures/dna_text.h"

/*
 * Brute Force
 */
void nsearch_hamming_brute_force_step(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t current_position,const uint64_t current_error) {
  char* const search_string = nsearch_schedule->search_string; // Use first
  uint8_t enc;
  for (enc=0;enc<DNA_RANGE;++enc) {
    const uint64_t next_error = current_error + (enc!=(nsearch_schedule->key[current_position]) ? 1 : 0);
    if (next_error <= nsearch_schedule->max_error) {
      search_string[current_position] = dna_decode(enc);
      search_string[current_position+1] = '\0';
      ++(nsearch_schedule->ns_nodes); // PROFILING
      const uint64_t next_position = current_position + 1;
      if (next_position < nsearch_schedule->key_length) {
        nsearch_hamming_brute_force_step(nsearch_schedule,next_position,next_error);
      } else {
        ++(nsearch_schedule->ns_nodes_success);
        fprintf(stdout,"%s\n",search_string);
      }
    }
  }
}
void nsearch_hamming_brute_force(
    fm_index_t* const fm_index,uint8_t* const key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  // Init
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_hamming,fm_index,
      NULL,key,key_length,max_error,intervals_result,mm_stack);
  // Search
  TIMER_START(&nsearch_schedule.ns_timer);
  nsearch_hamming_brute_force_step(&nsearch_schedule,0,0);
  TIMER_STOP(&nsearch_schedule.ns_timer);
  nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
/*
 * Perform scheduled search
 */
void nsearch_hamming_perform_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule,const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,const uint64_t position,
    const uint64_t local_error,const uint64_t global_error) {
  const uint64_t next_position = position + 1;
  uint8_t enc;
  for (enc=0;enc<DNA_RANGE;++enc) {
    // Compute distance function & constraints
    const uint64_t next_local_error = local_error + ((enc!=nsearch_schedule->key[position]) ? 1 : 0);
    if (next_local_error > nsearch_operation->max_local_error) continue;
    if (next_local_error+global_error > nsearch_operation->max_global_error) continue;
    if (next_position==nsearch_operation->end) {
      if (next_local_error < nsearch_operation->min_local_error) continue;
      if (next_local_error+global_error < nsearch_operation->min_global_error) continue;
    }
    // Search character
    const uint64_t search_string_pos = position-nsearch_operation->begin;
    nsearch_schedule->search_string[nsearch_operation->begin+search_string_pos] = dna_decode(enc);
    nsearch_schedule->search_string[nsearch_operation->begin+search_string_pos+1] = '\0';
    ++(nsearch_schedule->ns_nodes); // PROFILING
    // Keep searching
    if (next_position < nsearch_operation->end) {
      nsearch_hamming_perform_scheduled_search(nsearch_schedule,pending_searches,
          nsearch_operation,next_position,next_local_error,global_error);
    } else { // (next_position == limit) End of the chunk-search
      // Check pending searches
      if (pending_searches==0) {
        // End of the search (Print search-string)
        ++(nsearch_schedule->ns_nodes_success);
        nsearch_schedule_print_search_string(stdout,nsearch_schedule);
      } else {
        // Search next chunk
        const uint64_t next_pending_searches = pending_searches-1;
        nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
        nsearch_hamming_perform_scheduled_search(nsearch_schedule,next_pending_searches,
            next_nsearch_operation,next_nsearch_operation->begin,0,next_local_error+global_error);
      }
    }
  }
}
/*
 * Hamming Neighborhood Search
 */
void nsearch_hamming(
    fm_index_t* const fm_index,uint8_t* const key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  // Push
  mm_stack_push_state(mm_stack);
  // Search
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_hamming,fm_index,
      NULL,key,key_length,max_error,intervals_result,mm_stack);
  nsearch_schedule_search(&nsearch_schedule);
  nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
  // Free
  mm_stack_pop_state(mm_stack);
}
/*
 * Neighborhood Search (Preconditioned by region profile)
 */
void nsearch_hamming_preconditioned(
    fm_index_t* const fm_index,region_profile_t* const region_profile,
    uint8_t* const key,const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  // Push
  mm_stack_push_state(mm_stack);
  // Search
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_hamming,fm_index,
      region_profile,key,key_length,max_error,intervals_result,mm_stack);
  nsearch_schedule_search_preconditioned(&nsearch_schedule);
  nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
  // Free
  mm_stack_pop_state(mm_stack);
}
