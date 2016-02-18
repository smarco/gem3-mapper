/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_levenshtein.h"
#include "neighborhood_search/nsearch_levenshtein_state.h"
#include "neighborhood_search/nsearch_partition.h"
#include "data_structures/dna_text.h"

/*
 * Brute Force
 */
uint64_t nsearch_levenshtein_brute_print_solution(
    nsearch_schedule_t* const nsearch_schedule,const uint64_t text_position) {
  char* const search_string = nsearch_schedule->search_string; // Use first
  ++nsearch_schedule->ns_nodes_success; // PROFILE
  search_string[text_position+1] = '\0';
  fprintf(stdout,"%s\n",search_string);
  return 1;
}
uint64_t nsearch_levenshtein_brute_force_step_all(
    nsearch_schedule_t* const nsearch_schedule,const uint64_t text_position) {
  // Expand node for all characters
  nsearch_levenshtein_state_t* const nsearch_state = &nsearch_schedule->pending_searches->nsearch_state;
  char* const search_string = nsearch_schedule->search_string; // Use first
  const uint64_t max_error = nsearch_schedule->max_error;
  uint64_t total_matches_found = 0;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    uint64_t min_val, align_distance;
    nsearch_levenshtein_state_compute_chararacter(nsearch_state,true,nsearch_schedule->key,
        0,nsearch_schedule->key_length,text_position,char_enc,&min_val,&align_distance);
    if (min_val > max_error) continue;
    // Query character
    search_string[text_position] = dna_decode(char_enc);
    ++(nsearch_schedule->ns_nodes); // PROFILING
    // Keep searching
    uint64_t matches_found = 0;
    if (min_val <= align_distance && text_position < nsearch_schedule->max_text_length) {
      matches_found = nsearch_levenshtein_brute_force_step_all(nsearch_schedule,text_position+1);
    }
    if (align_distance <= max_error) {
      matches_found += nsearch_levenshtein_brute_print_solution(nsearch_schedule,text_position);
    }
    total_matches_found += matches_found;
  }
  if (total_matches_found==0) {
    ++nsearch_schedule->ns_nodes_closed;
    COUNTER_ADD(&nsearch_schedule->ns_nodes_closed_depth,text_position);
  }
  return total_matches_found;
}
uint64_t nsearch_levenshtein_brute_force_step_supercondensed(
    nsearch_schedule_t* const nsearch_schedule,const uint64_t text_position) {
  // Expand node for all characters
  nsearch_levenshtein_state_t* const nsearch_state = &nsearch_schedule->pending_searches->nsearch_state;
  char* const search_string = nsearch_schedule->search_string; // Use first
  const uint64_t max_error = nsearch_schedule->max_error;
  uint64_t total_matches_found = 0;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    uint64_t min_val, align_distance;
    nsearch_levenshtein_state_compute_chararacter(nsearch_state,true,nsearch_schedule->key,
        0,nsearch_schedule->key_length,text_position,char_enc,&min_val,&align_distance);
    if (min_val > max_error) continue;
    // Query character
    search_string[text_position] = dna_decode(char_enc);
    ++(nsearch_schedule->ns_nodes); // PROFILING
    // Keep searching
    if (align_distance <= max_error) {
      total_matches_found += nsearch_levenshtein_brute_print_solution(nsearch_schedule,text_position);
      // dp_matrix_print(stderr,&nsearch_state->dp_matrix,true,nsearch_schedule->key,0,2,nsearch_schedule->key,0,2);
    } else if (text_position < nsearch_schedule->max_text_length) {
      total_matches_found += nsearch_levenshtein_brute_force_step_supercondensed(nsearch_schedule,text_position+1);
    }
  }
  if (total_matches_found==0) {
    ++nsearch_schedule->ns_nodes_closed;
    COUNTER_ADD(&nsearch_schedule->ns_nodes_closed_depth,text_position);
  }
  return total_matches_found;
}
void nsearch_levenshtein_brute_force(
    fm_index_t* const fm_index,uint8_t* const key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  // Init
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_levenshtein,
      fm_index,NULL,key,key_length,max_error,intervals_result,mm_stack);
  // Search
  nsearch_levenshtein_state_t* const nsearch_state = &nsearch_schedule.pending_searches->nsearch_state;
  TIMER_START(&nsearch_schedule.ns_timer);
  // nsearch_levenshtein_state_prepare_full(nsearch_state,0,key_length,max_error);
  // nsearch_levenshtein_brute_force_step_full(&nsearch_schedule,0);
  nsearch_state->key_begin = 0;
  nsearch_state->key_end = key_length;
  nsearch_levenshtein_state_prepare_supercondensed(nsearch_state,max_error);
  nsearch_levenshtein_brute_force_step_supercondensed(&nsearch_schedule,0);
  TIMER_STOP(&nsearch_schedule.ns_timer);
  nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
/*
 * Perform Levenshtein Scheduled Search
 */
uint64_t nsearch_levenshtein_perform_scheduled_search_next_operation(
    nsearch_schedule_t* const nsearch_schedule,const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,const uint64_t global_error,
    const uint64_t align_distance) {
  // Parameters
  nsearch_levenshtein_state_t* const nsearch_state = &nsearch_operation->nsearch_state;
  // Check min-error constraints
  if (align_distance+global_error < nsearch_operation->min_local_error) return 0;
  if (align_distance+global_error < nsearch_operation->min_global_error) return 0;
  // Check pending searches
  if (pending_searches==0) {
    const bool forward_search = (nsearch_operation->search_direction==direction_forward);
    // nsearch_levenshtein_print_search_trace(stderr,nsearch_schedule,pending_searches);
    nsearch_levenshtein_state_print_search_text(stdout,nsearch_state,forward_search);
    ++nsearch_schedule->ns_nodes_success; // PROFILE
    return 1; // hi - lo;
  } else {
    // Next chunk
    const uint64_t next_pending_searches = pending_searches-1;
    nsearch_operation_t* const next_nsearch_operation = nsearch_schedule->pending_searches + next_pending_searches;
    nsearch_levenshtein_state_t* const next_nsearch_state = &next_nsearch_operation->nsearch_state;
    // Reconstruct DP
    const bool forward_search = (next_nsearch_operation->search_direction==direction_forward);
    next_nsearch_state->key_begin = next_nsearch_operation->begin;
    next_nsearch_state->key_end = next_nsearch_operation->end;
    nsearch_levenshtein_state_prepare_supercondensed(next_nsearch_state,next_nsearch_operation->max_local_error);
    nsearch_levenshtein_state_compute_sequence(nsearch_state,next_nsearch_state,forward_search);
    // Keep searching
    nsearch_levenshtein_perform_scheduled_search(
        nsearch_schedule,next_pending_searches,next_nsearch_operation,align_distance);
//    if (align_distance <= next_nsearch_operation->max_local_error) {
//      return nsearch_levenshtein_perform_scheduled_search_next_operation(
//          nsearch_schedule,next_pending_searches,next_nsearch_operation,align_distance);
//    } else if (text_length <= nsearch_schedule->max_text_length) {
//      return nsearch_levenshtein_perform_scheduled_search(
//          nsearch_schedule,next_pending_searches,next_nsearch_operation);
//    }
  }
  return 0;
}
uint64_t nsearch_levenshtein_perform_scheduled_search(
    nsearch_schedule_t* const nsearch_schedule,const uint64_t pending_searches,
    nsearch_operation_t* const nsearch_operation,const uint64_t global_error) {
  // Parameters
  nsearch_levenshtein_state_t* const nsearch_state = &nsearch_operation->nsearch_state;
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  const uint64_t text_position = nsearch_state->local_text_length;
  uint64_t total_matches_found = 0;
  // Search all characters (error >= key_length)
  const uint64_t key_length = nsearch_state->key_end-nsearch_state->key_begin;
  if (nsearch_operation->max_local_error >= key_length) {
    if (key_length > 1) {
      fprintf(stderr,"<<ERROR>>\n");
      GEM_NOT_IMPLEMENTED();
    }
    GEM_NOT_IMPLEMENTED();
//    uint8_t char_enc;
//    for (char_enc=0;char_enc<DNA_RANGE;++char_enc) { // TODO TODO TODO TODO TODO
//      const uint64_t align_distance ;
//      nsearch_state->local_text[text_position] = char_enc;
//      nsearch_state->local_text_length = text_position+1;
//      // Check max-error constraints
//      if (global_error+min_val > nsearch_operation->max_local_error) continue;
//      if (global_error+min_val > nsearch_operation->max_global_error) continue;
//    }
    return total_matches_found;
  }
  // Search all characters (expand node)
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    uint64_t min_val, align_distance;
    nsearch_levenshtein_state_compute_chararacter(
        nsearch_state,forward_search,nsearch_schedule->key,nsearch_state->key_begin,
        nsearch_state->key_end,text_position,char_enc,&min_val,&align_distance);
    nsearch_state->local_text[text_position] = char_enc;
    nsearch_state->local_text_length = text_position+1;
    // Check max-error constraints
    if (global_error+min_val > nsearch_operation->max_local_error) continue;
    if (global_error+min_val > nsearch_operation->max_global_error) continue;
    // Search character
    nsearch_state->local_text[text_position] = char_enc;
    nsearch_state->local_text_length = text_position+1;
    ++(nsearch_schedule->ns_nodes); // PROFILING
    // nsearch_levenshtein_state_print_text(stderr,nsearch_state,forward_search); // ACGGTTAC
    // nsearch_levenshtein_print_pair_key_text(stderr,nsearch_schedule,nsearch_operation);
    nsearch_levenshtein_print_search_trace(stderr,nsearch_schedule,pending_searches);
    // Keep searching
    uint64_t matches_found = 0;
    if (align_distance <= nsearch_operation->max_local_error &&
        global_error+align_distance <= nsearch_operation->max_global_error) {
      matches_found = nsearch_levenshtein_perform_scheduled_search_next_operation(
          nsearch_schedule,pending_searches,nsearch_operation,global_error,align_distance);
    }
    if ((matches_found==0 || pending_searches > 0) && text_position < nsearch_schedule->max_text_length) {
      matches_found += nsearch_levenshtein_perform_scheduled_search(
          nsearch_schedule,pending_searches,nsearch_operation,global_error);
    }
    total_matches_found += matches_found;
  }
  return total_matches_found;
}
/*
 * Hamming Neighborhood Search
 */
void nsearch_levenshtein(
    fm_index_t* const fm_index,uint8_t* const key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  // Push
  mm_stack_push_state(mm_stack);
  // Search
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_levenshtein,
      fm_index,NULL,key,key_length,max_error,intervals_result,mm_stack);
  nsearch_schedule_search(&nsearch_schedule);
  nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
  // Free
  mm_stack_pop_state(mm_stack);
}
/*
 * Neighborhood Search (Preconditioned by region profile)
 */
void nsearch_levenshtein_preconditioned(
    fm_index_t* const fm_index,region_profile_t* const region_profile,
    uint8_t* const key,const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  // Push
  mm_stack_push_state(mm_stack);
  // Search
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_levenshtein,
      fm_index,region_profile,key,key_length,max_error,intervals_result,mm_stack);
  nsearch_schedule_search_preconditioned(&nsearch_schedule);
  nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
  // Free
  mm_stack_pop_state(mm_stack);
}
/*
 * Display
 */
void nsearch_levenshtein_print_text(
    FILE* const stream,nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches) {
  nsearch_operation_t* nsearch_operation = NULL;
  nsearch_levenshtein_state_t* nsearch_state = NULL;
  const int64_t first_idx_op = pending_searches;
  const int64_t last_idx_op = nsearch_schedule->num_pending_searches-1;
  int64_t i;
  // Print searched string
  for (i=last_idx_op;i>=first_idx_op;--i) {
    nsearch_operation = nsearch_schedule->pending_searches + i;
    nsearch_state = &nsearch_operation->nsearch_state;
    fprintf(stream,"[%lu] key[%lu,%lu) ",last_idx_op-i,nsearch_state->key_begin,nsearch_state->key_end);
    nsearch_levenshtein_state_print_search_text(stream,nsearch_state,nsearch_operation->search_direction==direction_forward);
  }
  // Print current state
  nsearch_levenshtein_state_print(stream,nsearch_state,
      nsearch_operation->search_direction==direction_forward,nsearch_schedule->key);
}
void nsearch_levenshtein_print_search_trace(
    FILE* const stream,nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches) {
  nsearch_operation_t* nsearch_operation = NULL;
  nsearch_levenshtein_state_t* nsearch_state = NULL;
  const int64_t first_idx_op = pending_searches;
  const int64_t last_idx_op = nsearch_schedule->num_pending_searches-1;
  int64_t i;
  // Print searched string
  for (i=last_idx_op;i>=first_idx_op;--i) {
    nsearch_operation = nsearch_schedule->pending_searches + i;
    nsearch_state = &nsearch_operation->nsearch_state;
    fprintf(stream,"[%lu] key[%lu,%lu) Local-text=",last_idx_op-i,nsearch_state->key_begin,nsearch_state->key_end);
    nsearch_levenshtein_state_print_local_text(stream,
        nsearch_state,nsearch_operation->search_direction==direction_forward);
  }
  // Print current state
  nsearch_levenshtein_state_print(stream,nsearch_state,
      nsearch_operation->search_direction==direction_forward,nsearch_schedule->key);
}
//void nsearch_levenshtein_print_pair_key_text(
//    FILE* const stream,nsearch_schedule_t* const nsearch_schedule,
//    nsearch_operation_t* const nsearch_operation) {
//  nsearch_levenshtein_state_t* const nsearch_state = &nsearch_operation->nsearch_state;
//  const bool reverse_print = nsearch_operation->search_direction==direction_backward;
//  fprintf(stream,"(");
//  dna_buffer_print(stream,nsearch_state->text+nsearch_state->key_begin,
//      nsearch_state->key_end-nsearch_state->key_begin,reverse_print);
//  fprintf(stream,",");
//  dna_buffer_print(stream,nsearch_schedule->key,nsearch_schedule->key_length,reverse_print);
//  fprintf(stream,")\n");
//}
