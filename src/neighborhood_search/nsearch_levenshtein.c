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

/*
 * Levenshtein Brute Force
 */
uint64_t nsearch_levenshtein_brute_print_solution(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t text_position) {
  char* const search_string = nsearch_schedule->search_string; // Use first
  ++nsearch_schedule->ns_nodes_success; // PROFILE
  search_string[text_position+1] = '\0';
  fprintf(stdout,"%s\n",search_string);
  return 1;
}
uint64_t nsearch_levenshtein_brute_force_step_full_neighbourhood(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation) {
  // Parameters
  char* const search_string = nsearch_schedule->search_string; // Use first
  const uint64_t max_error = nsearch_schedule->max_error;
  const uint8_t* const key = nsearch_schedule->key;
  const uint64_t key_length = nsearch_schedule->key_length;
  const uint64_t text_position = nsearch_operation->text_position;
  // Expand node for all characters
  uint64_t total_matches_found = 0;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    uint64_t min_val, align_distance;
    nsearch_levenshtein_state_compute_chararacter(&nsearch_operation->nsearch_state,true,
        key,key_length,text_position,char_enc,max_error,&min_val,&align_distance);
    if (min_val > max_error) continue;
    // Query character
    search_string[text_position] = dna_decode(char_enc);
    nsearch_operation->text_position = text_position + 1;
    ++(nsearch_schedule->ns_nodes); // PROFILING
    // Keep searching
    uint64_t matches_found = 0;
    if (min_val <= align_distance && text_position < nsearch_schedule->max_text_length) {
      matches_found = nsearch_levenshtein_brute_force_step_full_neighbourhood(nsearch_schedule,nsearch_operation);
    }
    if (align_distance <= max_error) {
//      dp_matrix_print(
//          stderr,&nsearch_operation->nsearch_state.dp_matrix,true,
//          key,0,key_length,search_string,0,text_position + 1);
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
uint64_t nsearch_levenshtein_brute_force_step_supercondensed_neighbourhood(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation) {
  // Parameters
  char* const search_string = nsearch_schedule->search_string; // Use first
  const uint64_t max_error = nsearch_schedule->max_error;
  const uint8_t* const key = nsearch_schedule->key;
  const uint64_t key_length = nsearch_schedule->key_length;
  const uint64_t text_position = nsearch_operation->text_position;
  // Expand node for all characters
  uint64_t total_matches_found = 0;
  uint8_t char_enc;
  for (char_enc=0;char_enc<DNA_RANGE;++char_enc) {
    // Compute DP-next
    uint64_t min_val, align_distance;
    nsearch_levenshtein_state_compute_chararacter(&nsearch_operation->nsearch_state,true,
        key,key_length,text_position,char_enc,max_error,&min_val,&align_distance);
    if (min_val > max_error) continue;
    // Query character
    search_string[text_position] = dna_decode(char_enc);
    nsearch_operation->text_position = text_position + 1;
    ++(nsearch_schedule->ns_nodes); // PROFILING
    // Keep searching
    if (align_distance <= max_error) {
      total_matches_found += nsearch_levenshtein_brute_print_solution(nsearch_schedule,text_position);
      // dp_matrix_print(stderr,&nsearch_state->dp_matrix,true,nsearch_schedule->key,0,2,nsearch_schedule->key,0,2);
    } else if (text_position < nsearch_schedule->max_text_length) {
      total_matches_found += nsearch_levenshtein_brute_force_step_supercondensed_neighbourhood(
          nsearch_schedule,nsearch_operation);
    }
  }
  if (total_matches_found==0) {
    ++nsearch_schedule->ns_nodes_closed;
    COUNTER_ADD(&nsearch_schedule->ns_nodes_closed_depth,text_position);
  }
  return total_matches_found;
}
void nsearch_levenshtein_brute_force_full(
    fm_index_t* const fm_index,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t max_error,
    interval_set_t* const intervals_result,
    mm_stack_t* const mm_stack) {
  // Init
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_levenshtein,
      fm_index,NULL,key,key_length,max_error,intervals_result,mm_stack);
  // Search
  TIMER_START(&nsearch_schedule.ns_timer);
  nsearch_operation_t* const nsearch_operation = nsearch_schedule.pending_searches;
  nsearch_levenshtein_state_prepare_full_neighbourhood(&nsearch_operation->nsearch_state);
  nsearch_levenshtein_brute_force_step_full_neighbourhood(&nsearch_schedule,nsearch_operation);
  TIMER_STOP(&nsearch_schedule.ns_timer);
  nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
void nsearch_levenshtein_brute_force_supercondensed(
    fm_index_t* const fm_index,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t max_error,
    interval_set_t* const intervals_result,
    mm_stack_t* const mm_stack) {
  // Init
  nsearch_schedule_t nsearch_schedule;
  nsearch_schedule_init(&nsearch_schedule,nsearch_model_levenshtein,
      fm_index,NULL,key,key_length,max_error,intervals_result,mm_stack);
  // Search
  TIMER_START(&nsearch_schedule.ns_timer);
  nsearch_operation_t* const nsearch_operation = nsearch_schedule.pending_searches;
  nsearch_levenshtein_state_prepare_supercondensed_neighbourhood(&nsearch_operation->nsearch_state);
  nsearch_levenshtein_brute_force_step_supercondensed_neighbourhood(&nsearch_schedule,nsearch_operation);
  TIMER_STOP(&nsearch_schedule.ns_timer);
  nsearch_schedule_print_profile(stderr,&nsearch_schedule); // PROFILE
}
/*
 * Levenshtein Neighborhood Search
 */
void nsearch_levenshtein(
    fm_index_t* const fm_index,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t max_error,
    interval_set_t* const intervals_result,
    mm_stack_t* const mm_stack) {
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
    fm_index_t* const fm_index,
    region_profile_t* const region_profile,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t max_error,
    interval_set_t* const intervals_result,
    mm_stack_t* const mm_stack) {
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
void nsearch_levenshtein_print_status(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t pending_searches) {
  nsearch_operation_t* nsearch_operation = NULL;
  const int64_t first_idx_op = pending_searches;
  const int64_t last_idx_op = nsearch_schedule->num_pending_searches-1;
  int64_t i;
  // Print searched string
  for (i=last_idx_op;i>=first_idx_op;--i) {
    nsearch_operation = nsearch_schedule->pending_searches + i;
    fprintf(stream,"[%lu] key[%lu,%lu) Local-text=",last_idx_op-i,
        nsearch_operation->global_key_begin,nsearch_operation->global_key_end);
    nsearch_operation_state_print_local_text(stream,nsearch_operation);
    nsearch_operation_state_print(stream,nsearch_operation,nsearch_schedule->key);
  }
//  // Print current state
//  nsearch_operation_state_print(stream,nsearch_operation,nsearch_schedule->key);
}
typedef struct {
  nsearch_operation_t* nsearch_operation;
} nsearch_levenshtein_trace_t;
int nsearch_levenshtein_trace_sort(
    const nsearch_levenshtein_trace_t* const a,
    const nsearch_levenshtein_trace_t* const b) {
  return a->nsearch_operation->local_key_begin - b->nsearch_operation->local_key_begin;
}
#define VECTOR_SORT_NAME                 nsearch_levenshtein_trace
#define VECTOR_SORT_TYPE                 nsearch_levenshtein_trace_t
#define VECTOR_SORT_CMP(a,b)             nsearch_levenshtein_trace_sort(a,b)
#include "utils/vector_sort.h"
void nsearch_levenshtein_print_trace(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule) {
  // Save stack state & allocate mem
  mm_stack_t* const mm_stack = nsearch_schedule->mm_stack;
  mm_stack_push_state(mm_stack);
  const uint64_t num_searches = nsearch_schedule->num_pending_searches;
  nsearch_levenshtein_trace_t* const nsearch_levenshtein_trace =
      mm_stack_calloc(mm_stack,num_searches,nsearch_levenshtein_trace_t,true);
  // Copy all operations
  uint64_t i;
  for (i=0;i<num_searches;++i) {
    nsearch_levenshtein_trace[i].nsearch_operation = nsearch_schedule->pending_searches + i;
  }
  // Sort by local-key position
  buffer_sort_nsearch_levenshtein_trace(nsearch_levenshtein_trace,num_searches);
  // Print global-text
  nsearch_operation_t* const last_nsearch_operation = nsearch_schedule->pending_searches;
  if (last_nsearch_operation->search_direction==direction_forward) {
    dna_buffer_print(stream,last_nsearch_operation->global_text,last_nsearch_operation->global_text_length,false);
    dna_buffer_print(stream,last_nsearch_operation->text,last_nsearch_operation->text_position,false);
  } else {
    dna_buffer_print(stream,last_nsearch_operation->text,last_nsearch_operation->text_position,true);
    dna_buffer_print(stream,last_nsearch_operation->global_text,last_nsearch_operation->global_text_length,true);
  }
  fprintf(stream,"\t");
  // Print global-text with region separators
  for (i=0;i<num_searches;++i) {
    nsearch_operation_t* const nsearch_operation = nsearch_levenshtein_trace[i].nsearch_operation;
    if (i>0) fprintf(stream,"|");
    dna_buffer_print(stream,
        nsearch_operation->text,nsearch_operation->text_position,
        nsearch_operation->search_direction!=direction_forward);
  }
  fprintf(stream,"\t");
  // Print search-ranges min/max-error
  for (i=0;i<num_searches;++i) {
    nsearch_operation_t* const nsearch_operation = nsearch_levenshtein_trace[i].nsearch_operation;
    fprintf(stream,"{%lu,%lu}",nsearch_operation->min_local_error,nsearch_operation->max_local_error);
  }
  fprintf(stream,"\t");
  // Print search-ranges trace
  for (i=0;i<num_searches;++i) {
    nsearch_operation_t* const nsearch_operation = nsearch_schedule->pending_searches + (num_searches-i-1);
    fprintf(stream,"[%lu,%lu)",nsearch_operation->local_key_begin,nsearch_operation->local_key_end);
  }
  fprintf(stream,"\n");
  // Free
  mm_stack_pop_state(mm_stack);
}
