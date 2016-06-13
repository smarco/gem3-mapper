/*
 * PROJECT: GEMMapper
 * FILE: nsearch_operation.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_operation.h"
#include "data_structures/dna_text.h"

/*
 * Setup
 */
void nsearch_operation_init(
    nsearch_operation_t* const nsearch_operation,
    const uint64_t max_key_length,
    const uint64_t max_text_length,
    mm_stack_t* const mm_stack) {
  // Compute dimensions
  const uint64_t num_rows = max_key_length + 1;     // (+1) Initial Conditions Row
  const uint64_t num_columns = max_text_length + 2; // (+1) Initial Conditions Column
                                                    // (+1) Auxiliary Column (check there is no valid match)
  // Allocate text
  nsearch_operation->text = mm_stack_calloc(mm_stack,num_columns,uint8_t,false);
  nsearch_operation->global_text = mm_stack_calloc(mm_stack,num_columns,uint8_t,false);
  // Init levenshtein state
  nsearch_levenshtein_state_init(&nsearch_operation->nsearch_state,num_rows,num_columns,mm_stack);
}
/*
 * Prepare Sequence
 */
void nsearch_operation_chained_prepare_sequence_forward(
    nsearch_operation_t* const nsearch_operation,
    uint8_t* const text,
    uint64_t* const text_length) {
  // Copy global-text
  const uint64_t global_text_length = nsearch_operation->global_text_length;
  uint8_t* const global_text = nsearch_operation->global_text;
  uint64_t i;
  for (i=0;i<global_text_length;++i) text[i] = global_text[i];
  // Copy local-text
  const uint64_t local_text_length = nsearch_operation->text_position;
  uint8_t* const local_text = nsearch_operation->text;
  for (i=0;i<local_text_length;++i) text[global_text_length+i] = local_text[i];
  // Set global-text length
  *text_length = global_text_length + local_text_length;
}
void nsearch_operation_chained_prepare_sequence_reverse(
    nsearch_operation_t* const nsearch_operation,
    uint8_t* const text,
    uint64_t* const text_length) {
  // Copy local-text
  const uint64_t local_text_length = nsearch_operation->text_position;
  const uint64_t local_text_last_idx = local_text_length-1;
  uint8_t* const local_text = nsearch_operation->text;
  uint64_t i;
  for (i=0;i<local_text_length;++i) text[i] = local_text[local_text_last_idx-i];
  // Copy global-text
  const uint64_t global_text_length = nsearch_operation->global_text_length;
  const uint64_t global_text_last_idx = global_text_length-1;
  uint8_t* const global_text = nsearch_operation->global_text;
  for (i=0;i<global_text_length;++i) text[local_text_length+i] = global_text[global_text_last_idx-i];
  // Set global-text length
  *text_length = global_text_length + local_text_length;
}
/*
 * Prepare Operation
 */
void nsearch_operation_chained_prepare_forward(
    nsearch_operation_t* const current_nsearch_operation,
    nsearch_operation_t* const next_nsearch_operation,
    uint8_t* const key,
    const uint64_t key_length) {
  // Parameters
  nsearch_levenshtein_state_t* const current_nsearch_state = &current_nsearch_operation->nsearch_state;
  nsearch_levenshtein_state_t* const next_nsearch_state = &next_nsearch_operation->nsearch_state;
  // Prepare context (Isolate active columns)
  // fprintf(stderr,"[GEM]> Operation.Chained (FORWARD) >>>>>>>>>>>>>>>>>>> \n");
  // nsearch_operation_state_print(stderr,current_nsearch_operation,key); // DEBUG
  nsearch_operation_chained_prepare_sequence_forward(current_nsearch_operation,
      next_nsearch_operation->global_text,&next_nsearch_operation->global_text_length);
  next_nsearch_operation->text_position = 0;
  // Chain search (Same direction)
  const uint64_t current_key_length =
      current_nsearch_operation->global_key_end -
      current_nsearch_operation->global_key_begin;
  const uint64_t current_text_length =
      current_nsearch_operation->text_position;
  const bool next_operation_towards_end =
      (next_nsearch_operation->search_direction==direction_forward) ?
          (next_nsearch_operation->global_key_begin == 0) :
          (next_nsearch_operation->global_key_end == key_length);
  nsearch_levenshtein_state_prepare_chained(
      current_nsearch_state,next_nsearch_state,
      current_key_length,current_text_length,next_operation_towards_end);
  // nsearch_operation_state_print(stderr,next_nsearch_operation,key); // DEBUG
  // fprintf(stderr,"[GEM]> Operation.Chained (FORWARD) <<<<<<<<<<<<<<<<<<< \n");
}
void nsearch_operation_compute_reverse(
    nsearch_operation_t* const nsearch_operation,
    nsearch_operation_t* const nsearch_operation_rev,
    uint8_t* const key,
    const uint64_t key_length) {
  // Parameters
  const uint8_t* const key_chunk = key + nsearch_operation->global_key_begin;
  const uint64_t key_chunk_length = nsearch_operation->global_key_end - nsearch_operation->global_key_begin;
  const bool forward_search = (nsearch_operation->search_direction != direction_forward);
  const uint64_t current_max_error = nsearch_operation->max_global_error;
  // Reverse sequence
//  fprintf(stderr,"[GEM]> Operation.Chained (REVERSE) >>>>>>>>>>>>>>>>>>> \n");
//  nsearch_operation_state_print(stderr,nsearch_operation,key); // DEBUG
  nsearch_operation_chained_prepare_sequence_reverse(nsearch_operation,
      nsearch_operation_rev->text,&nsearch_operation_rev->text_position);
  // Reverse DP (recompute DP & isolate active columns)
  const bool next_operation_towards_end =
      (nsearch_operation->search_direction==direction_backward) ?
          (nsearch_operation->global_key_begin == 0) :
          (nsearch_operation->global_key_end == key_length);
  if (next_operation_towards_end) {
    nsearch_levenshtein_state_prepare_supercondensed_neighbourhood(&nsearch_operation_rev->nsearch_state);
  } else {
    nsearch_levenshtein_state_prepare_full_neighbourhood(&nsearch_operation_rev->nsearch_state);
  }
  const uint8_t* const text = nsearch_operation_rev->text;
  const uint64_t text_length = nsearch_operation_rev->text_position;
  nsearch_levenshtein_state_compute_text(
      &nsearch_operation_rev->nsearch_state,forward_search,
      key_chunk,key_chunk_length,text,text_length,current_max_error);
  // Copy necessary fields
  nsearch_operation_rev->search_direction = (forward_search) ? direction_forward : direction_backward;
  nsearch_operation_rev->min_global_error = nsearch_operation->min_global_error;
  nsearch_operation_rev->global_key_begin = nsearch_operation->global_key_begin;
  nsearch_operation_rev->global_key_end = nsearch_operation->global_key_end;
//  // DEBUG
//  dp_matrix_print(
//      stderr,&nsearch_operation_rev->nsearch_state.dp_matrix,
//      forward_search,key,nsearch_operation->global_key_begin,
//      nsearch_operation->global_key_end,text,0,text_length);
//  fprintf(stderr,"[GEM]> Operation.Chained (REVERSE) <<<<<<<<<<<<<<<<<<< \n");
}
/*
 * Utils
 */
int nsearch_operation_state_global_text_cmp(
    nsearch_operation_t* const nsearch_operation,
    char* const global_text,
    mm_stack_t* const mm_stack) {
  mm_stack_push_state(mm_stack);
  const uint64_t global_text_length = strlen(global_text);
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  // Copy global text
  uint8_t* const operation_global_text = mm_stack_calloc(mm_stack,global_text_length,uint8_t,true);
  uint64_t i, j;
  for (i=0;i<nsearch_operation->global_text_length;++i) {
    operation_global_text[i] = nsearch_operation->global_text[i];
  }
  const uint8_t* const local_text = nsearch_operation->text;
  const uint64_t local_text_length = nsearch_operation->text_position;
  for (j=0;j<local_text_length;++j,++i) {
    operation_global_text[i] = local_text[j];
  }
  // Reverse if needed
  const uint64_t operation_global_text_length = i;
  if (!forward_search) {
    for (j=0;j<operation_global_text_length/2;++j) {
      SWAP(operation_global_text[j],operation_global_text[operation_global_text_length-j-1]);
    }
  }
  // Return comparison
  int eq = global_text_length == operation_global_text_length;
  for (j=0;eq && j<operation_global_text_length;++j) {
    eq = operation_global_text[j] == dna_encode(global_text[j]);
  }
  mm_stack_pop_state(mm_stack);
  return !eq;
}
/*
 * Display
 */
void nsearch_operation_print(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation) {
  fprintf(stream,">>LKey[%lu,%lu)/GKey[%lu,%lu) LError{%lu,%lu}/GError(%lu,%lu)\n",
      nsearch_operation->local_key_begin,nsearch_operation->local_key_end,
      nsearch_operation->global_key_begin,nsearch_operation->global_key_end,
      nsearch_operation->min_local_error,nsearch_operation->max_local_error,
      nsearch_operation->min_global_error,nsearch_operation->max_global_error);
}
void nsearch_operation_state_print(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation,
    const uint8_t* const key) {
  const bool forward_search = (nsearch_operation->search_direction == direction_forward);
  fprintf(stream,"[GEM]> Levenshtein.State\n");
  fprintf(stream,"  => Search %s\n",forward_search ? "forward" : "reverse");
  fprintf(stream,"  => Local-Key [%lu,%lu)\n",nsearch_operation->local_key_begin,nsearch_operation->local_key_end);
  fprintf(stream,"  => Global-Key [%lu,%lu)\n",nsearch_operation->global_key_begin,nsearch_operation->global_key_end);
  fprintf(stream,"  => Global-Text ");
  nsearch_operation_state_print_global_text(stream,nsearch_operation);
  fprintf(stream,"  => Local-Text  ");
  nsearch_operation_state_print_local_text(stream,nsearch_operation);
  fprintf(stream,"  => DP-Matrix\n");
  dp_matrix_print(
      stream,&nsearch_operation->nsearch_state.dp_matrix,forward_search,
      key,nsearch_operation->global_key_begin,nsearch_operation->global_key_end,
      nsearch_operation->text,0,nsearch_operation->text_position);
}
void nsearch_operation_state_print_global_text(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation) {
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  if (forward_search) {
    dna_buffer_print(stream,nsearch_operation->global_text,nsearch_operation->global_text_length,false);
    dna_buffer_print(stream,nsearch_operation->text,nsearch_operation->text_position,false);
  } else {
    dna_buffer_print(stream,nsearch_operation->text,nsearch_operation->text_position,true);
    dna_buffer_print(stream,nsearch_operation->global_text,nsearch_operation->global_text_length,true);
  }
  fprintf(stream,"\n");
}
void nsearch_operation_state_print_local_text(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation) {
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  const uint8_t* const local_text = nsearch_operation->text;
  const uint64_t local_text_length = nsearch_operation->text_position;
  dna_buffer_print(stream,local_text,local_text_length,!forward_search);
  fprintf(stream,"\n");
}
