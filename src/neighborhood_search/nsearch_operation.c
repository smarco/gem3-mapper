/*
 * PROJECT: GEMMapper
 * FILE: nsearch_operation.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "text/dna_text.h"
#include "neighborhood_search/nsearch_operation.h"

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
  // Init levenshtein state
  nsearch_levenshtein_state_init(&nsearch_operation->nsearch_state,num_rows,num_columns,mm_stack);
}
/*
 * Prepare Operation Chained
 */
void nsearch_operation_chained_prepare_sequence_forward(
    nsearch_operation_t* const nsearch_operation,
    uint8_t* const text,
    uint64_t* const text_length) {
  // Copy text
  const uint64_t operation_text_length = nsearch_operation->text_position;
  uint8_t* const operation_text = nsearch_operation->text;
  memcpy(text,operation_text,operation_text_length);
  *text_length = operation_text_length; // Set length
}
void nsearch_operation_chained_prepare_sequence_reverse(
    nsearch_operation_t* const nsearch_operation,
    uint8_t* const text,
    uint64_t* const text_length) {
  // Copy reversed text
  const uint64_t operation_text_length = nsearch_operation->text_position;
  const uint64_t operation_text_last_idx = operation_text_length-1;
  uint8_t* const operation_text = nsearch_operation->text;
  uint64_t i;
  for (i=0;i<operation_text_length;++i) text[i] = operation_text[operation_text_last_idx-i];
  *text_length = operation_text_length; // Set global-text length
}
bool nsearch_operation_chained_prepare(
    nsearch_operation_t* const current_nsearch_operation,
    nsearch_operation_t* const next_nsearch_operation,
    uint8_t* const key,
    const uint64_t key_length,
    const bool reverse_sequence,
    mm_stack_t* const mm_stack) {
  // Prepare Sequence
  if (reverse_sequence) {
    nsearch_operation_chained_prepare_sequence_reverse(current_nsearch_operation,
        next_nsearch_operation->text,&next_nsearch_operation->text_position);
  } else {
    nsearch_operation_chained_prepare_sequence_forward(current_nsearch_operation,
        next_nsearch_operation->text,&next_nsearch_operation->text_position);
  }
  // Recompute DP
  nsearch_levenshtein_state_t* const next_nsearch_state = &next_nsearch_operation->nsearch_state;
  const uint8_t* const key_chunk = key + next_nsearch_operation->global_key_begin;
  const uint64_t key_chunk_length = next_nsearch_operation->global_key_end - next_nsearch_operation->global_key_begin;
  const bool forward_search = (next_nsearch_operation->search_direction == direction_forward);
  const uint8_t* const text = next_nsearch_operation->text;
  const uint64_t text_length = next_nsearch_operation->text_position;
  const uint64_t next_max_error = next_nsearch_operation->max_global_error;
  uint64_t min_align_distance, min_align_distance_column;
  nsearch_levenshtein_state_compute_text_banded(
      next_nsearch_state,forward_search,key_chunk,
      key_chunk_length,text,text_length,next_max_error,
      &min_align_distance,&min_align_distance_column,mm_stack);
  // Check supercondensed operation-chain
  return !(key_chunk_length==key_length &&
           min_align_distance<=next_max_error &&
           min_align_distance_column<key_length);
}
/*
 * Utils
 */
bool nsearch_operation_state_text_eq(
    nsearch_operation_t* const nsearch_operation,
    char* const text,
    mm_stack_t* const mm_stack) {
  mm_stack_push_state(mm_stack);
  const uint64_t text_length = strlen(text);
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  // Allocate text
  uint8_t* const enc_text = mm_stack_calloc(mm_stack,text_length,uint8_t,true);
  // Encode
  uint64_t i;
  if (forward_search) {
    for (i=0;i<text_length;++i) {
      enc_text[i] = dna_encode(text[i]);
    }
  } else {
    for (i=0;i<text_length;++i) {
      enc_text[text_length-i-1] = dna_encode(text[i]);
    }
  }
  // Return comparison
  bool eq = (text_length == nsearch_operation->text_position);
  for (i=0;eq && i<text_length;++i) {
    eq = (nsearch_operation->text[i] == enc_text[i]);
  }
  mm_stack_pop_state(mm_stack);
  return eq;
}
/*
 * Display
 */
void nsearch_operation_print(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation) {
  fprintf(stream,">>LKey[%lu,%lu)/GKey[%lu,%lu) Error{%lu/%lu,%lu}\n",
      nsearch_operation->local_key_begin,nsearch_operation->local_key_end,
      nsearch_operation->global_key_begin,nsearch_operation->global_key_end,
      nsearch_operation->min_local_error,nsearch_operation->min_global_error,
      nsearch_operation->max_global_error);
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
  dna_buffer_print(stream,nsearch_operation->text,nsearch_operation->text_position,!forward_search);
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
