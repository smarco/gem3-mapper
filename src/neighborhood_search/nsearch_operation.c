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
void nsearch_operation_chained_prepare_sequence(
    nsearch_operation_t* const current_operation,
    nsearch_operation_t* const next_operation) {
  // Copy global-text
  uint8_t* const next_global_text = next_operation->global_text;
  const uint64_t global_text_length = current_operation->global_text_length;
  uint8_t* const global_text = current_operation->global_text;
  uint64_t i;
  for (i=0;i<global_text_length;++i) next_global_text[i] = global_text[i];
  // Copy local-text
  const uint64_t local_text_length = current_operation->text_position - current_operation->context_length;
  uint8_t* const local_text = current_operation->text + current_operation->context_length;
  for (i=0;i<local_text_length;++i) next_global_text[global_text_length+i] = local_text[i];
  // Set global-text length
  next_operation->global_text_length = global_text_length + local_text_length;
}
void nsearch_operation_chained_prepare_context(
    nsearch_operation_t* const current_operation,
    nsearch_operation_t* const next_operation,
    const uint64_t first_active_column) {
  // Copy context
  const uint64_t last_active_column = current_operation->text_position;
  const uint64_t context_length = last_active_column - first_active_column + 1;
  uint8_t* const local_text = current_operation->text + current_operation->context_length;
  uint8_t* const next_local_text = next_operation->text;
  uint64_t i;
  for (i=0;i<context_length;++i) next_local_text[i] = local_text[first_active_column+i-1];
  next_operation->context_length = context_length;
  next_operation->text_position = context_length;
}
void nsearch_operation_chained_prepare_reverse_sequence(
    nsearch_operation_t* const current_operation,
    nsearch_operation_t* const next_operation) {
  // Parameters
  uint8_t* const next_global_text = next_operation->global_text;
  // Copy local-text
  const uint64_t local_text_length = current_operation->text_position - current_operation->context_length;
  const uint64_t local_text_last_idx = local_text_length-1;
  uint8_t* const local_text = current_operation->text + current_operation->context_length;
  uint64_t i;
  for (i=0;i<local_text_length;++i) next_global_text[i] = local_text[local_text_last_idx-i];
  // Copy global-text
  const uint64_t global_text_length = current_operation->global_text_length;
  const uint64_t global_text_last_idx = global_text_length-1;
  uint8_t* const global_text = current_operation->global_text;
  for (i=0;i<global_text_length;++i) next_global_text[local_text_length+i] = global_text[global_text_last_idx-i];
  // Set global-text length
  next_operation->global_text_length = global_text_length + local_text_length;
}
void nsearch_operation_chained_prepare_reverse_context(
    nsearch_operation_t* const next_operation,
    const uint64_t first_active_column,
    const uint64_t last_active_column) {
  // Copy context
  const uint64_t context_length = last_active_column - first_active_column + 1;
  uint8_t* const next_local_text = next_operation->text;
  uint64_t i;
  for (i=0;i<context_length;++i) next_local_text[i] = next_local_text[first_active_column+i];
  next_operation->context_length = context_length;
  next_operation->text_position = context_length;
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
    const bool forward_search,
    const uint8_t* const key) {
  fprintf(stream,"[GEM]> Levenshtein.State\n");
  fprintf(stream,"  => Search %s\n",forward_search ? "forward" : "reverse");
  fprintf(stream,"  => Local-Key [%lu,%lu)\n",nsearch_operation->local_key_begin,nsearch_operation->local_key_end);
  fprintf(stream,"  => Global-Key [%lu,%lu)\n",nsearch_operation->global_key_begin,nsearch_operation->global_key_end);
  fprintf(stream,"  => Global-Text ");
  nsearch_operation_state_print_global_text(stream,nsearch_operation,forward_search);
  fprintf(stream,"  => Local-Text  ");
  nsearch_operation_state_print_local_text(stream,nsearch_operation,forward_search);
  fprintf(stream,"  => DP-Matrix\n");
  dp_matrix_print(
      stream,&nsearch_operation->nsearch_state.dp_matrix,forward_search,
      key,nsearch_operation->global_key_begin,nsearch_operation->global_key_end,
      nsearch_operation->text,0,nsearch_operation->text_position);
}
void nsearch_operation_state_print_global_text(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search) {
  if (forward_search) {
    dna_buffer_print(stream,nsearch_operation->global_text,nsearch_operation->global_text_length,false);
    dna_buffer_print(stream,nsearch_operation->text+nsearch_operation->context_length,
        nsearch_operation->text_position-nsearch_operation->context_length,false);
  } else {
    dna_buffer_print(stream,nsearch_operation->text+nsearch_operation->context_length,
        nsearch_operation->text_position-nsearch_operation->context_length,true);
    dna_buffer_print(stream,nsearch_operation->global_text,nsearch_operation->global_text_length,true);
  }
  fprintf(stream,"\n");
}
void nsearch_operation_state_print_local_text(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search) {
  const uint8_t* const local_text = nsearch_operation->text+nsearch_operation->context_length;
  const uint64_t local_text_length = nsearch_operation->text_position-nsearch_operation->context_length;
  dna_buffer_print(stream,local_text,local_text_length,!forward_search);
  fprintf(stream,"\n");
}
