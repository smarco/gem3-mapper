/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein_state.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_levenshtein_state.h"
#include "data_structures/dna_text.h"

/*
 * DP-Matrix
 */
void nsearch_levenshtein_state_init(
    nsearch_levenshtein_state_t* const nsearch_levenshtein_state,
    const uint64_t num_columns,const uint64_t column_length,
    mm_stack_t* const mm_stack) {
  // Init dp_matrix
  dp_matrix_t* const dp_matrix = &nsearch_levenshtein_state->dp_matrix;
  dp_matrix->column_length = column_length;
  dp_matrix->num_columns = num_columns;
  // Allocate columns
  dp_column_t* const columns = mm_stack_calloc(mm_stack,num_columns,dp_column_t,false);
  dp_matrix->columns = columns;
  uint64_t i;
  for (i=0;i<num_columns;++i) {
    dp_matrix->columns[i].cells = mm_stack_calloc(mm_stack,column_length,uint64_t,false);
  }
  // Allocate text
  nsearch_levenshtein_state->local_text = mm_stack_calloc(mm_stack,num_columns,uint8_t,false);
  nsearch_levenshtein_state->global_text = mm_stack_calloc(mm_stack,num_columns,uint8_t,false);
}
/*
 * Prepare DP
 */
void nsearch_levenshtein_state_prepare_full(
    nsearch_levenshtein_state_t* const nsearch_state,const uint64_t key_begin,
    const uint64_t key_end,const uint64_t max_error) {
  // Parameters
  const uint64_t key_length = key_end - key_begin;
  const uint64_t max_text_length = key_length + max_error;
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  dp_column_t* const columns = dp_matrix->columns;
  // Initialize columns
  const uint64_t column_start_band = max_error+1;
  uint64_t i;
  columns[0].cells[0] = NS_SET_PRIORITY(NS_ENCODE_DISTANCE(0),1);
  for (i=1;i<=key_length;++i) columns[0].cells[i] = NS_SET_PRIORITY(NS_ENCODE_DISTANCE(i),1);
  for (i=0;i<=max_text_length;++i) {
    uint64_t band_high_offset = max_error + i + 1;
    dp_column_t* base_column = columns + i;
    dp_column_t* next_column = base_column + 1;
    next_column->band_high_offset = MIN(band_high_offset,key_length);
    if (i < column_start_band) {
      next_column->band_low_offset = 1;
      next_column->cells[0] = NS_SET_PRIORITY(NS_ENCODE_DISTANCE(i+1),1);
    } else {
      next_column->band_low_offset = i - column_start_band + 1;
      next_column->cells[next_column->band_low_offset-1] = NS_INF;
    }
    base_column->cells[next_column->band_high_offset] = NS_INF;
    next_column->cells[key_length] = NS_INF;
  }
}
void nsearch_levenshtein_state_prepare_supercondensed(
    nsearch_levenshtein_state_t* const nsearch_state,const uint64_t max_error) {
  // Parameters
  const uint64_t key_length = nsearch_state->key_end - nsearch_state->key_begin;
  const uint64_t max_text_length = key_length + max_error;
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  dp_column_t* const columns = dp_matrix->columns;
  // Initialize columns
  const uint64_t column_start_band = max_error+1;
  uint64_t i;
  columns[0].cells[0] = NS_SET_PRIORITY(NS_ENCODE_DISTANCE(0),1);
  for (i=1;i<=key_length;++i) columns[0].cells[i] = NS_SET_PRIORITY(NS_ENCODE_DISTANCE(i),1);
  for (i=0;i<=max_text_length;++i) {
    uint64_t band_high_offset = max_error + i + 1;
    dp_column_t* base_column = columns + i;
    dp_column_t* next_column = base_column + 1;
    next_column->band_high_offset = MIN(band_high_offset,key_length);
    if (i < column_start_band) {
      next_column->band_low_offset = 1;
      next_column->cells[0] = NS_SET_PRIORITY(NS_ENCODE_DISTANCE(0),0);
    } else {
      next_column->band_low_offset = i - column_start_band + 1;
      next_column->cells[next_column->band_low_offset-1] = NS_INF;
    }
    base_column->cells[next_column->band_high_offset] = NS_INF;
    next_column->cells[key_length] = NS_INF;
  }
}
/*
 * Compute DP
 */
void nsearch_levenshtein_state_compute_chararacter(
    nsearch_levenshtein_state_t* const nsearch_state,const bool forward_search,
    const uint8_t* const key,const uint64_t key_begin,const uint64_t key_end,
    const uint64_t text_offset,const uint8_t text_char_enc,
    uint64_t* const min_val,uint64_t* const align_distance) {
  // Parameters
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  // Index column (current_column=text_offset offset by 1)
  dp_column_t* const base_column = dp_matrix->columns + text_offset;
  dp_column_t* const next_column = base_column + 1;
  // Fill columns
  const uint64_t band_low_offset = next_column->band_low_offset;
  const uint64_t band_high_offset = next_column->band_high_offset;
  uint64_t column_min = NS_INF, i;
  for (i=band_low_offset;i<=band_high_offset;++i) {
    // Compute cells
    const uint64_t key_idx = forward_search ? key_begin+(i-1) : key_end-i;
    const uint64_t del = next_column->cells[i-1] + 2;
    const uint64_t sub = base_column->cells[i-1] + (text_char_enc==key[key_idx] ? 0 : 2);
    const uint64_t ins = base_column->cells[i] + 2;
    // Compute min
    const uint64_t min2 = MIN(del,sub);
    const uint64_t min3 = MIN(min2,ins);
    next_column->cells[i] = min3;
    if (NS_HAS_PRIORITY(min3,1)) column_min = MIN(column_min,min3);
  }
  // Set min/align values
  const uint64_t column_last = next_column->cells[key_end-key_begin];
  *min_val = NS_DECODE_DISTANCE(column_min);
  *align_distance = NS_HAS_PRIORITY(column_last,1) ? NS_DECODE_DISTANCE(column_last) : NS_INF;
}
void nsearch_levenshtein_state_compute_sequence(
    nsearch_levenshtein_state_t* const nsearch_state,
    nsearch_levenshtein_state_t* const next_nsearch_state,
    const bool forward_search) {
  const uint64_t local_text_length = nsearch_state->local_text_length;
  uint8_t* const local_text = nsearch_state->local_text;
  const uint64_t global_text_length = nsearch_state->global_text_length;
  uint8_t* const global_text = nsearch_state->global_text;
  uint8_t* const next_global_text = next_nsearch_state->global_text;
  uint64_t i;
  if (forward_search) {
    for (i=0;i<global_text_length;++i) next_global_text[i] = global_text[i];
    for (i=0;i<local_text_length;++i) next_global_text[global_text_length+i] = local_text[i];
  } else {
    const uint64_t local_last_idx = local_text_length-1;
    for (i=0;i<local_text_length;++i) next_global_text[i] = local_text[local_last_idx-i];
    const uint64_t global_last_idx = global_text_length-1;
    for (i=0;i<global_text_length;++i) next_global_text[local_text_length+i] = global_text[global_last_idx-i];
  }
  next_nsearch_state->local_text_length = 0;
  next_nsearch_state->global_text_length = global_text_length+local_text_length;
}
/*
 * Display
 */
void nsearch_levenshtein_state_print(
    FILE* const stream,nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,const uint8_t* const key) {
  fprintf(stream,"[GEM]> Levenshtein.State\n");
  fprintf(stream,"  => Search %s\n",forward_search ? "forward" : "reverse");
  fprintf(stream,"  => Key [%lu,%lu)\n",nsearch_state->key_begin,nsearch_state->key_end);
  fprintf(stream,"  => Global-Text ");
  nsearch_levenshtein_state_print_search_text(stream,nsearch_state,forward_search);
  fprintf(stream,"  => Local-Text  ");
  nsearch_levenshtein_state_print_local_text(stream,nsearch_state,forward_search);
  fprintf(stream,"  => DP-Matrix\n");
  dp_matrix_print(stream,&nsearch_state->dp_matrix,forward_search,
      key,nsearch_state->key_begin,nsearch_state->key_end,
      nsearch_state->local_text,0,nsearch_state->local_text_length);
}
void nsearch_levenshtein_state_print_search_text(
    FILE* const stream,nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search) {
  if (forward_search) {
    dna_buffer_print(stream,nsearch_state->global_text,nsearch_state->global_text_length,false);
    dna_buffer_print(stream,nsearch_state->local_text,nsearch_state->local_text_length,false);
  } else {
    dna_buffer_print(stream,nsearch_state->local_text,nsearch_state->local_text_length,true);
    dna_buffer_print(stream,nsearch_state->global_text,nsearch_state->global_text_length,true);
  }
  fprintf(stream,"\n");
}
void nsearch_levenshtein_state_print_local_text(
    FILE* const stream,nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search) {
  const uint8_t* const local_text = nsearch_state->local_text;
  const uint64_t local_text_length = nsearch_state->local_text_length;
  dna_buffer_print(stream,local_text,local_text_length,!forward_search);
  fprintf(stream,"\n");
}


