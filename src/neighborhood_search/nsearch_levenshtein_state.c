/*
 * PROJECT: GEMMapper
 * FILE: nsearch_levenshtein_state.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_levenshtein_state.h"

/*
 * Setup
 */
void nsearch_levenshtein_state_init(
    nsearch_levenshtein_state_t* const nsearch_levenshtein_state,
    const uint64_t num_rows,
    const uint64_t num_columns,
    mm_stack_t* const mm_stack) {
  // Init dp_matrix
  dp_matrix_t* const dp_matrix = &nsearch_levenshtein_state->dp_matrix;
  dp_matrix->num_rows = num_rows;
  dp_matrix->num_columns = num_columns;
  // Allocate columns
  dp_column_t* const columns = mm_stack_calloc(mm_stack,num_columns,dp_column_t,false);
  dp_matrix->columns = columns;
  uint64_t i;
  for (i=0;i<num_columns;++i) {
    dp_matrix->columns[i].cells = mm_stack_calloc(mm_stack,num_rows,uint64_t,false);
  }
}
/*
 * Prepare DP
 */
void nsearch_levenshtein_state_prepare_full_neighbourhood(
    nsearch_levenshtein_state_t* const nsearch_state) {
  // Parameters
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  dp_column_t* const columns = dp_matrix->columns;
  const uint64_t num_columns = dp_matrix->num_columns;
  const uint64_t num_rows = dp_matrix->num_rows;
  // Init first row
  uint64_t i;
  columns[0].cells[0] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(0),1);
  for (i=1;i<num_columns;++i) {
    columns[i].cells[0] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(i),1);
  }
  // Init first column
  for (i=1;i<num_rows;++i) {
    columns[0].cells[i] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(i),1);
  }
}
void nsearch_levenshtein_state_prepare_supercondensed_neighbourhood(
    nsearch_levenshtein_state_t* const nsearch_state) {
  // Parameters
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  dp_column_t* const columns = dp_matrix->columns;
  const uint64_t num_columns = dp_matrix->num_columns;
  const uint64_t num_rows = dp_matrix->num_rows;
  // Init first row
  uint64_t i;
  columns[0].cells[0] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(0),1);
  for (i=1;i<num_columns;++i) {
    columns[i].cells[0] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(0),0);
  }
  // Init first column
  for (i=1;i<num_rows;++i) {
    columns[0].cells[i] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(i),1);
  }
}
void nsearch_levenshtein_state_prepare_chained(
    nsearch_levenshtein_state_t* const current_nsearch_state,
    nsearch_levenshtein_state_t* const next_nsearch_state,
    const uint64_t current_key_length,
    const uint64_t current_text_length,
    const bool next_operation_towards_end) {
  // Parameters
  dp_matrix_t* const current_dp_matrix = &current_nsearch_state->dp_matrix;
  dp_column_t* const current_columns = current_dp_matrix->columns;
  dp_matrix_t* const next_dp_matrix = &next_nsearch_state->dp_matrix;
  dp_column_t* const next_columns = next_dp_matrix->columns;
  // Init first column :: Copy first column
  const uint64_t num_rows = next_dp_matrix->num_rows;
  const uint64_t num_columns = next_dp_matrix->num_columns;
  uint64_t i, cell_value;
  for (i=0;i<=current_key_length;++i) {
    next_columns[0].cells[i] = current_columns[current_text_length].cells[i];
  }
  // Init first column :: All
  cell_value = current_columns[current_text_length].cells[current_key_length];
  for (;i<=num_rows;++i) {
    cell_value = NS_DISTANCE_ENCODED_ADD(cell_value,1); // Next
    next_columns[0].cells[i] = cell_value;
  }
  // Init first row :: Remaining cells
  if (next_operation_towards_end) {
    for (i=1;i<num_columns;++i) {
      next_columns[i].cells[0] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(0),0);
    }
  } else {
    cell_value = current_columns[current_text_length].cells[0];
    for (i=1;i<num_columns;++i) {
      cell_value = NS_DISTANCE_ENCODED_ADD(cell_value,1); // Next
      next_columns[i].cells[0] = cell_value;
    }
  }
}
/*
 * Accessors
 */
uint64_t nsearch_levenshtein_state_get_align_distance(
    nsearch_levenshtein_state_t* const nsearch_state,
    const uint64_t key_length,
    const uint64_t text_length) {
  // Parameters
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  const uint64_t align_distance = dp_matrix->columns[text_length].cells[key_length];
  return NS_DISTANCE_HAS_PRIORITY(align_distance,1) ? NS_DISTANCE_DECODE(align_distance) : NS_DISTANCE_INF;
}
/*
 * Compute DP
 */
void nsearch_levenshtein_state_compute_text(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error) {
  uint64_t i, dummy;
  for (i=0;i<text_length;++i) {
    const uint8_t text_char_enc = text[i];
    nsearch_levenshtein_state_compute_chararacter(nsearch_state,
        forward_search,key,key_length,i,text_char_enc,max_error,&dummy,&dummy);
  }
}
void nsearch_levenshtein_state_compute_chararacter(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t text_position,
    const uint8_t text_char_enc,
    const uint64_t max_error,
    uint64_t* const min_val,
    uint64_t* const align_distance) {
  // Parameters
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  // Index column (current_column=text_offset offset by 1)
  dp_column_t* const base_column = dp_matrix->columns + text_position;
  dp_column_t* const next_column = base_column + 1;
  // Fill columns
  uint64_t column_min = NS_DISTANCE_INF, i;
  for (i=1;i<=key_length;++i) {
    // Compute cells
    const uint64_t key_idx = forward_search ? (i-1) : (key_length-i);
    const uint64_t del = next_column->cells[i-1] + 2;
    const uint64_t sub = base_column->cells[i-1] + (text_char_enc==key[key_idx] ? 0 : 2);
    const uint64_t ins = base_column->cells[i] + 2;
    // Compute min
    const uint64_t min2 = MIN(del,sub);
    const uint64_t min3 = MIN(min2,ins);
    next_column->cells[i] = min3;
    if (NS_DISTANCE_HAS_PRIORITY(min3,1)) column_min = MIN(column_min,min3); // Minimum with LOW priority
  }
  // Set min/align values
  const uint64_t column_last = next_column->cells[key_length];
  *min_val = NS_DISTANCE_DECODE(column_min);
  *align_distance = NS_DISTANCE_HAS_PRIORITY(column_last,1) ? NS_DISTANCE_DECODE(column_last) : NS_DISTANCE_INF;
}

