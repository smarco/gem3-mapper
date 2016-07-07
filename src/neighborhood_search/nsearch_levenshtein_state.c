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
 * Prepare
 */
void nsearch_levenshtein_state_prepare(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool supercondensed) {
  // Parameters
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  dp_column_t* const columns = dp_matrix->columns;
  const uint64_t num_columns = dp_matrix->num_columns;
  const uint64_t num_rows = dp_matrix->num_rows;
  // Init
  nsearch_state->supercondensed = supercondensed;
  // Init first row
  uint64_t i;
  columns[0].cells[0] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(0),1);
  if (supercondensed) {
    for (i=1;i<num_columns;++i) {
      columns[i].cells[0] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(0),0);
    }
  } else {
    for (i=1;i<num_columns;++i) {
      columns[i].cells[0] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(i),1);
    }
  }
  // Init first column
  for (i=1;i<num_rows;++i) {
    columns[0].cells[i] = NS_DISTANCE_SET_PRIORITY(NS_DISTANCE_ENCODE(i),1);
  }
}
/*
 * Accessors
 */
uint64_t nsearch_levenshtein_state_get_global_align_distance(
    nsearch_levenshtein_state_t* const nsearch_state,
    const uint64_t key_length,
    const uint64_t text_length,
    const uint64_t max_error) {
  // Check band limits
  if (text_length < key_length-max_error) return NS_DISTANCE_INF;
  // Fetch aling-distance
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  const uint64_t align_distance = dp_matrix->columns[text_length].cells[key_length];
  return NS_DISTANCE_HAS_PRIORITY(align_distance,1) ? NS_DISTANCE_DECODE(align_distance) : NS_DISTANCE_INF;
}
uint64_t nsearch_levenshtein_state_get_local_align_distance(
    nsearch_levenshtein_state_t* const nsearch_state,
    const uint64_t local_key_length,
    const uint64_t global_key_length,
    const uint64_t global_text_length,
    const uint64_t max_error) {
  // Check band limits
  if (global_text_length < global_key_length-max_error) return NS_DISTANCE_INF;
  // Check align-distance
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  const uint64_t align_distance =
      NS_DISTANCE_DECODE(dp_matrix->columns[global_text_length].cells[global_key_length]);
  // Backtrace DP until the local-key border is surpassed
  int64_t h = global_text_length, v = global_key_length;
  int64_t v_left = local_key_length;
  while (h > 0 && v_left > 0) {
    // Fetch columns
    dp_column_t* const current_column = dp_matrix->columns + h;
    dp_column_t* const prev_column = current_column - 1;
    // Fetch cells
    const uint64_t current_cell = current_column->cells[v];
    const uint64_t sub = prev_column->cells[v-1];
    const uint64_t del = current_column->cells[v-1] + 2;
    const uint64_t ins = prev_column->cells[v] + 2;
    // Backtrace
    if (current_cell == del) {
      --v; --v_left;
    } else if (current_cell == ins) {
      --h;
    } else if (current_cell == sub || current_cell == sub+2) {
      --h; --v; --v_left;
    } else {
      GEM_INVALID_CASE();
    }
  }
  // Check current cell
  dp_column_t* const current_column = dp_matrix->columns + h;
  return align_distance - NS_DISTANCE_DECODE(current_column->cells[v]);
}
/*
 * Compute DP
 */
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
/*
 * Compute DP-Banded
 */
void nsearch_levenshtein_state_compute_text_banded(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error,
    uint64_t* const min_align_distance,
    uint64_t* const min_align_distance_column) {
  uint64_t i, align_distance, dummy;
  *min_align_distance = NS_DISTANCE_INF;
  for (i=0;i<text_length;++i) {
    const uint8_t text_char_enc = text[i];
    nsearch_levenshtein_state_compute_chararacter_banded(nsearch_state,
        forward_search,key,key_length,i,text_char_enc,max_error,&dummy,&align_distance);
    if (align_distance < *min_align_distance) {
      *min_align_distance = align_distance;
      *min_align_distance_column = i+1;
    }
  }
}
void nsearch_levenshtein_state_compute_chararacter_banded(
    nsearch_levenshtein_state_t* const nsearch_state,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t text_position,
    const uint8_t text_char_enc,
    const uint64_t max_error,
    uint64_t* const min_val,
    uint64_t* const align_distance) {
  // Index column (current_column=text_offset offset by 1)
  dp_matrix_t* const dp_matrix = &nsearch_state->dp_matrix;
  dp_column_t* const base_column = dp_matrix->columns + text_position;
  dp_column_t* const next_column = base_column + 1;
  // Compute banded offsets
  const uint64_t row_limit = key_length+1;
  const uint64_t current_column = text_position+1;
  uint64_t v_banded_begin, v_banded_end, i;
  v_banded_end = MIN(current_column+max_error+1,row_limit);
  if (nsearch_state->supercondensed || current_column < (max_error+1)) {
    v_banded_begin = 1;
  } else {
    v_banded_begin = current_column - max_error;
    if (v_banded_begin-1 > 0) {
      next_column->cells[v_banded_begin-1] = NS_DISTANCE_INF;
    }
  }
  // Fill columns
  uint64_t column_min = NS_DISTANCE_INF;
  for (i=v_banded_begin;i<v_banded_end;++i) {
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
  next_column->cells[v_banded_end] = NS_DISTANCE_INF;
  // Set min/align values
  *min_val = NS_DISTANCE_DECODE(column_min);
  if (v_banded_end == row_limit) {
    const uint64_t column_last = next_column->cells[key_length];
    *align_distance = NS_DISTANCE_HAS_PRIORITY(column_last,1) ?
        NS_DISTANCE_DECODE(column_last) : NS_DISTANCE_INF;
  } else {
    *align_distance = NS_DISTANCE_INF;
  }
}

