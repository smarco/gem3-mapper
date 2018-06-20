/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
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
    mm_allocator_t* const mm_allocator) {
  // Compute dimensions
  const uint64_t num_rows = max_key_length + 1;     // (+1) Initial Conditions Row
  const uint64_t num_columns = max_text_length + 2; // (+1) Initial Conditions Column
                                                    // (+1) Auxiliary Column (check there is no valid match)
  // Allocate text
  nsearch_operation->text = mm_allocator_calloc(mm_allocator,num_columns,uint8_t,false);
  // Init levenshtein state
  nsearch_levenshtein_state_init(&nsearch_operation->nsearch_state,num_rows,num_columns,mm_allocator);
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
    mm_allocator_t* const mm_allocator) {
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
      &min_align_distance,&min_align_distance_column,mm_allocator);
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
    mm_allocator_t* const mm_allocator) {
  mm_allocator_push_state(mm_allocator);
  const uint64_t text_length = strlen(text);
  const bool forward_search = (nsearch_operation->search_direction==direction_forward);
  // Allocate text
  uint8_t* const enc_text = mm_allocator_calloc(mm_allocator,text_length,uint8_t,true);
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
  mm_allocator_pop_state(mm_allocator);
  return eq;
}
/*
 * Display
 */
void nsearch_operation_print(
    FILE* const stream,
    nsearch_operation_t* const nsearch_operation) {
  fprintf(stream,">>LKey[%"PRIu64",%"PRIu64")/GKey[%"PRIu64",%"PRIu64") Error{%"PRIu64"/%"PRIu64",%"PRIu64"}\n",
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
  fprintf(stream,"  => Local-Key [%"PRIu64",%"PRIu64")\n",nsearch_operation->local_key_begin,nsearch_operation->local_key_end);
  fprintf(stream,"  => Global-Key [%"PRIu64",%"PRIu64")\n",nsearch_operation->global_key_begin,nsearch_operation->global_key_end);
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
