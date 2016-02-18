/*
 * PROJECT: GEMMapper
 * FILE: dp_matrix.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/dp_matrix.h"
#include "data_structures/dna_text.h"

/*
 * DP-Matrix Traceback
 */
void dp_matrix_traceback(
    FILE* const stream,
    const dp_matrix_t* const dp_matrix,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t column_position) {
  dp_column_t* const columns = dp_matrix->columns;
  int64_t h = column_position;
  int64_t v = key_length;
  tab_fprintf(stream,"");
  while (h > 0 && v > 0) {
    const bool match = dna_encode(text[h-1])==key[key_length-v];
    const uint64_t current = NS_DECODE_DISTANCE(columns[h].cells[v]);
    const uint64_t del = NS_DECODE_DISTANCE(columns[h].cells[v-1]) + 1;
    const uint64_t sub = NS_DECODE_DISTANCE(columns[h-1].cells[v-1]) + (match ? 0 : 1);
    // const uint64_t ins = (dp_columns[h-1][v] & NSC_PROPER_CELL_EXTRACT_MASK) + 1;
    if (current == del) {
      fprintf(stream,"D"); --v;
    } else if (current == sub) {
      if (!match) {
        fprintf(stream,"%c",text[h]);
      } else {
        fprintf(stream,".");
      }
      --v; --h;
    } else {
      fprintf(stream,"I"); --h;
    }
  }
  while (v-- > 0) fprintf(stream,"D");
  while (h-- > 0) fprintf(stream,"I");
  fprintf(stream,"\n");
}
/*
void neighborhood_search_debug_match(
    FILE* const stream,
    const dp_matrix_t* const dp_matrix,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t column_position,
    const uint64_t min_val,
    const uint64_t num_matches_found) {
  dp_column_t* const next_column = dp_matrix->columns + column_position;
  uint64_t i;
  for (i=0;i<=key_length;++i) {
    fprintf(stream,"%"PRId64" ",next_column->cells[i]>1000 ? -1 : (int64_t)next_column->cells[i]);
  }
  fprintf(stream,"\n");
  dp_matrix_traceback(stream,dp_matrix,key,key_length,text,column_position);
  // dp_matrix_print(stream,dp_matrix,10,10);
  fprintf(stream,"\n[%02"PRIu64"](%03"PRIu64")> %.*s\n",
      min_val,num_matches_found,(int)column_position,text);
}*/
/*
 * Display
 */
void dp_column_print_summary(
    FILE* const stream,
    const dp_matrix_t* const dp_matrix,
    const uint64_t column_position,
    const uint64_t lo,
    const uint64_t hi) {
  const uint64_t column_length = dp_matrix->column_length;
  uint64_t* const cells = dp_matrix->columns[column_position].cells;
  uint64_t i, min = UINT64_MAX;
  for (i=0;i<=column_length;++i) min = MIN(min,cells[i]);
  tab_fprintf(stream,">> [%"PRIu64"](#%"PRIu64"){min=%"PRIu64",last=%"PRIu64"}\n",
      column_position,hi-lo,min,cells[column_length]);
}
void dp_matrix_print(
    FILE* const stream,
    const dp_matrix_t* const dp_matrix,
    const bool forward_search,
    const uint8_t* const key,
    const uint64_t key_begin,
    const uint64_t key_end,
    const uint8_t* const text,
    const uint64_t text_begin,
    const uint64_t text_end) {
  int64_t h, v;
  // Print Text
  fprintf(stream,"       ");
  for (h=text_begin;h<text_end;++h) fprintf(stream,"   %c",dna_decode(text[h]));
  fprintf(stream,"\n");
  // Print table
  const uint64_t num_columns = text_end-text_begin;
  const uint64_t num_rows = key_end-key_begin;
  for (v=0;v<=num_rows;++v) {
    if (v > 0) {
      fprintf(stream,"%c   ",forward_search ? dna_decode(key[key_begin+v-1]) : dna_decode(key[key_end-v]));
    } else {
      fprintf(stream,"    ");
    }
    for (h=0;h<=num_columns;++h) {
      const uint64_t cell = dp_matrix->columns[h].cells[v];
      if (cell > 1000) {
        fprintf(stream,"inf ");
      } else {
        fprintf(stream,"%c%02lu ",NS_HAS_PRIORITY(cell,0)? '*' : ' ',NS_DECODE_DISTANCE(cell));
      }
    }
    fprintf(stream,"\n");
  }
}
