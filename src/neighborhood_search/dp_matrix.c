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
#include "neighborhood_search/dp_matrix.h"

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
    const uint32_t current = NS_DISTANCE_DECODE(columns[h].cells[v]);
    const uint32_t del = NS_DISTANCE_DECODE(columns[h].cells[v-1]) + 1;
    const uint32_t sub = NS_DISTANCE_DECODE(columns[h-1].cells[v-1]) + (match ? 0 : 1);
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
 * Display
 */
void dp_column_print_summary(
    FILE* const stream,
    const dp_matrix_t* const dp_matrix,
    const uint64_t column_position,
    const uint64_t lo,
    const uint64_t hi) {
  const uint64_t num_rows = dp_matrix->num_rows;
  uint32_t* const cells = dp_matrix->columns[column_position].cells;
  uint32_t min = NS_DISTANCE_INF;
  uint64_t i;
  for (i=0;i<=num_rows;++i) min = MIN(min,cells[i]);
  tab_fprintf(stream,">> [%"PRIu64"](#%"PRIu64"){min=%"PRIu32",last=%"PRIu32"}\n",
      column_position,hi-lo,min,cells[num_rows]);
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
  for (h=text_begin;h<text_end;++h) {
    fprintf(stream,"   %c",dna_decode(text[h]));
  }
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
      const uint32_t cell = dp_matrix->columns[h].cells[v];
      if (cell > 1000) {
        fprintf(stream,"inf ");
      } else {
        fprintf(stream,"%c%02u ",NS_DISTANCE_HAS_PRIORITY(cell,0)? '*' : ' ',NS_DISTANCE_DECODE(cell));
      }
    }
    fprintf(stream,"\n");
  }
}
