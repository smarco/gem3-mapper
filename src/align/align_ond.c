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
 * DESCRIPTION:
 *   Alignment module using O(nd)-algorithm to compute levenshtein distance/alignment
 *   (Myers' O(nd)-algorithm to compute levenshtein distance/alignment)
 */

#include "align/alignment.h"
#include "align/align_ond.h"
#include "matches/matches.h"
#include "matches/matches_cigar.h"

/*
 * O(ND) Compute LCS
 */
void align_ond_compute_lcs_distance(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    uint64_t* const lcs_distance,
    uint64_t* const match_end_column,
    mm_allocator_t* const mm_allocator) {
  // Allocation
  mm_allocator_push_state(mm_allocator);
  const int32_t max_lcs = key_length + text_length;
  int32_t* const V = mm_allocator_calloc(mm_allocator,max_lcs,int32_t,true);
  // Init
  int32_t d, k;
  V[1] = 0;
  for (d=0;d<=max_lcs;d++) {
    for (k=-d;k<=d;k+=2) {
      // Select adjacent neighbor (up or left)
      const bool lower_neighbor_line = ((k == -d) || (k != d && V[k-1] < V[k+1]));
      const int32_t neighbor = lower_neighbor_line ? k + 1 : k - 1;
      // Compute snake
      const int32_t begin_column = V[neighbor];
      int32_t end_column = lower_neighbor_line ? begin_column : begin_column + 1;
      int32_t end_row = end_column - k;
      // Follow diagonal (matching)
      while (end_column < text_length && end_row < key_length && text[end_column] == key[end_row]) {
        end_column++;
        end_row++;
      }
      // Store end point
      V[k] = end_column;
      // Check for solution
      if (end_column >= text_length && end_row >= key_length) {
        *lcs_distance = d;
        *match_end_column = end_column;
        mm_allocator_pop_state(mm_allocator);
        return;
      }
    }
  }
  // Not found
  *lcs_distance = ALIGN_DISTANCE_INF;
  *match_end_column = ALIGN_COLUMN_INF;
  mm_allocator_pop_state(mm_allocator);
}
/*
 * O(ND) Align
 */
void align_ond_compute_contours(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    const int32_t max_distance,
    align_ond_contours_t* const align_ond_contours,
    mm_allocator_t* const mm_allocator) {
  // Allocation
  const int32_t num_contours = max_distance+2; // (+1) init (+1) max-distance
  align_ond_contours->contour = mm_allocator_calloc(mm_allocator,num_contours,int32_t*,false);
  int32_t i;
  align_ond_contours->contour[0] = mm_allocator_calloc(mm_allocator,2,int32_t,false);
  for (i=1;i<num_contours;++i) {
    const int32_t d = i-1;
    align_ond_contours->contour[i] = mm_allocator_calloc(mm_allocator,2*d+1,int32_t,false) + d /* offset */;
  }
  // Compute contours
  int32_t d, k;
  align_ond_contours->contour[0][1] = 0; // Init
  for (d=0;d<=max_distance;++d) {
    int32_t* const prev_contour = align_ond_contours->contour[d];
    int32_t* const current_contour = align_ond_contours->contour[d+1];
    for (k=-d;k<=d;k+=2) {
      // Select adjacent neighbor (up or left)
      const bool upper_neighbor_line = ((k == -d) || (k != d && prev_contour[k-1] < prev_contour[k+1]));
      const int32_t neighbor = upper_neighbor_line ? k + 1 : k - 1;
      // Compute snake
      const int32_t begin_column = prev_contour[neighbor];
      int32_t end_column = upper_neighbor_line ? begin_column : begin_column + 1;
      int32_t end_row = end_column - k;
      // Follow diagonal (matching)
      while (end_column < text_length && end_row < key_length && text[end_column] == key[end_row]) {
        end_column++;
        end_row++;
      }
      // Store end point
      current_contour[k] = end_column;
      // Check for solution
      if (end_row >= key_length && end_column >= text_length) { // &&
        align_ond_contours->lcs_distance = d;
        align_ond_contours->match_end_column = end_column;
        return;
      }
    }
    // DEBUG
    align_ond_print_contour(stderr,current_contour,-d,d,d);
  }
  // No solution found (within limits)
  align_ond_contours->lcs_distance = ALIGN_DISTANCE_INF;
  align_ond_contours->match_end_column = ALIGN_COLUMN_INF;
}
void align_ond_backtrace_contours(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    align_ond_contours_t* const align_ond_contours,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector) {
  // Parameters
  int32_t h = align_ond_contours->match_end_column, v = key_length;
  int32_t d = align_ond_contours->lcs_distance;
  // Allocate CIGAR string memory (worst case)
  match_alignment->cigar_offset = vector_get_used(cigar_vector); // Set CIGAR offset
  vector_reserve_additional(cigar_vector,MIN(key_length,2*align_ond_contours->lcs_distance+1)); // Reserve
  cigar_element_t* cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer;
  cigar_buffer->type = cigar_null; // Trick
  // Retrieve the alignment
  int64_t match_effective_length = key_length;
  while (true) {
    // Current contour & position
    int32_t* const prev_contour = align_ond_contours->contour[d];
    const int32_t* const current_contour = align_ond_contours->contour[d+1];
    const int32_t k = h - v;
    // Retrieve adjacent neighbor (up or left)
    const bool upper_neighbor_line = ((k==-d) || (k!=d && prev_contour[k-1] < prev_contour[k+1]));
    const int32_t neighbor = upper_neighbor_line ? k + 1 : k - 1;
    // Retrieve snake points
    const int32_t end_column = current_contour[k];
    const int32_t begin_column = prev_contour[neighbor];
    const int32_t begin_row = begin_column - neighbor;
    const int32_t mid_column = upper_neighbor_line ? begin_column : begin_column + 1;
    // Store CIGAR matches
    if (end_column > mid_column) {
      const uint64_t column_range = end_column-mid_column;
      uint64_t i;
      for (i=1;i<column_range;++i) {
        if (key[v-i] != text[h-i]) {
          fprintf(stderr,"Error\n");
        }
      }
      matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_match,column_range);
    }
    // Update position & check end
    h = begin_column;
    v = begin_row;
    --d;
    if (v <= 0 || h <= 0) {
      if (v > 0) {
        matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_del,v);
        match_effective_length -= v;
      }
      if (h > 0) {
        //match_alignment->match_position += h;
        matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_ins,h);
      }
      break;
    }
    // Store CIGAR Indels
    if (upper_neighbor_line) {
      matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_del,1);
      --match_effective_length;
    } else {
      matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_ins,1);
      ++match_effective_length;
    }
  }
  // Set effective length
  match_alignment->effective_length = match_effective_length;
  // Set CIGAR buffer used
  if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
  const uint64_t num_cigar_elements = cigar_buffer - cigar_buffer_base;
  match_alignment->cigar_length = num_cigar_elements; // Set CIGAR length
  // Reverse CIGAR Elements
  if (num_cigar_elements > 0) {
    const uint64_t middle_point = num_cigar_elements/2;
    uint64_t i;
    for (i=0;i<middle_point;++i) {
      SWAP(cigar_buffer_base[i],cigar_buffer_base[num_cigar_elements-i-1]);
    }
  }
  // Set used
  vector_add_used(cigar_vector,num_cigar_elements);
}
void align_ond_match(
    const uint8_t* const key,
    const int32_t key_length,
    const uint8_t* const text,
    const int32_t text_length,
    const int32_t max_distance,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator) {
  // Compute contours
  mm_allocator_push_state(mm_allocator); // Save allocator state
  align_ond_contours_t align_ond_contours;
  align_ond_compute_contours(
      key,key_length,text,text_length,
      max_distance,&align_ond_contours,mm_allocator);
  // Set distance
  match_alignment->score = align_ond_contours.lcs_distance;
  if (align_ond_contours.lcs_distance == ALIGN_DISTANCE_INF) {
    match_alignment->cigar_length = 0;
    mm_allocator_pop_state(mm_allocator); // Free
    return;
  }
  // Backtrace and generate CIGAR
  align_ond_backtrace_contours(
      key,key_length,text,text_length,
      &align_ond_contours,match_alignment,cigar_vector);
  // Free
  mm_allocator_pop_state(mm_allocator);
}
/*
 * Display
 */
void align_ond_print_contour(
    FILE* const stream,
    const int32_t* const contour,
    const int32_t begin_contour,
    const int32_t end_contour,
    const int32_t distance) {
  int32_t i;
  fprintf(stream,"[%d]",distance);
  for (i=begin_contour;i<=end_contour;++i) fprintf(stream,"\t%d",contour[i]);
  fprintf(stream,"\n");
}
