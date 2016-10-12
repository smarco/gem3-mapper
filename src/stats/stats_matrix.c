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

#include "stats/stats_matrix.h"

/*
 * Setup
 */
stats_matrix_t* stats_matrix_new(
    stats_vector_t* const dimension_x,
    stats_vector_t* const dimension_y) {
  // Allocate handler
  stats_matrix_t* const stats_matrix = mm_alloc(stats_matrix_t);
  // Init Dimensions
  stats_matrix->dimension_x = dimension_x;
  stats_matrix->dimension_y = dimension_y;
  // Return
  return stats_matrix;
}
void stats_matrix_clear(stats_matrix_t* const stats_matrix) {
  stats_vector_iterator_t* const iterator_x = stats_vector_iterator_new(stats_matrix->dimension_x);
  while (!stats_vector_iterator_eoi(iterator_x)) {
    const uint64_t counter = stats_vector_iterator_get_count(iterator_x);
    if (counter > 0) stats_vector_clear((stats_vector_t*)counter);
    // Next
    stats_vector_iterator_next(iterator_x);
  }
  stats_vector_iterator_delete(iterator_x);
}
void stats_matrix_delete(stats_matrix_t* const stats_matrix) {
  stats_vector_iterator_t* const iterator_x = stats_vector_iterator_new(stats_matrix->dimension_x);
  while (!stats_vector_iterator_eoi(iterator_x)) {
    const uint64_t counter = stats_vector_iterator_get_count(iterator_x);
    if (counter > 0) stats_vector_delete((stats_vector_t*)counter);
    // Next
    stats_vector_iterator_next(iterator_x);
  }
  stats_vector_iterator_delete(iterator_x);
  // Delete X/Y
  stats_vector_delete(stats_matrix->dimension_x);
  stats_vector_delete(stats_matrix->dimension_y);
}
/*
 * Increment/Add bucket counter
 */
stats_vector_t* stats_matrix_get_y_dimension(stats_matrix_t* const stats_matrix,const uint64_t value_x) {
  // Index X-Dimension and get Y-Stats_vector
  uint64_t* const dimension_y_placeholder = stats_vector_get_counter(stats_matrix->dimension_x,value_x);
  // Use the pointer to Y stored as uint64_t
  if (*dimension_y_placeholder == 0) {
    *dimension_y_placeholder = (uint64_t) stats_vector_new_from_template(stats_matrix->dimension_y);
    ((stats_vector_t*)(*dimension_y_placeholder))->out_of_range_bucket_size = UINT64_MAX;
  }
  return (stats_vector_t*) (*dimension_y_placeholder);
}
void stats_matrix_inc(
    stats_matrix_t* const stats_matrix,
    const uint64_t value_x,
    const uint64_t value_y) {
  // Increment the x-Dimension
  stats_vector_inc(stats_matrix_get_y_dimension(stats_matrix,value_x),value_y);
}
void stats_matrix_add(
    stats_matrix_t* const stats_matrix,
    const uint64_t value_x,
    const uint64_t value_y,
    const uint64_t amount) {
  // Add the x-Dimension
  stats_vector_add(stats_matrix_get_y_dimension(stats_matrix,value_x),value_y,amount);
}
/*
 * Bucket counters getters (Individual buckets)
 */
uint64_t stats_matrix_get_count(
    stats_matrix_t* const stats_matrix,
    const uint64_t value_x,
    const uint64_t value_y) {
  return stats_vector_get_count(stats_matrix_get_y_dimension(stats_matrix,value_x),value_y);
}
/*
 * Display (Printers)
 */
void stats_matrix_display(
    FILE* const stream,
    stats_matrix_t* const stats_matrix,
    const bool display_percentage,
    void (*print_label)(uint64_t)) {
  // Print Y-Labels
  fprintf(stream,"\t");
  stats_vector_print_ranges(stream,stats_matrix->dimension_y);
  // Print X-Label & (X,(Y-Values))
  stats_vector_iterator_t* const iterator = stats_vector_iterator_new(stats_matrix->dimension_x);
  while (!stats_vector_iterator_eoi(iterator)) {
    // Get current counter
    uint64_t lo_range, hi_range;
    const uint64_t counter = stats_vector_iterator_get_count(iterator);
    if (counter > 0) {
      // Print Label
      if (print_label==NULL) {
        stats_vector_iterator_get_range(iterator,&lo_range,&hi_range);
        if (hi_range-lo_range>1) {
          fprintf(stream,"[%"PRIu64",%"PRIu64") => ",lo_range,hi_range);
        } else {
          fprintf(stream,"[%"PRIu64"] => ",lo_range);
        }
      } else {
        print_label(stats_vector_iterator_get_index(iterator)); // Print label
        fprintf(stream," => ");
      }
      // Print Y-Dimension
      stats_vector_t* const stats_vector = (stats_vector_t*)counter;
      stats_vector_print_values(stream,stats_vector,display_percentage);
    }
    // Next
    stats_vector_iterator_next(iterator);
  }
  stats_vector_iterator_delete(iterator);
}
