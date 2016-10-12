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

#include "utils/interval_set.h"

/*
 * Constants
 */
#define INTERVAL_SET_NUN_INITIAL_INTERVALS 50

/*
 * Setup
 */
void interval_set_init(interval_set_t* const interval_set) {
  interval_set->intervals = vector_new(INTERVAL_SET_NUN_INITIAL_INTERVALS,interval_t);
}
void interval_set_clear(interval_set_t* const interval_set) {
  vector_clear(interval_set->intervals);
}
void interval_set_destroy(interval_set_t* const interval_set) {
  vector_delete(interval_set->intervals);
}
/*
 * Counting
 */
uint64_t interval_set_count_intervals(interval_set_t* const interval_set) {
  return vector_get_used(interval_set->intervals);
}
uint64_t interval_set_count_intervals_length(interval_set_t* const interval_set) {
  uint64_t count = 0;
  INTERVAL_SET_ITERATE(interval_set,interval) {
    count += interval->hi - interval->lo;
  }
  return count;
}
uint64_t interval_set_count_intervals_length_thresholded(
    interval_set_t* const interval_set,
    const uint64_t max_error) {
  uint64_t count = 0;
  INTERVAL_SET_ITERATE(interval_set,interval) {
    if (interval->distance <= max_error) count += interval->hi - interval->lo;
  }
  return count;
}
/*
 * Adding
 */
void interval_set_add(
    interval_set_t* const interval_set,
    const uint64_t lo,
    const uint64_t hi,
    const uint64_t distance,
    const uint64_t length) {
  // Allocate
  interval_t* interval;
  vector_alloc_new(interval_set->intervals,interval_t,interval);
  // Add
  interval->lo = lo;
  interval->hi = hi;
  interval->distance = distance;
  interval->length = length;
}
/*
 * Set Operators
 */
void interval_set_union(
    interval_set_t* const interval_set_a,
    interval_set_t* const interval_set_b) {
  // Appends to @interval_set_a the intervals contained into @interval_set_b (union set)
  const uint64_t total_size = vector_get_used(interval_set_a->intervals) + vector_get_used(interval_set_b->intervals);
  vector_reserve(interval_set_a->intervals,total_size,false);
  interval_t* int_set_a = vector_get_free_elm(interval_set_a->intervals,interval_t);
  INTERVAL_SET_ITERATE(interval_set_b,int_set_b) {
    int_set_a->lo = int_set_b->lo;
    int_set_a->hi = int_set_b->hi;
    int_set_a->distance = int_set_b->distance;
    int_set_a->length = int_set_b->length;
    ++int_set_a;
  }
  vector_set_used(interval_set_a->intervals,total_size);
}
void interval_set_subtract(
    interval_set_t* const result_set,
    interval_set_t* const exclusion_set) {
  // Subtracts to @result_set the intervals contained in @exclusion_set (difference set)
  const uint64_t exclusion_set_size = vector_get_used(exclusion_set->intervals);
  uint64_t result_set_size = vector_get_used(result_set->intervals);
  uint64_t i, j;
  for (i=0;i<result_set_size;++i) {
    interval_t* int_res = vector_get_mem(result_set->intervals,interval_t) + i;
    interval_t* int_excl = vector_get_mem(exclusion_set->intervals,interval_t);
    for (j=0;j<exclusion_set_size;++j,++int_excl) {
      const uint64_t hi1 = int_res->hi;
      const uint64_t lo1 = int_res->lo;
      const uint64_t hi2 = int_excl->hi;
      const uint64_t lo2 = int_excl->lo;
      if (hi1 <= lo2 || hi2 <= lo1) { // Disjoint intervals
        continue;
      } else {
        if (lo2 <= lo1 && hi1 <= hi2) { // Full subtraction
          int_res->lo = int_res->hi; // Close the interval
        } else if (lo1 < lo2 && hi2 < hi1) { // Exclusion inside result
          // Add the chunk to the end
          vector_reserve(result_set->intervals,result_set_size+1,false);
          interval_t* const int_res_chunk = vector_get_mem(result_set->intervals,interval_t) + result_set_size;
          int_res_chunk->lo = hi2;
          int_res_chunk->hi = hi1;
          int_res_chunk->distance = int_res->distance;
          int_res_chunk->length = int_res->length;
          // Shrink the interval
          int_res->hi = lo2;
          ++result_set_size;
        } else if (lo2 == lo1 && lo1 < hi2) { // Exclusion overlaps left side of result
          int_res->lo = hi2;
        } else /* if (lo2 < hi1 && hi1 == hi2) */ { // Exclusion overlaps right side of result
          int_res->hi = lo2;
        }
      }
    }
  }
  vector_set_used(result_set->intervals,result_set_size);
}
