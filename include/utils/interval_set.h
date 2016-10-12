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

#ifndef INTERVAL_SET_H_
#define INTERVAL_SET_H_

#include "utils/essentials.h"

typedef struct {
  /* BWT Interval */
  uint64_t lo;
  uint64_t hi;
  /* Matching Properties */
  uint64_t distance;
  uint64_t length;
} interval_t;

typedef struct {
  vector_t* intervals; // (interval_t*)
} interval_set_t;

/*
 * Setup
 */
void interval_set_init(interval_set_t* const interval_set);
void interval_set_clear(interval_set_t* const interval_set);
void interval_set_destroy(interval_set_t* const interval_set);

/*
 * Counting
 */
uint64_t interval_set_count_intervals(interval_set_t* const interval_set);
uint64_t interval_set_count_intervals_length(interval_set_t* const interval_set);
uint64_t interval_set_count_intervals_length_thresholded(
    interval_set_t* const interval_set,
    const uint64_t max_error);

/*
 * Adding
 */
void interval_set_add(
    interval_set_t* const interval_set,
    const uint64_t lo,
    const uint64_t hi,
    const uint64_t distance,
    const uint64_t length);

/*
 * Set Operators
 */
void interval_set_union(
    interval_set_t* const interval_set_a,
    interval_set_t* const interval_set_b);
void interval_set_subtract(
    interval_set_t* const result_set,
    interval_set_t* const exclusion_set);

/*
 * Macro iterator
 */
#define INTERVAL_SET_ITERATE(interval_set,element) \
  VECTOR_ITERATE(interval_set->intervals,element,interval_set_##position,interval_t)


#endif /* INTERVAL_SET_H_ */
