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

#ifndef STATS_MATRIX_H_
#define STATS_MATRIX_H_

#include "utils/essentials.h"
#include "stats/stats_vector.h"

/*
 * Stats Matrix f(x,y) := (StatsVector x StatsVector)
 */
typedef struct {
  stats_vector_t* dimension_x; // x dimension
  stats_vector_t* dimension_y; // Template for the y dimension
} stats_matrix_t;

/*
 * Setup
 */
stats_matrix_t* stats_matrix_new(
    stats_vector_t* const dimension_x,
    stats_vector_t* const dimension_y);
void stats_matrix_clear(stats_matrix_t* const stats_matrix);
void stats_matrix_delete(stats_matrix_t* const stats_matrix);

/*
 * Increment/Add bucket counter
 */
void stats_matrix_inc(
    stats_matrix_t* const stats_matrix,
    const uint64_t value_x,
    const uint64_t value_y);
void stats_matrix_add(
    stats_matrix_t* const stats_matrix,
    const uint64_t value_x,
    const uint64_t value_y,
    const uint64_t amount);

/*
 * Bucket counters getters (Individual buckets)
 */
uint64_t stats_matrix_get_count(
    stats_matrix_t* const stats_matrix,
    const uint64_t value_x,
    const uint64_t value_y);

/*
 * Display (Printers)
 */
void stats_matrix_display(
    FILE* const stream,
    stats_matrix_t* const stats_matrix,
    const bool display_percentage,
    void (*print_label)(uint64_t));

#endif /* STATS_MATRIX_H_ */
