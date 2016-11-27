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

#ifndef STATS_VECTOR_H_
#define STATS_VECTOR_H_

#include "system/commons.h"
#include "utils/hash.h"
#include "stats/stats_counter.h"

/*
 * Constants
 */
#define STATS_VECTOR_OUT_OF_RANGE UINT64_MAX

/*
 * Stats Vector
 */
typedef enum { STATS_VECTOR_CUSTOMED_RANGE, STATS_VECTOR_STEP_RANGE, STATS_VECTOR_RAW } stats_vector_type;
typedef struct {
  /* Aggregated counter */
  gem_counter_t aggregated_counter; // TODO TODO TODO TODO TODO TODO TODO TODO TODO
  /* Stats Vector */
  stats_vector_type type;
  uint64_t* counters;
  uint64_t num_counters;
  uint64_t max_index; // const uint64_t max_index = (num_counters-1)
  uint64_t max_value; // Range = [min_value,max_value). Maximum value (not inclusive).
  uint64_t step; /* STATS_VECTOR_STEP_RANGE */
  uint64_t* customed_range_values; /* STATS_VECTOR_CUSTOMED_RANGE */
  /* Values out of range */
  uint64_t out_of_range_bucket_size;
  ihash_t* out_values;
} stats_vector_t;
typedef struct {
  /* Stats vector */
  stats_vector_t* stats_vector;
  /* Iteration state */
  uint64_t index;
  uint64_t last_index; // Sliced iteration
  bool eoi; // End-of-iteration
  bool counters_iteration; // Iteration state. Over counters or ihash
  /* Hash State */
  ihash_iterator_t* ihash_iterator;
} stats_vector_iterator_t;

/*
 * Customed Range Stat-Vector. Eg.
 *   @range[] = { 0, 1, 2, 3 };
 *     @stats_vector_t[0] => [0,1)
 *     @stats_vector_t[1] => [1,2)
 *     @stats_vector_t[2] => [2,3)
 */
stats_vector_t* stats_vector_customed_range_new(
    uint64_t* const customed_range_values,const uint64_t num_ranges,
    const uint64_t out_of_range_bucket_size);
/*
 * Step Range Stat-Vector. Eg.
 *   @min_value = 0;
 *   @max_value = 6;
 *   @step = 3
 *     @stats_vector_t[0] => {0,1,2}
 *     @stats_vector_t[1] => {3,4,5}
 *     @stats_vector_t[2] => {6,7,8}
 */
stats_vector_t* stats_vector_step_range_new(
    const uint64_t max_value,
    const uint64_t step,
    const uint64_t out_of_range_bucket_size);
/*
 * Raw Stat-Vector. Eg.
 *   @num_values = 3
 *     @stats_vector_t[0] => {0}
 *     @stats_vector_t[1] => {1}
 *     @stats_vector_t[2] => {2}
 */
stats_vector_t* stats_vector_raw_new(
    const uint64_t num_values,
    const uint64_t out_of_range_bucket_size);
/*
 * Create from template copy (just structure, not actual data)
 */
stats_vector_t* stats_vector_new_from_template(stats_vector_t* const stats_vector_template);

/* Setup */
void stats_vector_clear(stats_vector_t* const stats_vector);
void stats_vector_delete(stats_vector_t* const stats_vector);

/*
 * Increment/Add bucket counter
 */
void stats_vector_inc(
    stats_vector_t* const stats_vector,
    const uint64_t value);
void stats_vector_add(
    stats_vector_t* const stats_vector,
    const uint64_t value,
    const uint64_t amount);

/*
 * Bucket counters getters (Individual buckets)
 */
uint64_t* stats_vector_get_counter(stats_vector_t* const stats_vector,const uint64_t value);
uint64_t stats_vector_get_count(stats_vector_t* const stats_vector,const uint64_t value);

/*
 * Bucket counters getters (Accumulated ranges)
 */
uint64_t stats_vector_get_accumulated_count(stats_vector_t* const stats_vector);
// TODO uint64_t stats_vector_get_range_accumulated_count(
//    stats_vector_t* const stats_vector,
//    const uint64_t value_from,
//    const uint64_t value_to);

/*
 * Inverse. Given the stats_vector_t index returns the corresponding value/range.
 */
void stats_vector_get_value_range(
    stats_vector_t* const stats_vector,
    const uint64_t index,
    uint64_t* const lo_value,
    uint64_t* const hi_value);

/*
 * Merge 2 stats-vector (adding bucket counting)
 */
// TODO void stats_vector_merge(stats_vector_t* const stats_dst,stats_vector_t* const stats_src);

/*
 * Display (Printers)
 */
void stats_vector_display(
    FILE* const stream,
    stats_vector_t* const stats_vector,
    const bool display_zeros,
    const bool display_percentage,
    void (*print_label)(uint64_t));
void stats_vector_print_ranges(
    FILE* const stream,
    stats_vector_t* const stats_vector);
void stats_vector_print_values(
    FILE* const stream,
    stats_vector_t* const stats_vector,
    const bool display_percentage);

/*
 * Iterator
 */
stats_vector_iterator_t* stats_vector_iterator_new(stats_vector_t* const stats_vector);
stats_vector_iterator_t* stats_vector_iterator_range_new(
    stats_vector_t* const stats_vector,
    const uint64_t value_from,
    const uint64_t value_to);
void stats_vector_iterator_delete(stats_vector_iterator_t* const sv_iterator);

bool stats_vector_iterator_eoi(stats_vector_iterator_t* const sv_iterator);
void stats_vector_iterator_next(stats_vector_iterator_t* const sv_iterator);
uint64_t stats_vector_iterator_get_index(stats_vector_iterator_t* const sv_iterator);
uint64_t stats_vector_iterator_get_count(stats_vector_iterator_t* const sv_iterator);
void stats_vector_iterator_get_range(
    stats_vector_iterator_t* const sv_iterator,
    uint64_t* const lo_value,
    uint64_t* const hi_value);

#endif /* STATS_VECTOR_H_ */
