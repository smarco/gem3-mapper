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

#ifndef MATCHES_COUNTERS_H_
#define MATCHES_COUNTERS_H_

#include "utils/essentials.h"

/*
 * Counters
 */
typedef struct {
  vector_t* counts;     // Global counters
  uint64_t total_count; // Sum of all (reflected in the counters)
} matches_counters_t;

/*
 * Setup
 */
matches_counters_t* matches_counters_new(void);
void matches_counters_clear(matches_counters_t* const counters);
void matches_counters_delete(matches_counters_t* const counters);

/*
 * Counters
 */
uint64_t matches_counters_get_num_counters(
    matches_counters_t* const matches_counters);
uint64_t* matches_counters_get_counts(
    matches_counters_t* const matches_counters);

uint64_t matches_counters_get_count(
    matches_counters_t* const matches_counters,
    const uint64_t distance);
uint64_t matches_counters_get_total_count(
    matches_counters_t* const matches_counters);

void matches_counters_add(
    matches_counters_t* const matches_counters,
    const uint64_t distance,
    const uint64_t num_matches);
void matches_counters_sub(
    matches_counters_t* const matches_counters,
    const uint64_t distance,
    const uint64_t num_matches);

/*
 * Utils
 */
uint64_t matches_counters_compact(matches_counters_t* const matches_counters);
void matches_counters_compute_matches_to_report(
    matches_counters_t* const matches_counters,
    const uint64_t min_reported_strata,
    const uint64_t max_reported_matches,
    uint64_t* const matches_to_report,
    uint64_t* const strata_to_report);

uint64_t matches_counters_count_first_subdominant(
    matches_counters_t* const matches_counters);
void matches_counters_count_delta_edit(
    matches_counters_t* const matches_counters,
    int64_t* const best_edit_distance,
    int64_t* const subdominant_edit_distance);

/*
 * Display
 */
void matches_counters_print(
    FILE* const stream,
    matches_counters_t* const matches_counter,
    const uint64_t mcs);

#endif /* MATCHES_COUNTERS_H_ */
