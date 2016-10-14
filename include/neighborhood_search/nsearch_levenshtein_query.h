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

#ifndef NSEARCH_LEVENSHTEIN_QUERY_H_
#define NSEARCH_LEVENSHTEIN_QUERY_H_

#include "utils/essentials.h"
#include "fm_index/fm_index_query.h"
#include "neighborhood_search/nsearch_schedule.h"

/*
 * NSearch query helper
 */
typedef struct {
  fm_2interval_t fm_2interval;
  fm_2erank_elms_t lo_2erank_elms;
  fm_2erank_elms_t hi_2erank_elms;
  uint64_t num_optimization_steps;
} nsearch_query_t;

/*
 * Setup
 */
void nsearch_query_init(
    nsearch_query_t* const nsearch_query,
    fm_index_t* const fm_index);

/*
 * Standard search query
 */
void nsearch_levenshtein_query(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const uint64_t current_position,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out);

/*
 * Scheduled search precomputed-query
 */
void nsearch_levenshtein_scheduled_precompute_query(
    nsearch_schedule_t* const nsearch_schedule,
    const bool forward_search,
    nsearch_query_t* const nsearch_query);

/*
 * Scheduled search query
 */
uint64_t nsearch_levenshtein_scheduled_query_exact(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const int64_t current_position,
    const uint8_t char_enc,
    nsearch_query_t* const nsearch_query);
uint64_t nsearch_levenshtein_scheduled_query(
    nsearch_schedule_t* const nsearch_schedule,
    nsearch_operation_t* const nsearch_operation,
    const bool forward_search,
    const int64_t current_position,
    const uint8_t char_enc,
    nsearch_query_t* const nsearch_query,
    nsearch_query_t* const next_nsearch_query);

#endif /* NSEARCH_LEVENSHTEIN_QUERY_H_ */
