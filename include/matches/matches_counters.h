/*
 * PROJECT: GEMMapper
 * FILE: matches_counters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
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
matches_counters_t* matches_counters_new();
void matches_counters_clear(matches_counters_t* const counters);
void matches_counters_delete(matches_counters_t* const counters);

/*
 * Counters
 */
uint64_t matches_counters_get_num_counters(matches_counters_t* const counters);
uint64_t* matches_counters_get_counts(matches_counters_t* const counters);

uint64_t matches_counters_get_count(matches_counters_t* const counters,const uint64_t distance);
uint64_t matches_counters_get_total_count(matches_counters_t* const counters);

void matches_counters_add(
    matches_counters_t* const counters,
    const uint64_t distance,
    const uint64_t num_matches);
void matches_counters_sub(
    matches_counters_t* const counters,
    const uint64_t distance,
    const uint64_t num_matches);

/*
 * Utils
 */
uint64_t matches_counters_compact(matches_counters_t* const counters);
void matches_counters_compute_matches_to_decode(
    matches_counters_t* const counters,
    const uint64_t min_reported_strata,
    const uint64_t min_reported_matches,
    const uint64_t max_reported_matches,
    uint64_t* const reported_strata,
    uint64_t* const last_stratum_reported_matches);

/*
 * Display
 */
void matches_counters_print(
    FILE* const stream,
    matches_counters_t* const matches_counter,
    const uint64_t mcs);

#endif /* MATCHES_COUNTERS_H_ */
