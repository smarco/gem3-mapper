/*
 * PROJECT: GEMMapper
 * FILE: matches_counters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_COUNTERS_H_
#define MATCHES_COUNTERS_H_

#include "essentials.h"

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
GEM_INLINE matches_counters_t* matches_counters_new();
GEM_INLINE void matches_counters_clear(matches_counters_t* const counters);
GEM_INLINE void matches_counters_delete(matches_counters_t* const counters);

/*
 * Counters
 */
GEM_INLINE uint64_t matches_counters_get_num_counters(matches_counters_t* const counters);
GEM_INLINE uint64_t* matches_counters_get_counts(matches_counters_t* const counters);

GEM_INLINE uint64_t matches_counters_get_count(matches_counters_t* const counters,const uint64_t distance);
GEM_INLINE uint64_t matches_counters_get_total_count(matches_counters_t* const counters);

GEM_INLINE void matches_counters_add(matches_counters_t* const counters,const uint64_t distance,const uint64_t num_matches);
GEM_INLINE void matches_counters_sub(matches_counters_t* const counters,const uint64_t distance,const uint64_t num_matches);

/*
 * Utils
 */
GEM_INLINE uint64_t matches_counters_compact(matches_counters_t* const counters);
GEM_INLINE void matches_counters_compute_matches_to_decode(
    matches_counters_t* const counters,const uint64_t min_decoded_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches,
    uint64_t* const strata_to_decode,uint64_t* const matches_to_decode_from_last_stratum);

#endif /* MATCHES_COUNTERS_H_ */
