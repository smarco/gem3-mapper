/*
 * PROJECT: GEMMapper
 * FILE: matches_counters.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches_counters.h"

/*
 * Constants
 */
#define MATCHES_INIT_COUNTERS          200

/*
 * Setup
 */
GEM_INLINE matches_counters_t* matches_counters_new() {
  // Allocate handler
  matches_counters_t* const counters = mm_alloc(matches_counters_t);
  // Init Counters
  counters->counts = vector_new(MATCHES_INIT_COUNTERS,uint64_t);
  counters->total_count = 0;
  // Return
  return counters;
}
GEM_INLINE void matches_counters_clear(matches_counters_t* const counters) {
  vector_clear(counters->counts);
  counters->total_count = 0;
}
GEM_INLINE void matches_counters_delete(matches_counters_t* const counters) {
  vector_delete(counters->counts);
  mm_free(counters);
}

/*
 * Counters
 */
GEM_INLINE uint64_t matches_counters_get_num_counters(matches_counters_t* const counters) {
  return vector_get_used(counters->counts);
}
GEM_INLINE uint64_t* matches_counters_get_counts(matches_counters_t* const counters) {
  return vector_get_mem(counters->counts,uint64_t);
}
GEM_INLINE uint64_t matches_counters_get_count(matches_counters_t* const counters,const uint64_t distance) {
  return *vector_get_elm(counters->counts,distance,uint64_t);
}
GEM_INLINE uint64_t matches_counters_get_total_count(matches_counters_t* const counters) {
  return counters->total_count;
}
GEM_INLINE void matches_counters_add(matches_counters_t* const counters,const uint64_t distance,const uint64_t num_matches) {
  vector_t* const counts = counters->counts;
  // Reserve Memory
  if (distance >= vector_get_used(counts)) {
    vector_reserve(counts,distance+1,true);
    vector_set_used(counts,distance+1);
  }
  // Add matches
  *vector_get_elm(counts,distance,uint64_t) += num_matches;
  counters->total_count += num_matches;
}
GEM_INLINE void matches_counters_sub(matches_counters_t* const counters,const uint64_t distance,const uint64_t num_matches) {
  uint64_t* const counts = vector_get_elm(counters->counts,distance,uint64_t);
  *counts -= num_matches;
  counters->total_count -= num_matches;
}
/*
 * Utils
 */
GEM_INLINE uint64_t matches_counters_compact(matches_counters_t* const counters) {
  const uint64_t* const counts = vector_get_mem(counters->counts,uint64_t);
  int64_t i = vector_get_used(counters->counts)-1;
  while (i>=0 && counts[i]==0) --i;
  vector_set_used(counters->counts,++i);
  return i;
}
GEM_INLINE void matches_counters_compute_matches_to_decode(
    matches_counters_t* const counters,const uint64_t min_decoded_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches,
    uint64_t* const strata_to_decode,uint64_t* const matches_to_decode_from_last_stratum) {
  // Compact counters (Shrink the counters to the last non-zero stratum)
  const uint64_t max_strata = matches_counters_compact(counters); // Strata is one based
  if (max_strata==0) return;
  const uint64_t* const counts = vector_get_mem(counters->counts,uint64_t);
  uint64_t current_stratum=0, total_matches=0, total_complete_strata=0;
  // Maximum stratum to decode (increased by @min_decoded_strata & @max_reported_matches)
  while (current_stratum < max_strata && (total_complete_strata < min_decoded_strata || total_matches < max_reported_matches)) {
    total_matches += counts[current_stratum];
    if (total_matches > max_reported_matches) {
      total_matches -= counts[current_stratum];
      break;
    }
    ++current_stratum;
    if (total_matches > 0) ++total_complete_strata;
  }
  // Maximum stratum to decode (increased by @min_reported_matches)
  if (current_stratum < max_strata && total_matches < min_reported_matches) {
    total_matches += counts[current_stratum];
    ++current_stratum;
    // Decode partial stratum
    *matches_to_decode_from_last_stratum = min_reported_matches - total_matches;
  } else {
    // Decode full stratum
    *matches_to_decode_from_last_stratum = UINT64_MAX;
  }
  *strata_to_decode = (total_matches!=0) ? current_stratum : 0;
}
