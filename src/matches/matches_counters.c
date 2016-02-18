/*
 * PROJECT: GEMMapper
 * FILE: matches_counters.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/matches_counters.h"

/*
 * Constants
 */
#define MATCHES_INIT_COUNTERS          200

/*
 * Setup
 */
matches_counters_t* matches_counters_new() {
  // Allocate handler
  matches_counters_t* const counters = mm_alloc(matches_counters_t);
  // Init Counters
  counters->counts = vector_new(MATCHES_INIT_COUNTERS,uint64_t);
  counters->total_count = 0;
  // Return
  return counters;
}
void matches_counters_clear(matches_counters_t* const counters) {
  vector_clear(counters->counts);
  counters->total_count = 0;
}
void matches_counters_delete(matches_counters_t* const counters) {
  vector_delete(counters->counts);
  mm_free(counters);
}

/*
 * Counters
 */
uint64_t matches_counters_get_num_counters(matches_counters_t* const counters) {
  return vector_get_used(counters->counts);
}
uint64_t* matches_counters_get_counts(matches_counters_t* const counters) {
  return vector_get_mem(counters->counts,uint64_t);
}
uint64_t matches_counters_get_count(matches_counters_t* const counters,const uint64_t distance) {
  return *vector_get_elm(counters->counts,distance,uint64_t);
}
uint64_t matches_counters_get_total_count(matches_counters_t* const counters) {
  return counters->total_count;
}
void matches_counters_add(
    matches_counters_t* const counters,
    const uint64_t distance,
    const uint64_t num_matches) {
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
void matches_counters_sub(
    matches_counters_t* const counters,
    const uint64_t distance,
    const uint64_t num_matches) {
  uint64_t* const counts = vector_get_elm(counters->counts,distance,uint64_t);
  *counts -= num_matches;
  counters->total_count -= num_matches;
}
/*
 * Utils
 */
uint64_t matches_counters_compact(matches_counters_t* const counters) {
  const uint64_t* const counts = vector_get_mem(counters->counts,uint64_t);
  int64_t i = vector_get_used(counters->counts)-1;
  while (i>=0 && counts[i]==0) --i;
  vector_set_used(counters->counts,++i);
  return i;
}
void matches_counters_compute_matches_to_decode(
    matches_counters_t* const counters,
    const uint64_t min_reported_strata,
    const uint64_t min_reported_matches,
    const uint64_t max_reported_matches,
    uint64_t* const reported_strata,
    uint64_t* const last_stratum_reported_matches) {
  // Compact counters (Shrink the counters to the last non-zero stratum)
  const uint64_t max_strata = matches_counters_compact(counters); // Strata is one based
  if (max_strata==0) return;
  const uint64_t* const counts = vector_get_mem(counters->counts,uint64_t);
  uint64_t current_stratum=0, total_matches=0, total_complete_strata=0;
  // Maximum stratum to decode (increased by @min_reported_strata)
  while (current_stratum < max_strata && total_complete_strata < min_reported_strata) {
    total_matches += counts[current_stratum];
    if (total_matches > 0) ++total_complete_strata;
    ++current_stratum;
  }
  // Maximum stratum to decode (increased by @max_reported_matches)
  while (current_stratum < max_strata && total_matches < max_reported_matches) {
    total_matches += counts[current_stratum];
    if (total_matches > max_reported_matches) {
      total_matches -= counts[current_stratum];
      break;
    }
    ++current_stratum;
  }
  // Maximum stratum to decode (increased by @min_reported_matches)
  if (current_stratum < max_strata && total_matches < min_reported_matches) {
    const uint64_t remaining_matches = min_reported_matches - total_matches;
    const uint64_t report_last_stratum = MIN(remaining_matches,counts[current_stratum]);
    *last_stratum_reported_matches = report_last_stratum; // Decode partial stratum
    ++current_stratum;
    total_matches += report_last_stratum;
  } else {
    *last_stratum_reported_matches = UINT64_MAX; // Decode full stratum
  }
  *reported_strata = (total_matches!=0) ? current_stratum : 0;
}
/*
 * Display
 */
void matches_counters_print_account_mcs(
    FILE* const stream,
    const uint64_t current_counter_pos,
    const uint64_t num_zeros) {
  uint64_t i = 0;
  if (current_counter_pos==0) {
    fprintf(stream,"0");
    i=1;
  }
  for (;i<num_zeros;++i) {
    fprintf(stream,":0");
  }
  fprintf(stream,"+0");
}
void matches_counters_print(
    FILE* const stream,
    matches_counters_t* const matches_counter,
    const uint64_t mcs) {
  const uint64_t num_counters = matches_counters_get_num_counters(matches_counter);
  // Zero counters
  if (gem_expect_false(num_counters==0)) {
    if (mcs==0 || mcs==ALL) {
      fprintf(stream,"0");
    } else {
      matches_counters_print_account_mcs(stream,0,mcs);
    }
    return;
  }
  // Print counters
  uint64_t i = 0;
  const uint64_t* counters = matches_counters_get_counts(matches_counter);
  while (i < num_counters) {
    // Print Counter
    if (i>0) fprintf(stream,(mcs==i?"+":":"));
    fprintf(stream,"%lu",*counters);
    i++; // Next (+1)
    ++counters;
  }
  // Account for MCS
  if (i<=mcs) matches_counters_print_account_mcs(stream,i,mcs-i);
}
