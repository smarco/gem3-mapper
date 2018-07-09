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

#include "matches/matches_counters.h"

/*
 * Constants
 */
#define MATCHES_INIT_COUNTERS          200

/*
 * Setup
 */
matches_counters_t* matches_counters_new(void) {
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
uint64_t matches_counters_get_num_counters(
    matches_counters_t* const matches_counters) {
  return vector_get_used(matches_counters->counts);
}
uint64_t* matches_counters_get_counts(
    matches_counters_t* const matches_counters) {
  return vector_get_mem(matches_counters->counts,uint64_t);
}
uint64_t matches_counters_get_count(
    matches_counters_t* const matches_counters,
    const uint64_t distance) {
  return *vector_get_elm(matches_counters->counts,distance,uint64_t);
}
uint64_t matches_counters_get_total_count(
    matches_counters_t* const matches_counters) {
  return matches_counters->total_count;
}
void matches_counters_add(
    matches_counters_t* const matches_counters,
    const uint64_t distance,
    const uint64_t num_matches) {
  vector_t* const counters = matches_counters->counts;
  // Reserve Memory
  if (distance >= vector_get_used(counters)) {
    vector_reserve(counters,distance+1,true);
    vector_set_used(counters,distance+1);
  }
  // Add matches
  *vector_get_elm(counters,distance,uint64_t) += num_matches;
  matches_counters->total_count += num_matches;
}
void matches_counters_sub(
    matches_counters_t* const matches_counters,
    const uint64_t distance,
    const uint64_t num_matches) {
  uint64_t* const counters = vector_get_elm(matches_counters->counts,distance,uint64_t);
  *counters -= num_matches;
  matches_counters->total_count -= num_matches;
}
/*
 * Utils
 */
uint64_t matches_counters_compact(matches_counters_t* const matches_counters) {
  const uint64_t* const counts = vector_get_mem(matches_counters->counts,uint64_t);
  int64_t i = vector_get_used(matches_counters->counts)-1;
  while (i>=0 && counts[i]==0) --i;
  vector_set_used(matches_counters->counts,++i);
  return i;
}
void matches_counters_compute_matches_to_report(
    matches_counters_t* const matches_counters,
    const uint64_t min_reported_strata,
    const uint64_t max_reported_matches,
    uint64_t* const matches_to_report,
    uint64_t* const strata_to_report) {
  // Compact counters (Shrink the counters to the last non-zero stratum)
  const uint64_t max_strata = matches_counters_compact(matches_counters); // Strata is one based
  if (max_strata==0) {
    *matches_to_report = 0;
    *strata_to_report = 0;
    return;
  }
  const uint64_t* const counts = vector_get_mem(matches_counters->counts,uint64_t);
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
  *matches_to_report = total_matches;
  *strata_to_report = (total_matches!=0) ? current_stratum : 0;
}
uint64_t matches_counters_count_first_subdominant(
    matches_counters_t* const matches_counters) {
  const uint64_t* const counters = vector_get_mem(matches_counters->counts,uint64_t);
  const uint64_t total_counters = vector_get_used(matches_counters->counts);
  uint64_t i;
  for (i=0;i<total_counters;++i) {
    if (counters[i]!=0) break;
  }
  if (i==total_counters) return 0;
  for (++i;i<total_counters;++i) {
    if (counters[i]!=0) return counters[i];
  }
  return 0;
}
void matches_counters_count_delta_edit(
    matches_counters_t* const matches_counters,
    int64_t* const best_edit_distance,
    int64_t* const subdominant_edit_distance) {
  const uint64_t* const counters = vector_get_mem(matches_counters->counts,uint64_t);
  const uint64_t total_counters = vector_get_used(matches_counters->counts);
  uint64_t i;
  *best_edit_distance = -1;
  *subdominant_edit_distance = -1;
  for (i=0;i<total_counters;++i) {
    if (counters[i]!=0) {
      *best_edit_distance = i;
      if (counters[i]>1) {
        *subdominant_edit_distance = i;
        return;
      }
      break;
    }
  }
  if (i==total_counters) return;
  for (++i;i<total_counters;++i) {
    if (counters[i]!=0) {
      *subdominant_edit_distance = i;
      return;
    }
  }
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
    fprintf(stream,"%"PRIu64,*counters);
    i++; // Next (+1)
    ++counters;
  }
  // Account for MCS
  if (i<=mcs) matches_counters_print_account_mcs(stream,i,mcs-i);
}
