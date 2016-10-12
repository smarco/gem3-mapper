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

#include "matches/matches_test.h"
#include "matches/classify/matches_classify.h"

/*
 * Matches Condition Tests
 */
bool matches_test_max_matches_reached(
    matches_t* const matches,
    const uint64_t mcs,
    const uint64_t key_length,
    search_parameters_t* const search_parameters) {
  // Parameters
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
  // Check matches
  if (matches_is_mapped(matches) && select_parameters->min_reported_strata_nominal==0) {
    const uint64_t num_matches = matches_get_num_match_traces(matches);
    if (num_matches >= select_parameters->max_reported_matches) {
      match_trace_t** const match_traces = matches_get_match_traces(matches);
      // All matches beyond top_match will have lower swg-score
      const match_trace_t* const top_match = match_traces[select_parameters->max_reported_matches-1];
      const uint64_t bounded_edit_distance =
          align_swg_score_compute_max_edit_bound(
              &search_parameters->swg_penalties,top_match->swg_score,key_length);
      // MCS sets the possibility of finding matches with distance max_error_reached+1
      if (bounded_edit_distance <= mcs) return true;
    }
  }
  return false;
}
bool matches_test_accuracy_reached(
    matches_t* const matches,
    const uint64_t mcs,
    const uint64_t key_length,
    search_parameters_t* const search_parameters,
    const uint64_t max_complete_error,
    uint64_t* const max_complete_error_required) {
  PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS);
  // Parameters
  const uint64_t delta = search_parameters->complete_strata_after_best_nominal; // (default = 1)
  // Classify
  const uint64_t min_edit_distance = matches_metrics_get_min_edit_distance(&matches->metrics);
  matches_classify(matches);
  switch (matches->matches_class) {
    case matches_class_unmapped: // Unmapped (not enough search depth)
      *max_complete_error_required = mcs + delta; // Adjust max-error
      PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_UNMAPPED);
      break; // Not done
    case matches_class_tie_perfect:
    case matches_class_tie:
      if (min_edit_distance <= mcs+1) { // (0:2+0:0) && (0+0:2:0) but not (0+0:0:2)
        if (!search_parameters->search_paired_parameters.paired_end_search) {
          PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_HIT);
          return true; // Done!
        }
      }
      *max_complete_error_required = mcs + delta; // Adjust max-error
      PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_TIE);
      break; // Not done
    case matches_class_mmap_d1:
      if (min_edit_distance+1 <= mcs) { // (0:1+1:0)
        if (!search_parameters->search_paired_parameters.paired_end_search) {
          PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_HIT);
          return true; // Done!
        }
      }
      PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_MMAP_D1);
      // no break
    case matches_class_mmap:
    case matches_class_unique:
      // Frontier-case 0:1+0 & Beyond-case 0:0+0:0:1
      if (min_edit_distance+1 < mcs) {
        PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_HIT);
        return true; // Done!
      }
      // Adjust max-error
      *max_complete_error_required = MIN(max_complete_error,min_edit_distance+1);
      if (*max_complete_error_required+1 <= mcs) {
        PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_HIT);
        return true; // Done!
      }
      // PROFILE
#ifdef GEM_PROFILE
      if (matches->matches_class==matches_class_mmap) PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_MMAP);
      if (matches->matches_class==matches_class_unique) PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_UNIQUE);
#endif
      break; // Not done
    default:
      GEM_INVALID_CASE();
      break;
  }
  // PROFILE
#ifdef GEM_PROFILE
  if (min_edit_distance+1 == mcs) PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_MAP_FRONTIER);
  if (min_edit_distance >= mcs) PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_MAP_INCOMPLETE);
#endif
  // Candidates
  const bool max_matches_reached = matches_test_max_matches_reached(matches,mcs,key_length,search_parameters);
  if (max_matches_reached) {
    PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_MAX_MATCHES);
    PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_HIT);
    return true; // Done!
  }
  // Return
  PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_MISS);
  return false;
}
