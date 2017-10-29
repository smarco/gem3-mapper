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

#include "matches/classify/matches_classify.h"

/*
 * Matches Classify
 */
void matches_classify(matches_t* const matches) {
  // Parameters
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  matches_classification_t* const classification = &matches->classification;
  // Unmapped
  if (num_matches == 0) {
    classification->matches_class = matches_class_unmapped;
    classification->delta_group = 0;
    classification->wdelta_group = 0;
    return;
  }
  // Unique
  matches_metrics_t* const metrics = &matches->metrics;
  if (num_matches == 1) {
    classification->matches_class = matches_class_unique;
    classification->delta_group = -1;
    classification->wdelta_group =
        (int64_t)matches->max_complete_stratum - (int64_t)metrics->min1_edit_distance;
    return;
  }
  // Multi-map
  match_trace_t* const primary_match = matches_get_primary_match(matches);
  classification->matches_class = matches_class_mmap;
  if (primary_match->swg_score != metrics->max1_swg_score ||
      primary_match->event_distance != metrics->min1_event_distance ||
      primary_match->edit_distance != metrics->min1_edit_distance) {
    classification->delta_group = 0; // Tie
  } else if (metrics->max1_swg_score == metrics->max2_swg_score ||
             metrics->min1_event_distance == metrics->min2_event_distance ||
             metrics->min1_edit_distance == metrics->min2_edit_distance) {
    classification->delta_group = 0; // Tie
  } else {
    // General Multi-map
    classification->delta_group =
        (int64_t)metrics->min2_edit_distance - (int64_t)metrics->min1_edit_distance;
  }
  classification->wdelta_group =
      (int64_t)MIN(matches->max_complete_stratum,metrics->min2_edit_distance) -
      (int64_t)metrics->min1_edit_distance;
}
/*
 * Matches Condition Tests
 */
bool matches_classify_pattern_viable(pattern_t* const pattern) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // All characters are wildcards
  return (key_length!=num_wildcards && key_length!=0);
}
bool matches_classify_max_matches_searched(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern) {
  // Parameters
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  // Hard max-matches limit
  //  const uint64_t max_searched_matches = search_parameters->select_parameters.max_searched_matches;
  //  // Check total number of matches found so far
  //  if (matches_get_num_match_traces(matches) >= max_searched_matches) {
  //    PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_MAX_MATCHES);
  //    PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_HIT);
  //    return true; // Done!
  //  }
  // Soft max-matches limit
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (select_parameters->min_reported_strata_nominal == 0 &&
      num_matches >= select_parameters->max_searched_matches) {
    match_trace_t** const match_traces = matches_get_match_traces(matches);
    // All matches beyond top_match will have lower swg-score
    const match_trace_t* const top_match = match_traces[select_parameters->max_searched_matches-1];
    const uint64_t bounded_edit_distance =
        align_swg_score_compute_max_edit_bound(
            &search_parameters->swg_penalties,
            top_match->swg_score,pattern->key_length);
    // MCS sets the possibility of finding matches with distance max_error_reached+1
    if (bounded_edit_distance <= matches->max_complete_stratum) {
      PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_MAX_MATCHES);
      PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_HIT);
      return true;
    }
  }
  return false;
}
bool matches_classify_min_depth_searched(
    matches_classification_t* const matches_classification) {
  PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS);
  // Select case
  switch (matches_classification->matches_class) {
    case matches_class_unmapped: // Unmapped (not enough search depth)
      return false;   // 0:0:0+0
    case matches_class_unique:
      if (matches_classification->wdelta_group > 1) {
        return true;  // 0:1:0+0:0
      } else {
        return false; // 0:1+0:0:0 || 0:0+0:1:0
      }
      break; // Not done
    case matches_class_mmap:
      if (matches_classification->delta_group==0) {
        if (matches_classification->wdelta_group >= 0) {
          return true;  // 0:0+2:0 || 0:0:2+0
        } else {
          return false; // 0+0:2:0
        }
      } else {
        if (matches_classification->wdelta_group <= 0) {
          return false;   // 0+0:1:0:1 || 0:0+1:0:1
        } else if (matches_classification->wdelta_group==1) {
          if (matches_classification->delta_group==1) {
            return true;  // 0:0:1+14:0
          } else {
            return false; // 0:0:1+0:1
          }
        } else { // wdelta_group > 1
          return true;    // 0:0:1:0+1
        }
      }
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}
/*
 * Alignment algorithm/paradigm fallback
 */
bool matches_classify_neighbourhood_fallback(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern) {
  // Classify
  matches_classify(matches);
  // Test search depth
  if (matches_classify_min_depth_searched(&matches->classification)) return false; // Done!
  // Test max-matches searched
  if (matches_classify_max_matches_searched(matches,search_parameters,pattern)) return false; // Done!
  // Otherwise, return true
  return true;
}
bool matches_classify_local_alignment_fallback(
    matches_t* const matches,
    const local_alignment_t alignment_local) {
  // Local alignment test
  return (alignment_local!=local_alignment_never && !matches_is_mapped(matches));
}
/*
 * Adjusting search maximum error
 */
uint64_t matches_classify_adjust_max_error_by_strata_after_best(
    matches_t* const matches,
    const uint64_t max_search_error,
    const uint64_t complete_strata_after_best) {
  /*
   * Control delta error adjustment
   *   If delta parameter is set (and is below the maximum number of mismatches),
   *   finds the minimum non zero stratum (mnzs) and adjusts
   *   the maximum number of mismatches to (mnzs+delta)
   */
  if (matches_is_mapped(matches)) {
    matches_metrics_t* const metrics = &matches->metrics;
    if (metrics->min1_edit_distance + complete_strata_after_best < max_search_error) {
      return metrics->min1_edit_distance + complete_strata_after_best;
    }
  }
  return max_search_error;
}
uint64_t matches_classify_compute_max_search_error(
    matches_t* const matches,
    pattern_t* const pattern,
    const uint64_t proper_length) {
  PROF_INC_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS);
  // Classify
  matches_classify(matches);
  // Max feasible error
  const uint64_t max_complete_stratum = matches->max_complete_stratum;
  const uint64_t max_pattern_error = MAX(pattern->key_length/(2*proper_length),1);
  const uint64_t max_feasible_error = (max_complete_stratum > 0) ?  2*max_complete_stratum-1 : max_pattern_error;
  // Compute ideal search error
  const int64_t wdelta_group = matches->classification.wdelta_group;
  uint64_t max_required_error;
  switch (matches->classification.matches_class) {
    case matches_class_unmapped: // Unmapped (not enough search depth)
      max_required_error = max_feasible_error;
      break;
    case matches_class_unique: {
      const uint64_t min_edit_distance = matches->metrics.min1_edit_distance;
      if (wdelta_group > 1) {
        max_required_error = max_complete_stratum - 1;
      } else {
        max_required_error = min_edit_distance + 2;
      }
      break;
    }
    case matches_class_mmap: {
      const uint64_t min_edit_distance = matches->metrics.min1_edit_distance;
      if (matches->classification.delta_group==0) {
        if (matches->classification.wdelta_group>=0) {
          max_required_error = max_complete_stratum - 1;  // 0:0+2:0 || 0:0:2+0
        } else {
          max_required_error = min_edit_distance - 1;     // 0+0:2:0
        }
      } else {
        if (matches->classification.wdelta_group<=0) {
          if (matches->classification.delta_group==1) {
            max_required_error = min_edit_distance;      // 0+0:1:1:0 || 0:0+1:1:0
          } else {
            max_required_error = min_edit_distance + 1;  // 0+0:1:0:1 || 0:0+1:0:1
          }
        } else if (matches->classification.wdelta_group==1) {
          if (matches->classification.delta_group==1) {
            max_required_error = max_complete_stratum - 1;  // 0:0:1+14:0
          } else {
            max_required_error = min_edit_distance + 1;     // 0:0:1+0:1
          }
        } else { // wdelta_group > 1
          max_required_error = max_complete_stratum - 1;    // 0:0:1:0+1
        }
      }
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Account for both limits
  return MIN(max_required_error,max_feasible_error);
}



