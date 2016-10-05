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
 * DESCRIPTION:
 *   Approximate-String-Matching (ASM) search control functions (regulate
 *   the depth of the search)
 */

#include "approximate_search/approximate_search_control.h"
#include "matches/classify/matches_classify.h"
#include "matches/classify/matches_classify_logit.h"
#include "matches/classify/matches_classify_logit_models.h"

/*
 * Search Limits
 */
void asearch_control_adjust_current_max_error(
    approximate_search_t* const search,
    matches_t* const matches) {
  const uint64_t current_max_complete_error = search->current_max_complete_error;
  const uint64_t delta = search->search_parameters->complete_strata_after_best_nominal;
  /*
   * Control delta error adjustment
   *   If delta parameter is set (and is below the maximum number of mismatches),
   *   finds the minimum non zero stratum (mnzs) and adjusts
   *   the maximum number of mismatches to (mnzs+delta)
   */
  if (delta < current_max_complete_error) {
    const uint64_t fms = matches_metrics_get_min_edit_distance(&matches->metrics);
    if (fms+delta < current_max_complete_error) {
      search->current_max_complete_error = fms+delta;
    }
  }
}
/*
 * Pattern test
 */
bool asearch_control_test_pattern(
    approximate_search_t* const search) {
  // Parameters
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // All characters are wildcards
  return (key_length!=num_wildcards && key_length!=0);
}
/*
 * Accuracy test
 */
bool asearch_control_test_accuracy__adjust_depth(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  const uint64_t mcs = search->region_profile.num_filtered_regions;
  const uint64_t delta = 1;
  // Test pattern
  if (search->pattern.num_wildcards > search->current_max_complete_error) return true; // Done!
  // Unmapped (not enough search depth)
  if (!matches_is_mapped(matches)) {
    search->current_max_complete_error = mcs + delta;  // Adjust max-error
    if (search->pattern.num_wildcards > search->current_max_complete_error) return true; // Done!
    return false;
  }
  // Ties (unsolvable)
  matches_classify(matches);
  if (matches->matches_class==matches_class_tie_perfect ||
      mcs >= search_parameters->complete_search_error_nominal+1) {
    return true; // Done!
  }
  // Frontier-case 0:1+0 & Beyond-case 0:0+0:0:1
  const uint64_t min_edit_distance = matches_metrics_get_min_edit_distance(&matches->metrics);
  if (min_edit_distance+1 >= mcs) {
    search->current_max_complete_error = MIN(search->current_max_complete_error,min_edit_distance+1);
    if (search->pattern.num_wildcards > search->current_max_complete_error) return true; // Done!
    return false;
  } else {
    return true; // Done!
  }
}
/*
 * Local-alignment
 */
bool asearch_control_test_local_alignment(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  // Local alignment test
  return (search_parameters->local_alignment!=local_alignment_never && !matches_is_mapped(matches));
}


///*
// * Search Control Stages
// */
//bool asearch_control_filter_ahead_candidates(
//    approximate_search_t* const search,
//    matches_t* const matches);
//
///*
// * Search Fulfilled & Predictors
// */
//void asearch_control_compute_predictors(
//    approximate_search_t* const search,
//    matches_t* const matches,
//    matches_predictors_t* const predictors);
//bool asearch_control_fulfilled(
//    approximate_search_t* const search,
//    matches_t* const matches);


///*
// * Search Control Stages
// */
//bool asearch_control_filter_ahead_candidates(
//    approximate_search_t* const search,
//    matches_t* const matches) {
//  // Parameters
//  const search_parameters_t* const search_parameters = search->search_parameters;
//  // Determines when the search is done following the mapping criteria
//  switch (search_parameters->mapping_mode) {
//    case mapping_adaptive_filtering_fast:
//      return true;
//    case mapping_adaptive_filtering_complete:
//      return search_parameters->complete_strata_after_best_nominal < search->current_max_complete_error;
//    default:
//      GEM_INVALID_CASE();
//      break;
//  }
//  return false;
//}
///*
// * Search Fulfilled & Predictors
// */
//void asearch_control_compute_predictors(
//    approximate_search_t* const search,
//    matches_t* const matches,
//    matches_predictors_t* const predictors) {
//  // TODO Update MCS
//  matches_predictors_compute(matches,predictors,&search->metrics,search->current_mcs);
//}
//bool asearch_control_fulfilled(
//    approximate_search_t* const search,
//    matches_t* const matches) {
//  // Parameters
//  const search_parameters_t* const search_parameters = search->search_parameters;
//  if (matches==NULL) return false;
//  // Determines when the search is done following the mapping criteria
//  switch (search_parameters->mapping_mode) {
//    case mapping_adaptive_filtering_fast: {
//      // TODO     asearch_control_adjust_max_differences_using_strata(search,matches);
//      if (matches->max_complete_stratum <= 1) return false;
//      switch (matches_class) {
//        case matches_class_tie_d0:
//        case matches_class_tie_d1:
//          return (matches->metrics.min1_edit_distance <= 1);
//        case matches_class_unique: {
//          matches_predictors_t predictors;
//          asearch_control_compute_predictors(search,matches,&predictors);
//          const double probability =
//              matches_classify_logit_unique(&predictors,&logit_model_single_end_default);
//          return (probability >= MATCHES_UNIQUE_CI);
//        }
//        case matches_class_unmapped:
//        case matches_class_mmap:
//        default:
//          return false;
//      }
//      return false;
//    }
//    case mapping_adaptive_filtering_complete:
//      // TODO
//    default:
//      GEM_INVALID_CASE();
//      break;
//  }
//  return false;
//}

