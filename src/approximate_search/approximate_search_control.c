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
#include "matches/matches_test.h"

/*
 * Search Limits
 */
void asearch_control_adjust_current_max_error(
    approximate_search_t* const search,
    matches_t* const matches) {
  /*
   * Control delta error adjustment
   *   If delta parameter is set (and is below the maximum number of mismatches),
   *   finds the minimum non zero stratum (mnzs) and adjusts
   *   the maximum number of mismatches to (mnzs+delta)
   */
  if (matches_is_mapped(matches)) {
    const uint64_t current_max_complete_error = search->current_max_complete_error;
    const uint64_t delta = search->search_parameters->complete_strata_after_best_nominal;
    const uint64_t min_edit_distance = matches_metrics_get_min_edit_distance(&matches->metrics);
    if (min_edit_distance+delta < current_max_complete_error) {
      search->current_max_complete_error = min_edit_distance+delta;
    }
  }
}
bool asearch_control_max_matches_reached(
    approximate_search_t* const search,
    matches_t* const matches) {
  search_parameters_t* const search_parameters = search->search_parameters;
  return matches_test_max_matches_reached(
      matches,search->region_profile.num_filtered_regions,
      search->pattern.key_length,search_parameters);
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
  // Test pattern
  if (search->pattern.num_wildcards > search->current_max_complete_error) return true; // Done!
  // Test accuracy
  search_parameters_t* const search_parameters = search->search_parameters;
  const uint64_t mcs = search->region_profile.num_filtered_regions; // (max_error_reached = mcs-1)
  const bool accuracy_reached = matches_test_accuracy_reached(
      matches,mcs,search->pattern.key_length,search_parameters,
      search->current_max_complete_error,&search->current_max_complete_error);
  if (accuracy_reached) return true; // Done!
  // Check error-scheduled
  if (search->pattern.num_wildcards > search->current_max_complete_error) return true; // Done!
  // Otherwise, return false
  return false;
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
  return (search_parameters->alignment_local!=local_alignment_never && !matches_is_mapped(matches));
}

