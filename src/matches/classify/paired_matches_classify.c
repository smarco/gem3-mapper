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

#include "matches/classify/paired_matches_classify.h"
#include "matches/classify/matches_classify.h"

/*
 * Paired-matches Classify
 */
void paired_matches_classify(paired_matches_t* const paired_matches) {
  // Parameters
  const uint64_t num_pairs = paired_matches_get_num_maps(paired_matches);
  matches_classification_t* const classification = &paired_matches->classification;
  // Unmapped
  if (num_pairs == 0) {
    classification->matches_class = matches_class_unmapped;
    classification->delta_group = 0;
    classification->wdelta_group = 0;
    return;
  }
  // Unique
  const uint64_t max_complete_stratum = paired_matches_get_max_complete_stratum(paired_matches);
  matches_metrics_t* const metrics = &paired_matches->metrics;
  if (num_pairs == 1) {
    classification->matches_class = matches_class_unique;
    classification->delta_group = -1;
    classification->wdelta_group = (int64_t)max_complete_stratum - metrics->min1_edit_distance;
    return;
  }
  // Multi-map
  classification->matches_class = matches_class_mmap;
  if (metrics->max1_swg_score == metrics->max2_swg_score ||
      metrics->min1_event_distance == metrics->min2_event_distance ||
      metrics->min1_edit_distance == metrics->min2_edit_distance) {
    classification->delta_group = 0; // Tie
  } else {
    // General Multi-map
    classification->delta_group =
        (int64_t)metrics->min2_edit_distance - (int64_t)metrics->min1_edit_distance;
  }
  classification->wdelta_group =
      (int64_t)MIN(max_complete_stratum,metrics->min2_edit_distance) -
      (int64_t)metrics->min1_edit_distance;
}
/*
 * Paired Matches Accuracy Tests
 */
bool paired_matches_classify_search_accomplished(
    paired_matches_t* const paired_matches) {
  // Classify
  paired_matches_classify(paired_matches);
  // Test search depth
  if (matches_classify_min_depth_searched(&paired_matches->classification)) return true; // Done!
//  // Test max-matches searched
//  if (matches_classify_max_matches_searched(matches,search_parameters,pattern)) return false; // Done!
  // Otherwise, return false
  return false;
}
/*
 * Subdominant End
 */
bool paired_matches_classify_subdominant_end(
    paired_matches_t* const paired_matches,
    matches_t* const candidate_matches,
    match_trace_t* const extended_match) {
  // Unmapped
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches == 0) return false;
  // The end needs to have less distance than the subdominant pair
  const uint64_t candidate_end_min_edit_distance = candidate_matches->max_complete_stratum;
  const uint64_t extended_end_min_edit_distance = extended_match->edit_distance;
  const uint64_t min_expected_edit_distance = candidate_end_min_edit_distance + extended_end_min_edit_distance;
  matches_metrics_t* const metrics = &paired_matches->metrics;
  const uint64_t subdominant_pair_edit_distance = metrics->min2_edit_distance;
  return subdominant_pair_edit_distance <= min_expected_edit_distance;
}
