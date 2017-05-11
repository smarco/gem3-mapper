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

#include "matches/paired_matches_accuracy.h"
#include "matches/classify/matches_classify.h"

/*
 * Paired Matches Accuracy Tests
 */
bool paired_matches_accuracy_reached(
    paired_matches_t* const paired_matches,
    search_parameters_t* const search_parameters) {
  // Parameters
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  //  // Check total matches found
  //  const uint64_t max_searched_paired_matches = select_parameters->max_searched_paired_matches;
  //  if (paired_matches_get_num_maps(paired_matches) >= max_searched_paired_matches) break;
  // Classify
  paired_matches_classify(paired_matches);
  switch (paired_matches->paired_matches_class) {
    case paired_matches_class_unmapped:
    case paired_matches_class_mmap_d1:
    case paired_matches_class_mmap:
    case paired_matches_class_unique:
      return false;
    case paired_matches_class_tie_perfect:
    case paired_matches_class_tie: {
    	const uint64_t max_reported_matches = search_parameters->select_parameters.max_reported_matches;
      const uint64_t max_complete_stratum =
          ((matches_end1->max_complete_stratum!=ALL) ? matches_end1->max_complete_stratum : 0) +
          ((matches_end2->max_complete_stratum!=ALL) ? matches_end2->max_complete_stratum : 0);
      const uint64_t num_maps = paired_matches_get_num_maps(paired_matches);
      paired_map_t* const paired_map = paired_matches_get_primary_map(paired_matches);
      return (paired_map->edit_distance < max_complete_stratum && num_maps >= max_reported_matches);
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}
/*
 * Subdominant End
 */
bool paired_matches_accuracy_subdominant_end(
    paired_matches_t* const paired_matches,
    search_parameters_t* const search_parameters,
    matches_t* const candidate_matches,
    match_trace_t* const extended_match) {
  // Parameters
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  // Basic cases
  const uint64_t max_reported_matches = select_parameters->max_reported_matches;
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches == 0 || num_matches < max_reported_matches) return false;
  // The end needs to have less distance than the last pair
  paired_map_t* const paired_map = paired_matches_get_maps(paired_matches)[max_reported_matches-1];
  if (extended_match->edit_distance >= paired_map->edit_distance) {
    return true; // Distance subdominant
  }
  // Check difference
  const uint64_t max_complete_stratum =
      ((candidate_matches->max_complete_stratum!=ALL) ? candidate_matches->max_complete_stratum : 0);
  const uint64_t match_trace_extended_max_distance = paired_map->edit_distance - extended_match->edit_distance;
  if (match_trace_extended_max_distance < max_complete_stratum) {
    return true; // Already found everything interesting
  }
  return false;
}
