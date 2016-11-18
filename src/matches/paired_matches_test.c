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

#include "matches/paired_matches_test.h"
#include "matches/classify/matches_classify.h"

/*
 * Paired Matches Condition Tests
 */
bool paired_matches_test_accuracy_reached(
    paired_matches_t* const paired_matches,
    search_parameters_t* const search_parameters) {
  // Parameters
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
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
    	const uint64_t max_reported_matches = search_parameters->select_parameters_report.max_reported_matches;
      const uint64_t max_complete_stratum =
          ((matches_end1->max_complete_stratum!=ALL) ? matches_end1->max_complete_stratum : 0) +
          ((matches_end2->max_complete_stratum!=ALL) ? matches_end2->max_complete_stratum : 0);
      const uint64_t num_paired_maps = paired_matches_get_num_maps(paired_matches);
      paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
      return (paired_map->edit_distance < max_complete_stratum && num_paired_maps >= max_reported_matches);
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}
