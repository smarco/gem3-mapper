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
 *   Approximate-String-Matching (ASM) using adaptive-filtering techniques (AF)
 */

#include "approximate_search/approximate_search_stages.h"
#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "matches/classify/matches_classify.h"

/*
 * Control
 */
asearch_stage_t as_filtering_control_begin(approximate_search_t* const search) {
  // Test pattern
  if (matches_classify_pattern_viable(&search->pattern)) {
    PROF_INC_COUNTER(GP_AS_FILTERING_ADATIVE_CALL);
    return asearch_stage_filtering_adaptive;
  } else {
    return asearch_stage_end;
  }
}
asearch_stage_t as_filtering_control_filtering_adaptive_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROF_ADD_COUNTER(GP_AS_FILTERING_ADATIVE_MCS,search->region_profile.num_filtered_regions);
  // Select state
  switch (search->processing_state) {
    case asearch_processing_state_no_regions:
      return asearch_stage_end;
    case asearch_processing_state_candidates_verified:
      if (matches_classify_local_alignment_fallback(
          matches,search->search_parameters->alignment_local)) {
        PROF_INC_COUNTER(GP_AS_LOCAL_ALIGN_CALL);
        return asearch_stage_local_alignment;
      } else {
        return asearch_stage_end;
      }
    default:
      GEM_INVALID_CASE();
      break;
  }
  return asearch_stage_end;
}
/*
 * Adaptive mapping [GEM-workflow 4.0]
 */
void approximate_search_filtering_adaptive(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Process proper search-stage
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search begin
        search->search_stage = as_filtering_control_begin(search);
        break;
      case asearch_stage_filtering_adaptive: // Exact-Filtering (Adaptive)
        // approximate_search_exact_filtering_adaptive_cutoff(search,matches);
        approximate_search_exact_filtering_adaptive(search,matches);
        search->search_stage = as_filtering_control_filtering_adaptive_next_state(search,matches);
        break;
      case asearch_stage_filtering_adaptive_finished:
        matches_update_mcs(matches,search->region_profile.num_filtered_regions); // (+ pattern->num_wildcards)
        matches_update_limited_exact_matches(matches,search->num_limited_exact_matches);
        search->search_stage = as_filtering_control_filtering_adaptive_next_state(search,matches);
        break;
      case asearch_stage_local_alignment: // Local alignments
        approximate_search_align_local(search,matches);
        search->search_stage = asearch_stage_end; // Next State
        break;
      case asearch_stage_end:
        approximate_search_end(search,matches);
        return;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
