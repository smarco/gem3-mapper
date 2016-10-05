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
 *   Approximate-String-Matching (ASM) using hybrid techniques.
 *   Combines adaptive-filtering (AF) with neighborhood-search (NS)
 *   to perform efficient complete searches
 */

#include "approximate_search/approximate_search_hybrid.h"
#include "approximate_search/approximate_search_control.h"
#include "approximate_search/approximate_search_stages.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_generate_candidates.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "filtering/region_profile/region_profile.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "matches/classify/matches_classify.h"

/*
 * Control
 */
asearch_stage_t as_hybrid_control_begin(approximate_search_t* const search) {
  // Test pattern
  if (asearch_control_test_pattern(search)) {
    PROF_INC_COUNTER(GP_AS_FILTERING_ADATIVE_CALL);
    return asearch_stage_filtering_adaptive;
  } else {
    return asearch_stage_end;
  }
}
asearch_stage_t as_hybrid_control_filtering_adaptive_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROF_ADD_COUNTER(GP_AS_FILTERING_ADATIVE_MCS,search->region_profile.num_filtered_regions);
  // Select state
  switch (search->processing_state) {
    case asearch_processing_state_no_regions:
      search->current_max_complete_error = 1;
      if (search->pattern.num_wildcards > search->current_max_complete_error) {
        return asearch_stage_end;
      } else {
        return asearch_stage_neighborhood;
      }
    case asearch_processing_state_candidates_verified:
      if (!asearch_control_test_accuracy__adjust_depth(search,matches)) {
#ifdef GEM_PROFILE
        PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL);
        if (!matches_is_mapped(matches)) {
          PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL_UNMAPPED);
        } else {
          const uint64_t mcs = search->region_profile.num_filtered_regions;
          const uint64_t min_edit_distance = matches_metrics_get_min_edit_distance(&matches->metrics);
          if (min_edit_distance+1 == mcs) PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL_MAP_FRONTIER);
          if (min_edit_distance >= mcs) PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL_MAP_INCOMPLETE);
        }
#endif
        return asearch_stage_neighborhood;
      } else {
        return asearch_stage_end;
      }
    default:
      GEM_INVALID_CASE();
      break;
  }
  return asearch_stage_end;
}
asearch_stage_t as_hybrid_control_neighborhood_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROF_ADD_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_MCS,search->current_max_complete_stratum);
  if (asearch_control_test_local_alignment(search,matches)) {
    PROF_INC_COUNTER(GP_AS_LOCAL_ALIGN_CALL);
    return asearch_stage_local_alignment;
  } else {
    return asearch_stage_end;
  }
}
/*
 * Hybrid mapping [GEM-workflow 5.0]
 */
void approximate_search_hybrid(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Process proper search-stage
  while (true) {
    switch (search->search_stage) {
      case asearch_stage_begin: // Search Start. Check basic cases
        search->search_stage = as_hybrid_control_begin(search);
        break;
      case asearch_stage_filtering_adaptive: // Exact-Filtering (Adaptive)
        approximate_search_exact_filtering_adaptive(search,matches);
        search->search_stage = as_hybrid_control_filtering_adaptive_next_state(search,matches);
        break;
      case asearch_stage_filtering_adaptive_finished:
        search->search_stage = as_hybrid_control_filtering_adaptive_next_state(search,matches);
        break;
      case asearch_stage_neighborhood:
        approximate_search_neighborhood(search,matches);
        search->search_stage = as_hybrid_control_neighborhood_next_state(search,matches);
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
/*
 * Approximate Complete-Search based on PURE filtering+NS-search
 */
void approximate_search_hybrid_complete_search(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Exact Search
  if (search->current_max_complete_error==0) {
    approximate_search_neighborhood_exact_search(search,matches);
    return;
  }
  // Compute the region profile
  approximate_search_region_profile_adaptive(search,region_profile_adaptive);
  // Generate exact-candidates
  region_profile_schedule_filtering_exact(region_profile);
  approximate_search_generate_candidates_exact(search);
  //  // Verify candidates
  //  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  //  filtering_candidates_process_candidates(filtering_candidates,pattern,true);
  //  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  //  // Align candidates
  //  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  // Check max-complete-error to be reached
  approximate_search_update_mcs(search,region_profile->num_filtered_regions + pattern->num_wildcards);
  if (search->current_max_complete_stratum < search->current_max_complete_error + 1) {
    search_parameters_t* const search_parameters = search->search_parameters;
    // Prepare region-profile (fill gaps)
    region_profile_fill_gaps(region_profile,pattern->key,pattern->key_length,
        search_parameters->allowed_enc,pattern->num_wildcards);
    region_profile_merge_small_regions(region_profile,search->archive->fm_index->proper_length);
    // Complete Search with NS-Search (preconditioned)
    approximate_search_neighborhood_search_partition_preconditioned(search,matches);
  }
  // Finish Search
  approximate_search_end(search,matches);
}
