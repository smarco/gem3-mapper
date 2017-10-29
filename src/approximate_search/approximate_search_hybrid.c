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
  if (matches_classify_pattern_viable(&search->pattern)) {
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
  search_parameters_t* const search_parameters = search->search_parameters;
  pattern_t* const pattern = &search->pattern;
  if (search_parameters->mapping_mode==mapping_hybrid_sensitive) {
    // Hybrid Sensitive
    if (!matches_classify_neighbourhood_fallback(matches,search_parameters,pattern)) {
      return asearch_stage_end; // Done!
    }
    // Compute search max-error
    const uint64_t proper_length = fm_index_get_proper_length(search->archive->fm_index);
    const uint64_t max_search_error =
        matches_classify_compute_max_search_error(matches,pattern,proper_length);
    if (max_search_error+1 <= matches->max_complete_stratum) return asearch_stage_end; // Done!
    if (pattern->num_wildcards > max_search_error) return asearch_stage_end; // Done!
    // NS
    PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL);
    search->max_search_error = max_search_error;
    return asearch_stage_neighborhood;
  } else {
    // Hybrid Complete
    const uint64_t complete_strata_after_best = search_parameters->complete_strata_after_best_nominal;
    search->max_search_error =
        matches_classify_adjust_max_error_by_strata_after_best(
            matches,search->max_search_error,complete_strata_after_best);
    // Check search depth
    const uint64_t max_complete_stratum = matches->max_complete_stratum;
    if (search->max_search_error+1 <= max_complete_stratum) return asearch_stage_end;
    if (pattern->num_wildcards > search->max_search_error) return asearch_stage_end;
    // NS
    search_parameters->nsearch_parameters.matches_accuracy_cutoff = false;
    PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL);
    return asearch_stage_neighborhood;
  }
}
asearch_stage_t as_hybrid_control_neighborhood_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROF_ADD_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_MCS,matches->max_complete_stratum);
  if (search->search_parameters->mapping_mode==mapping_hybrid_sensitive) {
    // Hybrid Sensitive
    if (matches_classify_local_alignment_fallback(
        matches,search->search_parameters->alignment_local)) {
      PROF_INC_COUNTER(GP_AS_LOCAL_ALIGN_CALL);
      return asearch_stage_local_alignment;
    } else {
      return asearch_stage_end;
    }
  } else {
    // Hybrid Complete
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
        matches_update_mcs(matches,search->region_profile.num_filtered_regions); // (+ pattern->num_wildcards)
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
