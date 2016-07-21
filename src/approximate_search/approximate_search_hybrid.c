/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_hybrid.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_hybrid.h"
#include "approximate_search/approximate_search_stages.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_generate_candidates.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "filtering/region_profile.h"
#include "filtering/region_profile_schedule.h"

/*
 * Control
 */
void as_hybrid_control_filtering_adaptive_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Select state
  switch (search->processing_state) {
    case asearch_processing_state_no_regions:
      search->search_stage = asearch_stage_neighborhood;
      break;
    case asearch_processing_state_candidates_verified:
      // Neighborhood Search
      if (!matches_is_mapped(matches)) {
        search->search_stage = asearch_stage_neighborhood;
      } else {
        search->search_stage = asearch_stage_end;
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void as_hybrid_control_neighborhood_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  search_parameters_t* const search_parameters = search->search_parameters;
  if (search_parameters->local_alignment==local_alignment_never || matches_is_mapped(matches)) {
    search->search_stage = asearch_stage_end;
  } else {  // local_alignment_if_unmapped
    search->search_stage = asearch_stage_local_alignment;
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
        approximate_search_begin(search);
        break;
      case asearch_stage_filtering_adaptive: // Exact-Filtering (Adaptive)
        approximate_search_exact_filtering_adaptive(search,matches);
        as_hybrid_control_filtering_adaptive_next_state(search,matches); // Next State
        break;
      case asearch_stage_neighborhood:
        approximate_search_neighborhood(search,matches);
        as_hybrid_control_neighborhood_next_state(search,matches); // Next State
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
  approximate_search_region_profile_adaptive(search,region_profile_adaptive,search->mm_stack);
  // Generate exact-candidates
  region_profile_schedule_filtering_exact(region_profile);
  approximate_search_generate_candidates_exact(search,matches);
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
        search_parameters->allowed_enc,pattern->num_wildcards,search->mm_stack);
    region_profile_merge_small_regions(region_profile,search->archive->fm_index->proper_length);
    // Complete Search with NS-Search (preconditioned)
    approximate_search_neighborhood_search_partition_preconditioned(search,matches);
  }
  // Finish Search
  approximate_search_end(search,matches);
}
