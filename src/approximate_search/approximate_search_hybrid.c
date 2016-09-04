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
#include "matches/matches_classify.h"

/*
 * Control
 */
void as_hybrid_control_begin(approximate_search_t* const search) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Basic Cases\n");
    tab_global_inc();
  }
  // Parameters
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // All characters are wildcards
  if (key_length==num_wildcards || key_length==0) {
    search->search_stage = asearch_stage_end;
    return;
  }
  // Otherwise, go to standard exact filtering
  search->search_stage = asearch_stage_filtering_adaptive;
  PROF_INC_COUNTER(GP_AS_FILTERING_ADATIVE_CALL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  return;
}
void as_hybrid_control_filtering_adaptive_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROF_ADD_COUNTER(GP_AS_FILTERING_ADATIVE_MCS,search->region_profile.num_filtered_regions);
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  const uint64_t mcs = search->region_profile.num_filtered_regions;
  const uint64_t delta = search->search_parameters->complete_strata_after_best_nominal;
  // Select state
  switch (search->processing_state) {
    case asearch_processing_state_no_regions:
      search->current_max_complete_error = delta;
      search->search_stage = asearch_stage_neighborhood;
      break;
    case asearch_processing_state_candidates_verified:
      // Neighborhood Search : Unmapped
      if (!matches_is_mapped(matches)) {
        PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL);
        PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL_UNMAPPED);
        // Adjust max-error
        search->current_max_complete_error = mcs + delta;
        search->search_stage = asearch_stage_neighborhood;
        break;
      }
      // Check match-class & error-reached
      matches_classify(matches);
      if (matches->matches_class==matches_class_tie_perfect ||
          mcs >= search_parameters->complete_search_error_nominal+1) {
        search->search_stage = asearch_stage_end;
        break;
      }
      // Neighborhood Search : Frontier 0:1+0 & Beyond 0:0+0:0:1
      const uint64_t min_edit_distance = matches_metrics_get_min_edit_distance(&matches->metrics);
      if (min_edit_distance+1 >= mcs) {
#ifdef GEM_PROFILE
        PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL);
        if (min_edit_distance+1 == mcs) PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL_MAP_FRONTIER);
        if (min_edit_distance >= mcs) PROF_INC_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL_MAP_INCOMPLETE);
#endif
        search->current_max_complete_error = MIN(search->current_max_complete_error,min_edit_distance+1);
        search->search_stage = asearch_stage_neighborhood;
        break;
      }
      search->search_stage = asearch_stage_end;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Check neighborhood-search & number of wildcards
  if (search->search_stage==asearch_stage_neighborhood &&
      search->pattern.num_wildcards > search->current_max_complete_error) {
    search->search_stage = asearch_stage_end;
  }
}
void as_hybrid_control_neighborhood_next_state(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROF_ADD_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_MCS,search->current_max_complete_stratum);
  search_parameters_t* const search_parameters = search->search_parameters;
  if (search_parameters->local_alignment==local_alignment_never || matches_is_mapped(matches)) {
    search->search_stage = asearch_stage_end;
  } else {  // local_alignment_if_unmapped
    PROF_INC_COUNTER(GP_AS_LOCAL_ALIGN_CALL);
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
        as_hybrid_control_begin(search);
        break;
      case asearch_stage_filtering_adaptive: // Exact-Filtering (Adaptive)
        approximate_search_exact_filtering_adaptive(search,matches);
        as_hybrid_control_filtering_adaptive_next_state(search,matches); // Next State
        break;
      case asearch_stage_filtering_adaptive_finished:
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
