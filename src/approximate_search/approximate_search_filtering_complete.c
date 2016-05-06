/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_complete.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_filtering_complete.h"
#include "approximate_search/approximate_search_neighborhood.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_generate_candidates.h"
#include "filtering/region_profile_schedule.h"
#include "filtering/filtering_candidates_process.h"
#include "filtering/filtering_candidates_verify.h"
#include "filtering/filtering_candidates_align.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Approximate Search based on Filtering Complete Search
 */
void approximate_search_filtering_complete(
    approximate_search_t* const search,
    matches_t* const matches) {
  // PROFILE_START(GP_AS_FILTERING_EXACT,PROFILE_LEVEL); // TODO
  // Parameters
  search_parameters_t* const parameters = search->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Exact Search
  if (search->max_complete_error==0) {
    approximate_search_neighborhood_exact_search(search,matches);
    return;
  }
  // Compute the region profile
  approximate_search_region_profile_adaptive(search,region_profile_adaptive_limited,search->mm_stack);
  if (search->processing_state==asearch_processing_state_no_regions) return;
  // Generate exact-candidates
  parameters->filtering_threshold = ALL; // No restriction
  region_profile_schedule_filtering_fixed(region_profile,ALL,REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
  approximate_search_generate_candidates_exact(search,matches);
  // Verify candidates
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_process_candidates(filtering_candidates,pattern);
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  // Align candidates
  filtering_candidates_align_candidates(filtering_candidates,
      pattern,search->emulated_rc_search,false,false,matches);
  // Update MCS (maximum complete stratum)
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
  // PROFILE_STOP(GP_AS_FILTERING_EXACT); // TODO
}

