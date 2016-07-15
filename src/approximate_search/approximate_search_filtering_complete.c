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
 * Approximate Complete-Search based on filtering
 */
void approximate_search_filtering_complete(
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
  approximate_search_region_profile_adaptive(search,region_profile_adaptive_limited,search->mm_stack);
  if (search->processing_state==asearch_processing_state_no_regions) return;
  // Generate exact-candidates
  region_profile_schedule_filtering_exact(region_profile);
  approximate_search_generate_candidates_exact(search,matches);
  // Verify candidates
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_process_candidates(filtering_candidates,pattern,true);
  filtering_candidates_verify_candidates(filtering_candidates,pattern);
  // Align candidates
  filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
  // Update MCS (maximum complete stratum)
  approximate_search_update_mcs(search,region_profile->num_filtered_regions + pattern->num_wildcards);
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->current_max_complete_stratum);
}
/*
 * Approximate Complete-Search based on filtering+NS-search
 */
void approximate_search_hybrid_complete_search(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Parameters
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
  // Complete Search with NS-Search (preconditioned)
  approximate_search_neighborhood_search_partition_preconditioned(search,matches);
  // Update MCS (maximum complete stratum)
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->current_max_complete_stratum);
}

