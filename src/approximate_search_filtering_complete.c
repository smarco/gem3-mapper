/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_complete.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_filtering_complete.h"
#include "approximate_search_neighborhood.h"
#include "approximate_search_region_profile.h"
#include "approximate_search_generate_candidates.h"
#include "region_profile_schedule.h"
#include "filtering_candidates_process.h"
#include "filtering_candidates_verify.h"
#include "filtering_candidates_align.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Approximate Search based on Filtering Complete Search
 */
void approximate_search_filtering_complete(approximate_search_t* const search,matches_t* const matches) {
  // PROFILE_START(GP_AS_FILTERING_EXACT,PROFILE_LEVEL); // TODO
  // Parameters
  as_parameters_t* const actual_parameters = search->as_parameters;
  search_parameters_t* const parameters = actual_parameters->search_parameters;
  archive_t* const archive = search->archive;
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
  if (search->processing_state==asearch_processing_state_exact_matches) {
    matches_add_interval_match(matches,search->lo_exact_matches,
        search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search);
    return;
  }
  // Generate exact-candidates
  parameters->filtering_threshold = ALL; // No restriction
  region_profile_schedule_filtering_fixed(region_profile,ALL,REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
  approximate_search_generate_exact_candidates(search,matches);
  // Verify candidates
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_process_candidates(filtering_candidates,archive,
      pattern,actual_parameters,true,search->mm_stack);
  filtering_candidates_verify_candidates(
      filtering_candidates,search->archive,search->text_collection,
      pattern,actual_parameters,matches,search->mm_stack);
  // Align candidates
  filtering_candidates_align_candidates(filtering_candidates,
      search->archive->text,search->archive->locator,search->text_collection,pattern,
      search->emulated_rc_search,actual_parameters,false,false,matches,search->mm_stack);
  // Update MCS (maximum complete stratum)
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
  // PROFILE_STOP(GP_AS_FILTERING_EXACT); // TODO
}

