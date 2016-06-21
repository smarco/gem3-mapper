/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_filtering_stages.h"
#include "approximate_search/approximate_search_control.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_generate_candidates.h"
#include "approximate_search/approximate_search_verify_candidates.h"
#include "filtering/region_profile.h"
#include "filtering/region_profile_adaptive.h"
#include "filtering/region_profile_schedule.h"
#include "filtering/filtering_candidates_process.h"
#include "filtering/filtering_candidates_verify.h"
#include "filtering/filtering_candidates_align.h"
#include "filtering/filtering_candidates_align_local.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Exact Filtering Adaptive
 */
void approximate_search_exact_filtering_adaptive(
    approximate_search_t* const search,
    const region_profile_strategy_t strategy,
    matches_t* const matches) {
  // Region-Minimal Profile (Reduce the number of candidates per region and maximize number of regions)
  approximate_search_region_profile_adaptive(search,strategy,search->mm_stack);
  if (search->processing_state==asearch_processing_state_no_regions) return; // Corner case
  // Generate candidates
  region_profile_schedule_filtering_exact(&search->region_profile);
  // Verify Candidates (if needed)
  const bool verify_candidates = (matches != NULL);
  if (verify_candidates) {
    // Generate candidates
    approximate_search_generate_candidates_exact(search,matches);
    // Process candidates
    filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern);
    // Verify candidates
    approximate_search_verify_candidates(search,matches);
    search->processing_state = asearch_processing_state_candidates_verified;
  } else {
    search->processing_state = asearch_processing_state_candidates_processed;
  }
}
void approximate_search_exact_filtering_adaptive_lightweight(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Exact)\n");
    tab_global_inc();
    tab_global_inc();
  }
  approximate_search_exact_filtering_adaptive(search,region_profile_adaptive,matches);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); tab_global_dec(); }
  PROFILE_STOP(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
}
void approximate_search_exact_filtering_adaptive_cutoff(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Fast Cut-off Exact)\n");
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  archive_t* const archive = search->archive;
  fm_index_t* const fm_index = archive->fm_index;
  pattern_t* const pattern = &search->pattern;
  const bool* const allowed_enc = search_parameters->allowed_enc;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  // Select Key (Regular/RL)
  uint8_t* key;
  uint64_t key_length;
  if (pattern->run_length) {
    key = pattern->rl_key;
    key_length = pattern->rl_key_length;
  } else {
    key = pattern->key;
    key_length = pattern->key_length;
  }
  // Iterate process of region-candidates-verification
  const region_profile_model_t* const profile_model = &search_parameters->region_profile_model;
  region_profile_generator_t generator;
  region_profile_generator_init(&generator,region_profile,fm_index,key,key_length,allowed_enc,false);
  region_profile->num_filtered_regions = 0;
  while (region_profile_generator_next_region(region_profile,&generator,profile_model)) {
    // Cut-off
    if (region_profile->num_filtering_regions >= region_profile->max_regions_allocated) break;
    // Generate candidates for the last region found
    PROFILE_START(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
    region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
    filtering_candidates_add_region_interval(
        filtering_candidates,search_parameters,pattern,
        last_region->lo,last_region->hi,last_region->begin,last_region->end,0);
    PROFILE_STOP(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
    // Verify candidates
    PROFILE_START(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING,PROFILE_LEVEL);
    filtering_candidates_process_candidates(filtering_candidates,pattern);
    filtering_candidates_verify_candidates(filtering_candidates,pattern);
    filtering_candidates_align_candidates(filtering_candidates,pattern,false,false,matches);
    PROFILE_STOP(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING,PROFILE_LEVEL);
    // Cut-off condition
    ++(region_profile->num_filtered_regions);
    approximate_search_update_mcs(search,region_profile->num_filtered_regions + pattern->num_wildcards);
    if (asearch_control_fulfilled(search,matches)) break;
  }
  // Check region-profile status
  if (region_profile->num_filtering_regions==0) {
    approximate_search_update_mcs(search,pattern->num_wildcards);
    search->processing_state = asearch_processing_state_no_regions;
  } else {
    search->processing_state = asearch_processing_state_candidates_verified;
  }
  PROFILE_STOP(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Filtering Verification (+ realign)
 */
void approximate_search_verify(
    approximate_search_t* const search,
    matches_t* const matches) {
  // Verify
  const uint64_t num_accepted_regions =
      filtering_candidates_verify_candidates(search->filtering_candidates,&search->pattern);
  // Realign
  if (num_accepted_regions > 0) {
    filtering_candidates_align_candidates(
        search->filtering_candidates,&search->pattern,false,false,matches);
  }
  // Update state
  search->processing_state = asearch_processing_state_candidates_verified;
}
/*
 * Unbound Filtering (+ realign)
 */
void approximate_search_align_local(
    approximate_search_t* const search,
    matches_t* const matches) {
  // DEBUG
  PROFILE_START(GP_AS_FILTERING_LOCAL_ALIGN,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Local-Alignment Filtering (local align)\n");
    tab_global_inc();
  }
  // Align-local discarded candidates
  filtering_candidates_align_local(search->filtering_candidates,&search->pattern,matches);
  // DEBUG
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  PROFILE_STOP(GP_AS_FILTERING_LOCAL_ALIGN,PROFILE_LEVEL);
}
/*
 * End of the search
 */
void approximate_search_end(
    approximate_search_t* const search,
    matches_t* const matches) {
  // DEBUG
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Search END\n");
  }
  // Check final state
  pattern_t* const pattern = &search->pattern;
  if (search->processing_state == asearch_processing_state_no_regions) {
    gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,"[GEM]>ASM::No-Regions\n"); }
    approximate_search_update_mcs(search,pattern->num_wildcards);
  }
  // Set MCS
  PROF_ADD_COUNTER(GP_AS_ADAPTIVE_MCS,search->current_max_complete_stratum);
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->current_max_complete_stratum);
}
