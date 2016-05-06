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
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  // Region-Minimal Profile (Reduce the number of candidates per region and maximize number of regions)
  approximate_search_region_profile_adaptive(search,strategy,search->mm_stack);
  if (search->processing_state==asearch_processing_state_no_regions) return; // Corner case
  // Generate candidates
  region_profile_schedule_filtering_fixed(&search->region_profile,ALL,
      REGION_FILTER_DEGREE_ZERO,search_parameters->filtering_threshold);
  approximate_search_generate_candidates_exact(search,matches);
  // Process candidates (just prepare to verification)
  filtering_candidates_process_candidates(search->filtering_candidates,&search->pattern);
  // Verify Candidates (if needed)
  const bool verify_candidates = (matches != NULL);
  if (verify_candidates) {
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
  approximate_search_exact_filtering_adaptive(search,region_profile_adaptive_lightweight,matches);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); tab_global_dec(); }
  PROFILE_STOP(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
}
void approximate_search_exact_filtering_adaptive_heavyweight(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Exact)\n");
    tab_global_inc();
  }
  approximate_search_exact_filtering_adaptive(search,region_profile_adaptive_heavyweight,matches);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
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
  const region_profile_model_t* const profile_model = &search_parameters->rp_lightweight;
  region_profile_generator_t generator;
  region_profile_generator_init(&generator,region_profile,fm_index,key,key_length,allowed_enc,false);
  region_profile->errors_allowed = 0;
  while (region_profile_generator_next_region(region_profile,&generator,profile_model)) {
    // Generate candidates for the last region found
    PROFILE_START(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
    region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
    filtering_candidates_add_region_interval(filtering_candidates,
        last_region->lo,last_region->hi,last_region->begin,last_region->end,0);
    PROFILE_STOP(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
    // Verify candidates
    PROFILE_START(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING,PROFILE_LEVEL);
    filtering_candidates_process_candidates(filtering_candidates,pattern);
    filtering_candidates_verify_candidates(filtering_candidates,pattern);
    filtering_candidates_align_candidates(filtering_candidates,
        pattern,search->emulated_rc_search,false,false,matches);
    asearch_control_adjust_max_differences_using_strata(search,matches);
    PROFILE_STOP(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING,PROFILE_LEVEL);
    // Cut-off condition
    ++(region_profile->errors_allowed);
    approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
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
 * Inexact Filtering Adaptive
 */
void approximate_search_inexact_filtering(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_INEXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Inexact)\n");
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const parameters = search->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Region-Delimit Profile (Maximize number of unique-regions and try to isolate standard/repetitive regions)
  approximate_search_region_profile_adaptive(search,region_profile_adaptive_lightweight,search->mm_stack);
  if (search->region_profile.num_filtering_regions <= 1) {
    search->processing_state = asearch_processing_state_no_regions;
    PROFILE_STOP(GP_AS_FILTERING_INEXACT,PROFILE_LEVEL);
    gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
    return;
  }
  if (search->processing_state == asearch_processing_state_no_regions) {
    PROFILE_STOP(GP_AS_FILTERING_INEXACT,PROFILE_LEVEL);
    gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
    return;
  }
  // Generate exact-candidates (Dynamic filtering incorporated)
  const uint64_t proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const uint64_t sensibility_misms_length = parameters->filtering_region_factor*proper_length;
  region_profile_schedule_filtering_adaptive(region_profile,search->max_complete_error,sensibility_misms_length);
  approximate_search_generate_candidates_inexact(search,true,true,matches);
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards); // Update MCS
  // Process candidates (just prepare to verification)
  filtering_candidates_process_candidates(search->filtering_candidates,pattern);
  // Verify candidates
  approximate_search_verify_candidates(search,matches);
  // Stats -> TODO to control
//  PROF_ADD_COUNTER(GP_AS_FILTERING_INEXACT_MAPPED,matches_is_mapped(matches)?1:0);
//  PROF_ADD_COUNTER(GP_AS_FILTERING_INEXACT_MCS,search->max_complete_stratum);
  PROFILE_STOP(GP_AS_FILTERING_INEXACT,PROFILE_LEVEL);
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
  if (num_accepted_regions > 0) {
    // Realign
    filtering_candidates_align_candidates(search->filtering_candidates,
        &search->pattern,search->emulated_rc_search,false,false,matches);
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
  filtering_candidates_align_local(search->filtering_candidates,
      &search->pattern,search->emulated_rc_search,matches);
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
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->max_complete_stratum);
  PROF_ADD_COUNTER(GP_AS_ADAPTIVE_MCS,search->max_complete_stratum);
}
