/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_filtering_stages.h"
#include "approximate_search_control.h"
#include "approximate_search_region_profile.h"
#include "approximate_search_generate_candidates.h"
#include "approximate_search_verify_candidates.h"
#include "region_profile.h"
#include "region_profile_adaptive.h"
#include "region_profile_schedule.h"
#include "filtering_candidates_process.h"
#include "filtering_candidates_verify.h"
#include "filtering_candidates_align.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Exact Filtering Adaptive
 */
void approximate_search_exact_filtering_adaptive(
    approximate_search_t* const search,
    const approximate_search_region_profile_strategy_t profiling_strategy,
    matches_t* const matches) {
  // Region-Minimal Profile (Reduce the number of candidates per region and maximize number of regions)
  approximate_search_region_profile_adaptive(search,profiling_strategy,search->mm_stack);
  if (search->processing_state==asearch_processing_state_no_regions) return; // Corner cases
  if (search->processing_state==asearch_processing_state_exact_matches) return; // Corner cases
  // Generate exact-candidates
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const search_parameters = actual_parameters->search_parameters;
  region_profile_schedule_filtering_fixed(&search->region_profile,ALL,
      REGION_FILTER_DEGREE_ZERO,search_parameters->filtering_threshold);
  approximate_search_generate_exact_candidates(search,matches);
  // Process candidates (just prepare to verification)
  const bool verify_candidates = (matches != NULL);
  filtering_candidates_process_candidates(search->filtering_candidates,
      search->archive,&search->pattern,actual_parameters,verify_candidates,search->mm_stack);
  // Verify Candidates (if needed)
  if (verify_candidates) {
    approximate_search_verify_candidates(search,matches);
    search->processing_state = asearch_processing_state_candidates_verified;
  } else {
    search->processing_state = asearch_processing_state_candidates_processed;
  }
}
void approximate_search_exact_filtering_adaptive_lightweight(
    approximate_search_t* const search,matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Exact)\n");
    tab_global_inc();
  }
  approximate_search_exact_filtering_adaptive(search,region_profile_adaptive_lightweight,matches);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  PROFILE_STOP(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
}
void approximate_search_exact_filtering_adaptive_recovery(
    approximate_search_t* const search,matches_t* const matches) {
  PROFILE_START(GP_AS_READ_RECOVERY,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Recovery Adaptive Filtering (Exact)\n");
    tab_global_inc();
  }
  approximate_search_exact_filtering_adaptive(search,region_profile_adaptive_recovery,matches);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  PROFILE_STOP(GP_AS_READ_RECOVERY,PROFILE_LEVEL);
}
void approximate_search_exact_filtering_adaptive_cutoff(
    approximate_search_t* const search,matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Fast Cut-off Exact)\n");
    tab_global_inc();
  }
  // Parameters
  const as_parameters_t* const as_parameters = search->as_parameters;
  const search_parameters_t* const search_parameters = as_parameters->search_parameters;
  archive_t* const archive = search->archive;
  fm_index_t* const fm_index = archive->fm_index;
  pattern_t* const pattern = &search->pattern;
  const bool* const allowed_enc = search_parameters->allowed_enc;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  mm_stack_t* const mm_stack = search->mm_stack;
  // Select proper key
  const uint8_t* key;
  uint64_t key_length;
  if (!archive->text->run_length) {
    key = pattern->key;
    key_length = pattern->key_length;
  } else {
    key = pattern->rl_key;
    key_length = pattern->rl_key_length;
  }
  // Iterate process of region-candidates-verification
  const region_profile_model_t* const profile_model = &search_parameters->rp_minimal;
  region_profile_generator_t generator;
  region_profile_generator_init(&generator,region_profile,fm_index,key,key_length,allowed_enc,false);
  region_profile->errors_allowed = 0;
  while (region_profile_generator_next_region(region_profile,&generator,profile_model)) {
    // Generate candidates for the last region found
    PROFILE_START(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
    region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
    filtering_candidates_add_interval(filtering_candidates,last_region->lo,
        last_region->hi,last_region->begin,last_region->end,0,mm_stack);
    PROFILE_STOP(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
    // Verify candidates
    PROFILE_START(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING,PROFILE_LEVEL);
    filtering_candidates_process_candidates(filtering_candidates,archive,pattern,as_parameters,true,mm_stack);
    filtering_candidates_verify_candidates(filtering_candidates,archive,
        search->text_collection,pattern,as_parameters,matches,mm_stack);
    filtering_candidates_align_candidates(
        filtering_candidates,archive->text,archive->locator,search->text_collection,
        pattern,search->emulated_rc_search,as_parameters,false,matches,mm_stack);
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
  } else if (region_profile_has_exact_matches(region_profile)) {
    const region_search_t* const first_region = region_profile->filtering_region;
    search->hi_exact_matches = first_region->hi;
    search->lo_exact_matches = first_region->lo;
    approximate_search_update_mcs(search,1);
    search->processing_state = asearch_processing_state_exact_matches;
  } else {
    // Update MCS (maximum complete stratum)
    search->max_matches_reached = matches_get_num_match_traces(matches) >= search_parameters->search_max_matches;
    if (search->max_matches_reached) approximate_search_update_mcs(search,0);
    // Next State
    search->processing_state = asearch_processing_state_candidates_verified;
  }
  PROFILE_STOP(GP_AS_FILTERING_EXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Boost Exact Filtering Adaptive
 */
void approximate_search_exact_filtering_boost(approximate_search_t* const search,matches_t* const matches) {
#ifdef GEM_PROFILE
  PROFILE_START(GP_AS_FILTERING_EXACT_BOOST,PROFILE_LEVEL);
  const bool already_mapped = matches_is_mapped(matches);
#endif
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Boost)\n");
    tab_global_inc();
  }
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  const uint64_t old_mcs = search->max_complete_stratum;
  // Region-Boost Profile
  approximate_search_region_profile_adaptive(search,region_profile_adaptive_boost,search->mm_stack);
  if (search->processing_state == asearch_processing_state_no_regions ||
      search->processing_state == asearch_processing_state_exact_matches) {
    PROFILE_STOP(GP_AS_FILTERING_EXACT_BOOST,PROFILE_LEVEL);
    gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
    return;
  }
  // Generate exact-candidates
  region_profile_schedule_filtering_fixed(&search->region_profile,
      ALL,REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
  approximate_search_generate_exact_candidates(search,matches);
  // Process candidates (just prepare to verification)
  filtering_candidates_process_candidates(search->filtering_candidates,
      search->archive,&search->pattern,actual_parameters,true,search->mm_stack);
  // Verify candidates
  approximate_search_verify_candidates(search,matches);
  if (!search->max_matches_reached && old_mcs > search->max_complete_stratum) {
    approximate_search_update_mcs(search,old_mcs);
  }
  // DEBUG
#ifdef GEM_PROFILE
  if (!already_mapped) PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_BOOST_MAPPED,matches_is_mapped(matches)?1:0);
  PROFILE_STOP(GP_AS_FILTERING_EXACT_BOOST,PROFILE_LEVEL);
#endif
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Inexact Filtering Adaptive
 */
void approximate_search_inexact_filtering(approximate_search_t* const search,matches_t* const matches) {
  PROFILE_START(GP_AS_FILTERING_INEXACT,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Inexact)\n");
    tab_global_inc();
  }
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
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
  if (search->processing_state == asearch_processing_state_no_regions ||
      search->processing_state == asearch_processing_state_exact_matches) {
    PROFILE_STOP(GP_AS_FILTERING_INEXACT,PROFILE_LEVEL);
    gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
    return;
  }
  // Generate exact-candidates (Dynamic filtering incorporated)
  const uint64_t proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const uint64_t sensibility_misms_length = parameters->filtering_region_factor*proper_length;
  region_profile_schedule_filtering_adaptive(region_profile,search->max_complete_error,sensibility_misms_length);
  approximate_search_generate_inexact_candidates(search,true,true,matches);
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards); // Update MCS
  // Process candidates (just prepare to verification)
  filtering_candidates_process_candidates(search->filtering_candidates,
      search->archive,pattern,actual_parameters,true,search->mm_stack);
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
void approximate_search_verify(approximate_search_t* const search,matches_t* const matches) {
  // Verify
  const uint64_t num_accepted_regions = filtering_candidates_verify_candidates(
      search->filtering_candidates,search->archive,search->text_collection,
      &search->pattern,search->as_parameters,matches,search->mm_stack);
  if (num_accepted_regions > 0) {
    // Realign
    filtering_candidates_align_candidates(search->filtering_candidates,
        search->archive->text,search->archive->locator,
        search->text_collection,&search->pattern,search->emulated_rc_search,
        search->as_parameters,true,matches,search->mm_stack);
  }
  // Update state
  search->processing_state = asearch_processing_state_candidates_verified;
}
/*
 * Unbound Filtering (+ realign)
 */
void approximate_search_unbounded_align(approximate_search_t* const search,matches_t* const matches) {
#ifdef GEM_PROFILE
  PROFILE_START(GP_AS_FILTERING_UNBOUNDED_ALIGN,PROFILE_LEVEL);
  const bool already_mapped = matches_is_mapped(matches);
#endif
  // DEBUG
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Local-Alignment Filtering (unbounded align)\n");
    tab_global_inc();
  }
  // Unbounded-Align discarded candidates
  filtering_candidates_align_unbounded(search->filtering_candidates,
      search->archive->text,search->archive->locator,search->text_collection,
      &search->pattern,search->emulated_rc_search,search->as_parameters,
      matches,search->mm_stack);
  // DEBUG
#ifdef GEM_PROFILE
  if (!already_mapped) PROF_ADD_COUNTER(GP_AS_FILTERING_UNBOUNDED_ALIGN_MAPPED,matches_is_mapped(matches)?1:0);
  PROFILE_STOP(GP_AS_FILTERING_UNBOUNDED_ALIGN,PROFILE_LEVEL);
#endif
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * End of the search
 */
void approximate_search_end(approximate_search_t* const search,matches_t* const matches) {
  // DEBUG
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Search END\n");
  }
  // Check final state
  pattern_t* const pattern = &search->pattern;
  if (search->processing_state == asearch_processing_state_no_regions) {
    gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,"[GEM]>ASM::No-Regions\n"); }
    approximate_search_update_mcs(search,pattern->num_wildcards);
  } else if (search->processing_state == asearch_processing_state_exact_matches) {
    gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,"[GEM]>ASM::Exact-Matches\n"); }
    matches_add_interval_match(matches,search->lo_exact_matches,
        search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search); // Add interval
    approximate_search_update_mcs(search,1);
  }
  // Set MCS
  matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->max_complete_stratum);
  PROF_ADD_COUNTER(GP_AS_ADAPTIVE_MCS,search->max_complete_stratum);
}
