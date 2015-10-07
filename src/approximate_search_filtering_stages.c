/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_stages.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_filtering_stages.h"
#include "approximate_search_filtering_base.h"
#include "approximate_search_filtering_control.h"
#include "region_profile.h"
#include "region_profile_adaptive.h"
#include "region_profile_schedule.h"

/*
 * Basic Cases
 */
GEM_INLINE void approximate_search_adaptive_mapping_basic_cases(approximate_search_t* const search) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Basic Cases\n");
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // All characters are wildcards
  if (key_length==num_wildcards) {
    search->hi_exact_matches = 0;
    search->lo_exact_matches = 0;
    approximate_search_update_mcs(search,key_length);
    search->search_state = asearch_end;
    return;
  }
  /*
   * Recovery
   *   If the number of wildcards (or errors required) is too high we try to recover as
   *   many matches as possible. We extract feasible regions from the read and filter
   *   them trying to recover anything out of bad quality reads.
   */
  if (num_wildcards > 0) {
    switch (search_parameters->mapping_mode) {
      case mapping_adaptive_filtering_fast:
      case mapping_adaptive_filtering_thorough:
        if (num_wildcards >= search->as_parameters->alignment_max_error_nominal) {
          search->search_state = asearch_read_recovery;
          return;
        }
        break;
      default:
        if (num_wildcards >= search->max_complete_error) {
          search->search_state = asearch_read_recovery;
          return;
        }
        break;
    }
  }
  // Exact search
  if (search->max_complete_error==0) {
    search->search_state = asearch_neighborhood;
    return;
  }
  // Very short reads (Neighborhood search)
  if (key_length <= RANK_MTABLE_SEARCH_DEPTH || key_length < search->archive->fm_index->proper_length) {
    search->search_state = asearch_neighborhood;
    return;
  }
  // Otherwise, go to standard exact filtering
  search->search_state = asearch_exact_filtering_adaptive;
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  return;
}
/*
 * Exact Filtering Fixed
 */
GEM_INLINE void approximate_search_exact_filtering_fixed(approximate_search_t* const search) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Fixed)\n");
    tab_global_inc();
  }
  // Region-Fixed Profile (Equal size regions)
  approximate_search_generate_region_profile_fixed(search,search->mm_stack);
//  if (search->search_state==asearch_no_regions || search->search_state==asearch_exact_matches) { // Check corner cases
//    asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
//    return;
//  }
//  // Generate exact-candidates
//  const as_parameters_t* const actual_parameters = search->as_parameters;
//  const search_parameters_t* const search_parameters = actual_parameters->search_parameters;
//  region_profile_schedule_filtering_fixed(&search->region_profile,ALL,
//      REGION_FILTER_DEGREE_ZERO,search_parameters->filtering_threshold);
//  approximate_search_generate_exact_candidates(search,matches);
//  // Process candidates (just prepare to verification)
//  filtering_candidates_process_candidates(search->filtering_candidates,
//      search->archive,&search->pattern,actual_parameters,verify_candidates,search->mm_stack);
//  // Verify Candidates (if needed)
//  if (verify_candidates) {
//    approximate_search_verify_candidates(search,matches);
//  } else {
//    search->search_state = asearch_verify_candidates;
//  }
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Exact Filtering Adaptive
 */
GEM_INLINE void approximate_search_exact_filtering_adaptive(
    approximate_search_t* const search,const region_profiling_strategy_t profiling_strategy,
    const bool verify_candidates,matches_t* const matches) {
  // Region-Minimal Profile (Reduce the number of candidates per region and maximize number of regions)
  approximate_search_generate_region_profile_adaptive(search,profiling_strategy,search->mm_stack);
  if (search->search_state==asearch_no_regions || search->search_state==asearch_exact_matches) { // Check corner cases
    asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
    return;
  }
  // Generate exact-candidates
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const search_parameters = actual_parameters->search_parameters;
  region_profile_schedule_filtering_fixed(&search->region_profile,ALL,
      REGION_FILTER_DEGREE_ZERO,search_parameters->filtering_threshold);
  approximate_search_generate_exact_candidates(search,matches);
  // Process candidates (just prepare to verification)
  filtering_candidates_process_candidates(search->filtering_candidates,
      search->archive,&search->pattern,actual_parameters,verify_candidates,search->mm_stack);
  // Verify Candidates (if needed)
  if (verify_candidates) {
    approximate_search_verify_candidates(search,matches);
  } else {
    search->search_state = asearch_verify_candidates;
  }
}
GEM_INLINE void approximate_search_exact_filtering_adaptive_lightweight(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Adaptive Filtering (Exact)\n");
    tab_global_inc();
  }
  approximate_search_exact_filtering_adaptive(search,region_profile_adaptive_lightweight,verify_candidates,matches);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  PROF_STOP(GP_AS_FILTERING_EXACT);
}
GEM_INLINE void approximate_search_exact_filtering_adaptive_recovery(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  PROF_START(GP_AS_READ_RECOVERY);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Recovery Adaptive Filtering (Exact)\n");
    tab_global_inc();
  }
  approximate_search_exact_filtering_adaptive(search,region_profile_adaptive_recovery,verify_candidates,matches);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
  PROF_STOP(GP_AS_READ_RECOVERY);
}
GEM_INLINE void approximate_search_exact_filtering_adaptive_cutoff(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT);
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
    PROF_START(GP_AS_GENERATE_CANDIDATES);
    region_search_t* const last_region = region_profile->filtering_region + (region_profile->num_filtering_regions-1);
    filtering_candidates_add_interval(filtering_candidates,last_region->lo,
        last_region->hi,last_region->begin,last_region->end,0,mm_stack);
    PROF_STOP(GP_AS_GENERATE_CANDIDATES);
    // Verify candidates
    PROF_START(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING);
    filtering_candidates_process_candidates(filtering_candidates,archive,pattern,as_parameters,true,mm_stack);
    filtering_candidates_verify_candidates(filtering_candidates,archive,
        search->text_collection,pattern,as_parameters,matches,mm_stack);
    filtering_candidates_align_candidates(
        filtering_candidates,archive->text,archive->locator,search->text_collection,
        pattern,search->emulated_rc_search,as_parameters,false,matches,mm_stack);
    approximate_search_adjust_max_differences_using_strata(search,matches);
    PROF_STOP(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING);
    // Cut-off condition
    ++(region_profile->errors_allowed);
    approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
    if (asearch_fulfilled(search,matches)) break;
  }
  // Check region-profile status
  if (region_profile->num_filtering_regions==0) {
    approximate_search_update_mcs(search,pattern->num_wildcards);
    search->search_state = asearch_no_regions;
    asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
  } else if (region_profile_has_exact_matches(region_profile)) {
    const region_search_t* const first_region = region_profile->filtering_region;
    search->hi_exact_matches = first_region->hi;
    search->lo_exact_matches = first_region->lo;
    approximate_search_update_mcs(search,1);
    search->search_state = asearch_exact_matches;
    asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
  } else {
    // Update MCS (maximum complete stratum)
    search->max_matches_reached = matches_get_num_match_traces(matches) >= search_parameters->search_max_matches;
    if (search->max_matches_reached) approximate_search_update_mcs(search,0);
    // Next State
    search->search_state = asearch_candidates_verified;
    asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
  }
  PROF_STOP(GP_AS_FILTERING_EXACT);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Boost Exact Filtering Adaptive
 */
GEM_INLINE void approximate_search_exact_filtering_boost(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT_BOOST);
#ifdef GEM_PROFILE
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
  approximate_search_generate_region_profile_adaptive(search,region_profile_adaptive_boost,search->mm_stack);
  if (search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
    asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
    PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
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
  search->search_state = asearch_exact_filtering_boost; // Correct State
  // Next State
  asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
#ifdef GEM_PROFILE
  if (!already_mapped) PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_BOOST_MAPPED,matches_is_mapped(matches)?1:0);
#endif
  PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Inexact Filtering Adaptive
 */
GEM_INLINE void approximate_search_inexact_filtering(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_INEXACT);
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
  approximate_search_generate_region_profile_adaptive(search,region_profile_adaptive_lightweight,search->mm_stack);
  if (search->region_profile.num_filtering_regions <= 1) {
    asearch_control_next_state(search,asearch_inexact_filtering,matches);
    PROF_STOP(GP_AS_FILTERING_INEXACT);
    gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
    return;
  }
  if (search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
    PROF_STOP(GP_AS_FILTERING_INEXACT);
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
  search->search_state = asearch_inexact_filtering; // Correct State
  // Next State
  asearch_control_next_state(search,asearch_inexact_filtering,matches);
  PROF_ADD_COUNTER(GP_AS_FILTERING_INEXACT_MAPPED,matches_is_mapped(matches)?1:0);
  PROF_ADD_COUNTER(GP_AS_FILTERING_INEXACT_MCS,search->max_complete_stratum);
  PROF_STOP(GP_AS_FILTERING_INEXACT);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Unbound Filtering
 */
GEM_INLINE void approximate_search_unbounded_align(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_UNBOUNDED_ALIGN);
#ifdef GEM_PROFILE
  const bool already_mapped = matches_is_mapped(matches);
#endif
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,"[GEM]>ASM::Local-Alignment Filtering (unbounded align)\n");
    tab_global_inc();
  }
  // Unbounded-Align discarded candidates
  filtering_candidates_align_unbounded(search->filtering_candidates,
      search->archive->text,search->archive->locator,search->text_collection,
      &search->pattern,search->emulated_rc_search,search->as_parameters,
      matches,search->mm_stack);
  // Next State
  asearch_control_next_state(search,asearch_unbounded_alignment,matches);
#ifdef GEM_PROFILE
  if (!already_mapped) PROF_ADD_COUNTER(GP_AS_FILTERING_UNBOUNDED_ALIGN_MAPPED,matches_is_mapped(matches)?1:0);
#endif
  PROF_STOP(GP_AS_FILTERING_UNBOUNDED_ALIGN);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_global_dec(); }
}
/*
 * Filtering Verification (+ realign)
 */
GEM_INLINE void approximate_search_verify(approximate_search_t* const search,matches_t* const matches) {
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
  if (search->search_state==asearch_verify_candidates) {
    search->search_state = asearch_candidates_verified;
  }
}
GEM_INLINE void approximate_search_verify_using_bpm_buffer(
    approximate_search_t* const search,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end,
    matches_t* const matches) {
  // Retrieve
  const uint64_t num_accepted_regions = filtering_candidates_bpm_buffer_retrieve(
      search->filtering_candidates,search->archive->text,search->text_collection,
      &search->pattern,bpm_gpu_buffer,candidate_offset_begin,candidate_offset_end,search->mm_stack);
  if (num_accepted_regions > 0) {
    // Realign
    PROF_START(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS);
    filtering_candidates_align_candidates(search->filtering_candidates,
        search->archive->text,search->archive->locator,
        search->text_collection,&search->pattern,search->emulated_rc_search,
        search->as_parameters,true,matches,search->mm_stack);
    PROF_STOP(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS);
  }
  // Update state
  if (search->search_state==asearch_verify_candidates) {
    search->search_state = asearch_candidates_verified;
  }
}
