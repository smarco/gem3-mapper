/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_base.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_filtering_base.h"
#include "approximate_search_filtering_control.h"
#include "search_parameters.h"

#include "neighborhood_search.h"
#include "region_profile.h"
#include "region_profile_adaptive.h"
#include "region_profile_boost.h"
#include "region_profile_fixed.h"
#include "region_profile_schedule.h"

/*
 * Debug
 */
#define DEBUG_REGION_PROFILE_PRINT  true //GEM_DEEP_DEBUG
#define DEBUG_REGION_SCHEDULE_PRINT GEM_DEEP_DEBUG

/*
 * Constants
 */
#define ZERO_ERRORS 0
#define ONE_ERROR   1
#define TWO_ERRORS  2

/*
 * Stats
 */
#ifndef GEM_NOPROFILE
GEM_INLINE void region_profile_lightweight_stats(region_profile_t* const region_profile) {
  PROF_INC_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT);
  uint64_t total_candidates = 0;
  PROF_ADD_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS_STANDARD,region_profile->num_standard_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS_UNIQUE,
      region_profile->num_filtering_regions-region_profile->num_standard_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_REGION_LENGTH,region->end-region->begin);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_TOTAL_CANDIDATES,total_candidates);
}
GEM_INLINE void region_profile_boost_stats(region_profile_t* const region_profile) {
  PROF_INC_COUNTER(GP_REGION_PROFILE_BOOST);
  uint64_t total_candidates = 0;
  PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS_STANDARD,region_profile->num_standard_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS_UNIQUE,
      region_profile->num_filtering_regions-region_profile->num_standard_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_REGION_LENGTH,region->end-region->begin);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_TOTAL_CANDIDATES,total_candidates);
}
GEM_INLINE void region_profile_delimit_stats(region_profile_t* const region_profile) {
  PROF_INC_COUNTER(GP_REGION_PROFILE_DELIMIT);
  uint64_t total_candidates = 0;
  PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS_STANDARD,region_profile->num_standard_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS_UNIQUE,
      region_profile->num_filtering_regions-region_profile->num_standard_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_REGION_LENGTH,region->end-region->begin);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_TOTAL_CANDIDATES,total_candidates);
}
#endif
/*
 * Adaptive region-profile generation
 */
GEM_INLINE void approximate_search_generate_region_profile_fixed(
    approximate_search_t* const search,mm_stack_t* const mm_stack) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  const archive_t* archive = search->archive;
  fm_index_t* const fm_index = search->archive->fm_index;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
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
  // Schedule the region profile
  const uint64_t proper_length = fm_index_get_proper_length(fm_index);
  const uint64_t region_length = proper_length; // FIXME tune
  region_profile_generate_fixed_schedule(region_profile,key,key_length,parameters->allowed_enc,region_length);
  // Search the region profile
  region_profile_generate_fixed_query(region_profile,fm_index,key);
  // DEBUG
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) {
    region_profile_print(stderr,region_profile,false);
    region_profile_print_fixed_regions(stdout,region_profile,key);
  }
  // Check Zero-Region & Exact-Matches
  if (region_profile->num_filtering_regions==0) {
    approximate_search_update_mcs(search,pattern->num_wildcards);
    search->search_state = asearch_no_regions;
    return;
  } else if (region_profile_has_exact_matches(region_profile)) {
    const region_search_t* const first_region = region_profile->filtering_region;
    search->hi_exact_matches = first_region->hi;
    search->lo_exact_matches = first_region->lo;
    approximate_search_update_mcs(search,1);
    search->search_state = asearch_exact_matches;
    return;
  }
  // Stats
#ifndef GEM_NOPROFILE
  // region_profile_fixed_stats(region_profile); // TODO
#endif
}
/*
 * Adaptive region-profile generation
 */
GEM_INLINE void approximate_search_generate_region_profile_adaptive(
    approximate_search_t* const search,const region_profiling_strategy_t strategy,mm_stack_t* const mm_stack) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  const archive_t* archive = search->archive;
  fm_index_t* const fm_index = search->archive->fm_index;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
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
  // Compute the region profile
  switch (strategy) {
    case region_profile_adaptive_lightweight:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_minimal,ALL,false);
      break;
    case region_profile_adaptive_boost:
      region_profile_generate_adaptive_boost(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_boost,mm_stack);
      break;
    case region_profile_adaptive_limited:
      region_profile_generate_adaptive_limited(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_minimal,search->max_complete_error+1);
      break;
    case region_profile_adaptive_delimit:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_delimit,UINT64_MAX,true);
      break;
    case region_profile_adaptive_recovery:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_boost,UINT64_MAX,true);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) { region_profile_print(stderr,region_profile,false); }
  // Check Zero-Region & Exact-Matches
  if (region_profile->num_filtering_regions==0) {
    approximate_search_update_mcs(search,pattern->num_wildcards);
    search->search_state = asearch_no_regions;
    return;
  } else if (region_profile_has_exact_matches(region_profile)) {
    const region_search_t* const first_region = region_profile->filtering_region;
    search->hi_exact_matches = first_region->hi;
    search->lo_exact_matches = first_region->lo;
    approximate_search_update_mcs(search,1);
    search->search_state = asearch_exact_matches;
    return;
  }
  // Stats
#ifndef GEM_NOPROFILE
  switch (strategy) {
    case region_profile_adaptive_lightweight:
      region_profile_lightweight_stats(region_profile);
      break;
    case region_profile_adaptive_boost:
      region_profile_boost_stats(region_profile);
      break;
    case region_profile_adaptive_delimit:
      region_profile_delimit_stats(region_profile);
      break;
    case region_profile_adaptive_limited:
    case region_profile_adaptive_recovery:
      region_profile_lightweight_stats(region_profile);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
#endif
}
/*
 * Generate Candidates
 */
GEM_INLINE uint64_t approximate_search_generate_region_candidates(
    fm_index_t* const fm_index,uint8_t* const key,region_search_t* const region,
    const uint64_t filtering_threshold,filtering_candidates_t* const filtering_candidates,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  bool perform_search = true;
  uint64_t candidates;
  // Filter up to n errors (n>=2)
  while (region->min >= REGION_FILTER_DEGREE_TWO) {
    PROF_START(GP_AS_GENERATE_CANDIDATES_SEARCH_D2);
    PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2);
    if (perform_search) {
      interval_set_clear(intervals_result);
      neighborhood_search(fm_index,key+region->begin,region->end-region->begin,region->min-1,intervals_result,mm_stack);
      perform_search = false;
    }
    PROF_STOP(GP_AS_GENERATE_CANDIDATES_SEARCH_D2);
    candidates = interval_set_count_intervals_length(intervals_result);
    if (candidates <= filtering_threshold) {
      PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2_HIT);
      PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2_HIT_CANDIDATES,candidates);
      filtering_candidates_add_interval_set(filtering_candidates,intervals_result,region->begin,region->end,mm_stack);
      return region->min; // Return filtered-degree (errors-allowed)
    }
    --(region->min);
  }
  // Filter up to 1 errors
  if (region->min == REGION_FILTER_DEGREE_ONE) {
    PROF_START(GP_AS_GENERATE_CANDIDATES_SEARCH_D1);
    PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1);
    if (perform_search) {
      interval_set_clear(intervals_result);
      neighborhood_search(fm_index,key+region->begin,region->end-region->begin,ONE_ERROR,intervals_result,mm_stack);
      perform_search = false;
    }
    PROF_STOP(GP_AS_GENERATE_CANDIDATES_SEARCH_D1);
    candidates = interval_set_count_intervals_length_thresholded(intervals_result,ONE_ERROR);
    if (candidates <= filtering_threshold) {
      PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1_HIT);
      PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1_HIT_CANDIDATES,candidates);
      filtering_candidates_add_interval_set_thresholded(filtering_candidates,
          intervals_result,region->begin,region->end,ONE_ERROR,mm_stack);
      return REGION_FILTER_DEGREE_ONE; // Return filtered-degree (errors-allowed)
    }
    --(region->min);
  }
  // Filter exact
  if (region->min == REGION_FILTER_DEGREE_ZERO) {
    PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT);
    PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT_CANDIDATES,region->hi-region->lo);
    filtering_candidates_add_interval(filtering_candidates,
        region->lo,region->hi,region->begin,region->end,ZERO_ERRORS,mm_stack);
    return REGION_FILTER_DEGREE_ZERO; // Return filtered-degree (errors-allowed)
  }
  /*
   * Otherwise, ignore the region
   */
  return 0; // Return filtered-degree (errors-allowed)
}
GEM_INLINE void approximate_search_generate_exact_candidates(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_GENERATE_CANDIDATES);
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  mm_stack_t* const mm_stack = search->mm_stack;
  // Pattern
  pattern_t* const pattern = &search->pattern;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // Region profile
  const uint64_t num_regions = region_profile->num_filtering_regions;
  uint64_t errors_allowed = 0; // Number of errors allowed/generated/applied so far
  // Generate candidates for each region
  PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS,num_regions);
  REGION_LOCATOR_ITERATE(region_profile,region,position) {
    PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_PROCESSED);
    // Generate exact-candidates for the region
    if (region->min == REGION_FILTER_DEGREE_ZERO) {
      PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT);
      PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT_CANDIDATES,region->hi-region->lo);
      filtering_candidates_add_interval(filtering_candidates,
          region->lo,region->hi,region->begin,region->end,ZERO_ERRORS,mm_stack);
      ++errors_allowed; // Increase errors-allowed
    }
  }
  // Set the minimum number of mismatches required
  approximate_search_update_mcs(search,errors_allowed + num_wildcards);
  region_profile->errors_allowed = errors_allowed;
  PROF_STOP(GP_AS_GENERATE_CANDIDATES);
}
GEM_INLINE void approximate_search_generate_inexact_candidates(
    approximate_search_t* const search,const bool dynamic_scheduling,
    const bool verify_ahead,matches_t* const matches) {
  PROF_START(GP_AS_GENERATE_CANDIDATES);
  /*
   * Filters all the regions up to the scheduled degree
   *  - Dynamic Scheduling: Assigns a filtering degree to each region as the search goes on
   *  - Verify Ahead: Verifies candidates each queried region. Thus, it can reduce the scope of the search
   */
  // Parameters
  const as_parameters_t* const as_parameters = search->as_parameters;
  const search_parameters_t* const parameters = as_parameters->search_parameters;
  const uint64_t filtering_threshold = parameters->filtering_threshold;
  const double proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const uint64_t sensibility_error_length = parameters->filtering_region_factor*proper_length;
  // Approximate Search Data
  fm_index_t* const fm_index = search->archive->fm_index;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  interval_set_t* const intervals_result = search->interval_set;
  mm_stack_t* const mm_stack = search->mm_stack;
  // Pattern
  pattern_t* const pattern = &search->pattern;
  const uint64_t num_wildcards = pattern->num_wildcards;
  uint8_t* const key = search->pattern.key;
  // Region profile
  const uint64_t num_regions = region_profile->num_filtering_regions;
  const uint64_t num_standard_regions = region_profile->num_standard_regions;
  const uint64_t num_unique_regions = num_regions - region_profile->num_standard_regions;
  uint64_t num_standard_regions_left = num_standard_regions;
  uint64_t num_unique_regions_left = num_unique_regions;
  uint64_t errors_allowed = 0; // Number of errors allowed/generated/applied so far
  // Dynamic Scheduling (Pre-Sort regions)
  if (dynamic_scheduling) {
    region_profile_sort_by_estimated_mappability(&search->region_profile);
    gem_cond_debug_block(DEBUG_REGION_SCHEDULE_PRINT) {
      region_profile_schedule_filtering_adaptive(region_profile,search->max_complete_error,sensibility_error_length);
      region_profile_schedule_print(region_profile,search->max_complete_error,sensibility_error_length);
    }
  }
  // Generate candidates for each region
  PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS,MIN(search->max_complete_error+1,num_regions));
  REGION_LOCATOR_ITERATE(region_profile,region,position) {
    PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_PROCESSED);
    // Dynamic Schedule
    if (dynamic_scheduling) {
      if (region->type == region_standard) {
        --num_standard_regions_left;
      } else {
        --num_unique_regions_left;
      }
      region_schedule_filtering_adaptive(
          region,num_standard_regions_left,num_unique_regions_left,
          search->max_complete_error,sensibility_error_length,errors_allowed);
    }
    // Generate candidates for the region
    errors_allowed += approximate_search_generate_region_candidates(fm_index,
        key,region,filtering_threshold,filtering_candidates,intervals_result,mm_stack);
    // Check candidates ahead (Dynamic scheduling+filtering)
    if (verify_ahead && asearch_filter_ahead_candidates(search,matches)) {
      PROF_PAUSE(GP_AS_GENERATE_CANDIDATES);
      PROF_START(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING);
      filtering_candidates_process_candidates(filtering_candidates,search->archive,pattern,as_parameters,true,mm_stack);
      filtering_candidates_verify_candidates(filtering_candidates,search->archive,
          search->text_collection,pattern,as_parameters,matches,mm_stack);
      filtering_candidates_align_candidates(
          filtering_candidates,search->archive->text,search->archive->locator,
          search->text_collection,pattern,search->emulated_rc_search,as_parameters,false,matches,mm_stack);
      approximate_search_adjust_max_differences_using_strata(search,matches);
      PROF_STOP(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING);
      PROF_CONTINUE(GP_AS_GENERATE_CANDIDATES);
    }
    // Check cut-off condition
    approximate_search_update_mcs(search,errors_allowed + num_wildcards);
    if (asearch_fulfilled(search,matches)) {
      PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_SKIPPED,num_regions-(position+1));
      break;
    }
  }
  // Set the minimum number of mismatches required
  region_profile->errors_allowed = errors_allowed;
  PROF_STOP(GP_AS_GENERATE_CANDIDATES);
}
/*
 * Verify Candidates
 */
GEM_INLINE void approximate_search_verify_candidates(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  // Verify Candidates
  filtering_candidates_verify_candidates(filtering_candidates,search->archive,
      search->text_collection,pattern,actual_parameters,matches,search->mm_stack);
  filtering_candidates_align_candidates(
      filtering_candidates,search->archive->text,search->archive->locator,search->text_collection,
      pattern,search->emulated_rc_search,actual_parameters,false,matches,search->mm_stack);
  // Adjust max-differences
  approximate_search_adjust_max_differences_using_strata(search,matches);
  // Update MCS (maximum complete stratum)
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  search->max_matches_reached = num_matches >= parameters->search_max_matches;
  if (search->max_matches_reached) approximate_search_update_mcs(search,0);
  // Next State
  search->search_state = asearch_candidates_verified;
}
