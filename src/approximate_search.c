/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   // TODO CSearch
 *     -- fmi_probe_matches
 */

#include "approximate_search.h"

#include "quality_model.h"
#include "rank_mtable.h"
#include "neighborhood_search.h"

#include "mapper_profile.h"

/*
 * Debug
 */
#define DEBUG_REGION_PROFILE_PRINT  GEM_DEEP_DEBUG
#define DEBUG_REGION_SCHEDULE_PRINT GEM_DEEP_DEBUG

/*
 * Explicit search values
 */
// Unique
#define ALL_MAPPINGS UINT64_MAX
#define UNIQUE 1
// Filtering policies
#define DYNAMIC_SCHEDULING true
#define STATIC_SCHEDULING false
#define DYNAMIC_FILTERING true
#define STATIC_FILTERING false
// Mismatches to perform a Neighborhood Search
#define ZERO_ERRORS 0
#define ONE_ERROR   1
#define TWO_ERRORS  2
// Fast-mapping modes
#define FM_NO_FAST_MODE 0
#define FM_ADAPTIVE (UINT64_MAX-1)
/*
 * Workflow return value (control)
 */
#define WF_CONTINUE 0
#define WF_STOP 1

/*
 * Setup
 */
GEM_INLINE void approximate_search_init(
    approximate_search_t* const search,archive_t* const archive,
    as_parameters_t* const as_parameters,const bool emulated_rc_search) {
  // Index Structures & Parameters
  search->archive = archive;
  search->as_parameters = as_parameters;
  search->emulated_rc_search = emulated_rc_search;
}
GEM_INLINE void approximate_search_configure(
    approximate_search_t* const search,filtering_candidates_t* const filtering_candidates,
    text_collection_t* text_collection,interval_set_t* const interval_set,mm_stack_t* const mm_stack) {
  // Set Auxiliary Structures (external)
  search->filtering_candidates = filtering_candidates; // Filtering Candidates
  search->text_collection = text_collection;           // Text-Collection
  search->interval_set = interval_set;                 // Interval Set
  search->mm_stack = mm_stack;                         // Set MM
}
GEM_INLINE void approximate_search_reset(approximate_search_t* const search) {
  // Reset Approximate Search State
  search->search_state = asearch_begin;
  search->verify_candidates = true;
  search->stop_before_neighborhood_search = false;
  const uint64_t max_search_error = search->as_parameters->max_search_error_nominal;
  const uint64_t max_effective_filtering_error = search->pattern.max_effective_filtering_error;
  search->max_differences = MIN(max_search_error,max_effective_filtering_error);
  search->max_complete_stratum = ALL;
  search->max_matches_reached = false;
  // Prepare filtering candidates
  filtering_candidates_clear(search->filtering_candidates);
  // Prepare region profile
  if (search->max_differences > 0) {
    region_profile_new(&search->region_profile,search->pattern.key_length,search->mm_stack);
  }
}
GEM_INLINE void approximate_search_destroy(approximate_search_t* const search) { /* NOP */ }
/*
 * Accessors
 */
GEM_INLINE uint64_t approximate_search_get_num_filtering_candidates(const approximate_search_t* const search) {
  if (search->search_state == asearch_exact_matches) {
    return search->hi_exact_matches - search->lo_exact_matches;
  } else {
    const filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
    return filtering_candidates_get_num_candidate_regions(filtering_candidates);
  }
}
GEM_INLINE uint64_t approximate_search_get_num_exact_filtering_candidates(const approximate_search_t* const search) {
  return (search->search_state == asearch_exact_matches) ? search->hi_exact_matches - search->lo_exact_matches : 0;
}
GEM_INLINE void approximate_search_update_mcs(approximate_search_t* const search,const uint64_t max_complete_stratum) {
  search->max_complete_stratum = max_complete_stratum;
}
/*
 * Stats
 */
GEM_INLINE void approximate_search_generate_region_profile_minimal_stats(approximate_search_t* const search) {
#ifndef GEM_NOPROFILE
  PROF_INC_COUNTER(GP_REGION_PROFILE_MINIMAL);
  region_profile_t* const region_profile = &search->region_profile;
  uint64_t total_candidates = 0;
  PROF_ADD_COUNTER(GP_REGION_PROFILE_MINIMAL_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_MINIMAL_NUM_REGIONS_STANDARD,region_profile->num_standard_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_MINIMAL_NUM_REGIONS_UNIQUE,region_profile->num_filtering_regions-region_profile->num_standard_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_MINIMAL_REGION_LENGTH,region->start-region->end);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_MINIMAL_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_MINIMAL_TOTAL_CANDIDATES,total_candidates);
#endif
}
GEM_INLINE void approximate_search_generate_region_profile_boost_stats(approximate_search_t* const search) {
#ifndef GEM_NOPROFILE
  PROF_INC_COUNTER(GP_REGION_PROFILE_BOOST);
  region_profile_t* const region_profile = &search->region_profile;
  uint64_t total_candidates = 0;
  PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS_STANDARD,region_profile->num_standard_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS_UNIQUE,region_profile->num_filtering_regions-region_profile->num_standard_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_REGION_LENGTH,region->start-region->end);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_BOOST_TOTAL_CANDIDATES,total_candidates);
#endif
}
GEM_INLINE void approximate_search_generate_region_profile_delimit_stats(approximate_search_t* const search) {
#ifndef GEM_NOPROFILE
  PROF_INC_COUNTER(GP_REGION_PROFILE_DELIMIT);
  region_profile_t* const region_profile = &search->region_profile;
  uint64_t total_candidates = 0;
  PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS_STANDARD,region_profile->num_standard_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS_UNIQUE,region_profile->num_filtering_regions-region_profile->num_standard_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_REGION_LENGTH,region->start-region->end);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_DELIMIT_TOTAL_CANDIDATES,total_candidates);
#endif
}
/*
 * Approximate Search Control
 */
GEM_INLINE void approximate_search_adjust_max_differences_using_strata(
    approximate_search_t* const search,matches_t* const matches) {
  const uint64_t max_differences = search->max_differences;
  const uint64_t delta = search->as_parameters->complete_strata_after_best_nominal;
  /*
   * If delta parameter is set (and is below the maximum number of mismatches),
   * finds the minimum non zero stratum (mnzs) and adjusts
   * the maximum number of mismatches to (mnzs+delta)
   */
  if (delta < max_differences) {
    const int64_t fms = matches_counters_get_min_distance(matches);
    if (fms>=0 && fms+delta < max_differences) {
      search->max_differences = fms+delta;
    }
  }
}
// Control Fulfilled-search
GEM_INLINE bool asearch_fulfilled(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  if (matches==NULL) return false;
  // Determines when the search is done following the mapping criteria
  switch (search_parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast: {
      // Unmapped
      if (matches_counters_get_total_count(matches)==0) return false;
      // Ties
      const uint64_t min_distance = matches_counters_get_min_distance(matches);
      if (matches_counters_get_count(matches,min_distance) > 1) return true;
      // Detect signal
      const uint64_t key_length = search->pattern.key_length;
      const swg_penalties_t* const swg_penalties = &search->as_parameters->search_parameters->swg_penalties;
      double pr = matches_classify_ambiguous(matches,search->max_complete_stratum,swg_penalties,key_length);
      if (pr <= 0.98) return false;
      pr = matches_classify_unique(matches,search->max_complete_stratum,swg_penalties,key_length);
      if (pr >= 0.999) return true;
      // Otherwise
      return false;
    }
    case mapping_adaptive_filtering_match: {
      // Unmapped
      if (matches_counters_get_total_count(matches)==0) return false;
      // Ties
      const uint64_t min_distance = matches_counters_get_min_distance(matches);
      if (matches_counters_get_count(matches,min_distance) > 1) return true;
      // Otherwise
      return false;
    }
    case mapping_adaptive_filtering_complete:
      return search->max_complete_stratum > search->max_differences;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return true;
}
// Control Filter-ahead
GEM_INLINE bool asearch_filter_ahead_candidates(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const search_parameters = actual_parameters->search_parameters;
  // Determines when the search is done following the mapping criteria
  switch (search_parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
      return true;
    case mapping_adaptive_filtering_match:
      return true;
    case mapping_adaptive_filtering_complete:
      return actual_parameters->complete_strata_after_best_nominal < search->max_differences;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return true;
}
// Control Filtering-Fast
GEM_INLINE void asearch_control_fast_next_state(
    approximate_search_t* const search,const approximate_search_state_t processing_step,matches_t* const matches) {
  // Filtering-Fast: Best efficiency, keeping good-quality results (search for a "good" match)
  search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  switch (processing_step) {
    case asearch_exact_filtering_adaptive:
      if (search->search_state==asearch_no_regions) return;
      if (search->search_state==asearch_exact_matches) return;
      search->search_state = (matches_get_num_match_traces(matches)==0) ?
          (search_parameters->local_alignment == local_alignment_never) ? asearch_end : asearch_local_alignment:
          (search_parameters->local_alignment == local_alignment_always) ? asearch_local_alignment : asearch_end;
      break;
    case asearch_local_alignment:
      search->search_state = asearch_end;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
// Control Filtering-Match
GEM_INLINE void asearch_control_match_next_state(
    approximate_search_t* const search,
    const approximate_search_state_t processing_step,matches_t* const matches) {
  // Filtering-Match: Spend more computing resources as to guarantee good-quality (increasing likelihood of results)
  search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  switch (processing_step) {
    case asearch_exact_filtering_adaptive:
      if (search->search_state==asearch_exact_matches) return;
      if (search->search_state==asearch_no_regions) {
        search->search_state = asearch_exact_filtering_boost;
        return;
      }
      if (asearch_fulfilled(search,matches)) {
        search->search_state = asearch_end;
      } else {
        search->search_state = asearch_exact_filtering_boost;
      }
      break;
    case asearch_exact_filtering_boost:
      search->search_state = (matches_get_num_match_traces(matches)==0) ?
          (search_parameters->local_alignment == local_alignment_never) ? asearch_end : asearch_local_alignment:
          (search_parameters->local_alignment == local_alignment_always) ? asearch_local_alignment : asearch_end;
      break;
    case asearch_local_alignment:
      search->search_state = asearch_end;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
// Control Complete
GEM_INLINE void asearch_control_complete_next_state(
    approximate_search_t* const search,const approximate_search_state_t processing_step,
    matches_t* const matches) {
  GEM_NOT_IMPLEMENTED();
}
// Control HUB
GEM_INLINE void asearch_control_next_state(
    approximate_search_t* const search,const approximate_search_state_t processing_step,
    matches_t* const matches) {
  // Parameters
  const search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  // Determines when the search is done following the mapping criteria
  switch (search_parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
      asearch_control_fast_next_state(search,processing_step,matches);
      break;
    case mapping_adaptive_filtering_match:
      asearch_control_match_next_state(search,processing_step,matches);
      break;
    case mapping_adaptive_filtering_complete:
      asearch_control_complete_next_state(search,processing_step,matches);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Region Profile
 */
typedef enum {
  region_profile_adaptive_lightweight, // Adaptive Profile stopping whenever any cut-off cond. is reached
  region_profile_adaptive_limited,     // Adaptive Profile limiting the max. length of a region (max-region-length)
  region_profile_adaptive_extensive,   // Adaptive Profile extensive extracting all possible regions
} region_profile_generation_strategy_t;
GEM_INLINE void approximate_search_generate_region_profile(
    approximate_search_t* const search,const region_profile_model_t* const profile_model,
    const region_profile_generation_strategy_t region_profile_generation_strategy,
    const uint64_t min_regions,const bool allow_zero_regions) {
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
  switch (region_profile_generation_strategy) {
    case region_profile_adaptive_lightweight:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,profile_model,search->max_differences+1,allow_zero_regions);
      break;
    case region_profile_adaptive_limited:
      region_profile_generate_adaptive_limited(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,profile_model,min_regions);
      break;
    case region_profile_adaptive_extensive:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,profile_model,UINT64_MAX,allow_zero_regions);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) { region_profile_print(stderr,region_profile,false,false); }
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
}
/*
 * Filtering Regions
 */
GEM_INLINE void approximate_search_add_regions_to_filter(
    filtering_candidates_t* const filtering_candidates,
    const region_profile_t* const region_profile,mm_stack_t* const mm_stack) {
  // Add all candidates from the hi/lo regions of the profile to filtering_candidates
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    const uint64_t num_candidates = region->hi-region->lo;
    if (gem_expect_false(num_candidates==0)) continue;
    // Add region candidates (zero degree)
    filtering_candidates_add_interval(filtering_candidates,
        region->lo,region->hi,region->start,region->end,ZERO_ERRORS,mm_stack);
  }
}
GEM_INLINE void approximate_search_region_profile_generate_candidates(
    approximate_search_t* const search,const bool dynamic_scheduling,
    const bool verify_ahead,const uint64_t sensibility_error_length,
    const uint64_t filtering_threshold,matches_t* const matches) {
  PROF_START(GP_AS_GENERATE_CANDIDATES);
  /*
   * Filters all the regions up to the scheduled degree
   *  - Dynamic Scheduling: Assigns a filtering degree to each region as the search goes on
   *  - Verify Ahead: Verifies candidates each queried region. Thus, it can reduce the scope of the search
   */
  // Data
  const as_parameters_t* const actual_parameters = search->as_parameters;
  mm_stack_t* const mm_stack = search->mm_stack;
  // Approximate Search Data
  fm_index_t* const fm_index = search->archive->fm_index;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  interval_set_t* const intervals_result = search->interval_set;
  // Pattern
  pattern_t* const pattern = &search->pattern;
  const uint64_t num_wildcards = pattern->num_wildcards;
  uint8_t* const key = search->pattern.key;
  // Region profile
  const uint64_t num_regions = region_profile->num_filtering_regions;
  const uint64_t num_standard_regions = region_profile->num_standard_regions;
  const uint64_t num_unique_regions = num_regions - region_profile->num_standard_regions;
  PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS,MIN(search->max_differences+1,num_regions));
  uint64_t num_standard_regions_left = num_standard_regions;
  uint64_t num_unique_regions_left = num_unique_regions;
  uint64_t candidates, errors_allowed = 0; // Number of errors allowed/generated/applied so far
  // Dynamic Scheduling (Pre-Sort regions)
  if (dynamic_scheduling) {
    region_profile_sort_by_estimated_mappability(&search->region_profile);
    gem_cond_debug_block(DEBUG_REGION_SCHEDULE_PRINT) {
      region_profile_schedule_filtering_adaptive(region_profile,search->max_differences,sensibility_error_length);
      region_profile_schedule_print(region_profile,search->max_differences,sensibility_error_length);
    }
  }
  // Generate candidates for each region
  REGION_LOCATOR_ITERATE(region_profile,region,position) {
    PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_PROCESSED);
    bool perform_search = true;
    // Dynamic Schedule
    if (dynamic_scheduling) {
      if (region->type == region_standard) {
        --num_standard_regions_left;
      } else {
        --num_unique_regions_left;
      }
      region_schedule_filtering_adaptive(
          region,num_standard_regions_left,num_unique_regions_left,
          search->max_differences,sensibility_error_length,errors_allowed);
    }
    // Filter up to n errors (n>=2)
    while (region->min >= REGION_FILTER_DEGREE_TWO) {
      PROF_START(GP_AS_GENERATE_CANDIDATES_SEARCH_D2);
      PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2);
      if (perform_search) {
        interval_set_clear(intervals_result);
        neighborhood_search(fm_index,key+region->end,region->start-region->end,region->min-1,intervals_result,mm_stack);
        perform_search = false;
      }
      PROF_STOP(GP_AS_GENERATE_CANDIDATES_SEARCH_D2);
      candidates = interval_set_count_intervals_length(intervals_result);
      if (candidates <= filtering_threshold) {
        PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2_HIT);
        PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2_HIT_CANDIDATES,candidates);
        filtering_candidates_add_interval_set(filtering_candidates,intervals_result,region->start,region->end,mm_stack);
        errors_allowed += region->min;
        break;
      } else {
        --(region->min);
      }
    }
    // Filter up to 1 errors
    if (region->min == REGION_FILTER_DEGREE_ONE) {
      PROF_START(GP_AS_GENERATE_CANDIDATES_SEARCH_D1);
      PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1);
      if (perform_search) {
        interval_set_clear(intervals_result);
        neighborhood_search(fm_index,key+region->end,region->start-region->end,ONE_ERROR,intervals_result,mm_stack);
        perform_search = false;
      }
      PROF_STOP(GP_AS_GENERATE_CANDIDATES_SEARCH_D1);
      candidates = interval_set_count_intervals_length_thresholded(intervals_result,ONE_ERROR);
      if (candidates <= filtering_threshold) {
        PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1_HIT);
        PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1_HIT_CANDIDATES,candidates);
        filtering_candidates_add_interval_set_thresholded(filtering_candidates,
            intervals_result,region->start,region->end,ONE_ERROR,mm_stack);
        errors_allowed += REGION_FILTER_DEGREE_ONE;
      } else {
        --(region->min);
      }
    }
    // Filter exact
    if (region->min == REGION_FILTER_DEGREE_ZERO) {
      PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT);
      PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT_CANDIDATES,region->hi-region->lo);
      filtering_candidates_add_interval(filtering_candidates,
          region->lo,region->hi,region->start,region->end,ZERO_ERRORS,mm_stack);
      errors_allowed += REGION_FILTER_DEGREE_ZERO;
    }
    /*
     * Otherwise, ignore the region
     */
    // Check candidates ahead (Dynamic scheduling+filtering)
    if (verify_ahead && asearch_filter_ahead_candidates(search,matches)) {
      PROF_PAUSE(GP_AS_GENERATE_CANDIDATES);
      PROF_START(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING);
      filtering_candidates_process_candidates(filtering_candidates,search->archive,pattern,actual_parameters,true,mm_stack);
      filtering_candidates_verify_candidates(filtering_candidates,search->archive,
          search->text_collection,pattern,actual_parameters,matches,mm_stack);
      filtering_candidates_align_candidates(filtering_candidates,search->archive->text,
          search->text_collection,pattern,search->emulated_rc_search,actual_parameters,false,matches,mm_stack);
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
GEM_INLINE void approximate_search_generate_candidates(
    approximate_search_t* const search,matches_t* const matches,
    const region_filter_type filter_type,const uint64_t filtering_degree,
    const bool dynamic_scheduling,const bool verify_ahead) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Process the proper region-filter scheme
  const uint64_t filtering_threshold = parameters->filtering_threshold;
  const double proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const uint64_t sensibility_error_length = parameters->filtering_region_factor*proper_length;
  if (filter_type == region_filter_fixed) {
    region_profile_schedule_filtering_fixed(region_profile,
        ALL,filtering_degree,filtering_threshold);
    approximate_search_region_profile_generate_candidates(
        search,false,verify_ahead,sensibility_error_length,filtering_threshold,matches);
  } else if (filter_type == region_filter_adaptive_exact) {
    region_profile_schedule_filtering_fixed(region_profile,
        search->max_differences+1,REGION_FILTER_DEGREE_ZERO,filtering_threshold);
    approximate_search_region_profile_generate_candidates(
        search,false,verify_ahead,sensibility_error_length,filtering_threshold,matches);
  } else { // filter_type == region_filter_adaptive_dynamic
    if (!dynamic_scheduling) {
      region_profile_schedule_filtering_adaptive(
          region_profile,search->max_differences,sensibility_error_length);
    }
    approximate_search_region_profile_generate_candidates(
        search,dynamic_scheduling,verify_ahead,sensibility_error_length,filtering_threshold,matches);
  }
  // Update MCS (maximum complete stratum) [Hint]
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
}
/*
 * Mapping workflows 4.0
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
  filtering_candidates_align_candidates(filtering_candidates,search->archive->text,search->text_collection,
      pattern,search->emulated_rc_search,actual_parameters,false,matches,search->mm_stack);
  // Adjust max-differences
  approximate_search_adjust_max_differences_using_strata(search,matches);
  // Update MCS (maximum complete stratum)
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  search->max_matches_reached = num_matches >= parameters->max_search_matches;
  if (search->max_matches_reached) approximate_search_update_mcs(search,0);
  // Next State
  search->search_state = asearch_candidates_verified;
}
/*
 * Neighborhood Generation (Inexact Search)
 */
GEM_INLINE void approximate_search_neighborhood_exact_search(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_EXACT_SEARCH);
  pattern_t* const pattern = &search->pattern;
  // FM-Index basic exact search
//  fm_index_bsearch(search->archive->fm_index,pattern->key,
//    pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_bsearch_pure(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
//  fm_index_reverse_bsearch_pure(search->archive->fm_index,pattern->key,
//      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
  fm_index_bsearch(search->archive->fm_index,pattern->key,
      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
  // Add to matches
  matches_add_interval_match(matches,search->lo_exact_matches,
      search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search);
  // Update MCS
  approximate_search_update_mcs(search,1);
  // Update next state
  search->search_state = asearch_end;
  PROF_STOP(GP_AS_EXACT_SEARCH);
}
GEM_INLINE void approximate_search_neighborhood_inexact_search(approximate_search_t* const search,matches_t* const matches) {
  // Parameters, pattern & interval-set
  const as_parameters_t* const actual_parameters = search->as_parameters;
  pattern_t* const pattern = &search->pattern;
  interval_set_t* const intervals_result = search->interval_set;
  // Basic search (Brute force mitigated by mrank_table)
  interval_set_clear(search->interval_set); // Clear
  neighborhood_search(search->archive->fm_index,pattern->key,pattern->key_length,
      actual_parameters->max_search_error_nominal,intervals_result,search->mm_stack);
  // Add results
  matches_add_interval_set(matches,intervals_result,pattern->key_length,search->emulated_rc_search);
  // Update MCS
  approximate_search_update_mcs(search,actual_parameters->max_search_error_nominal + 1);
  // Update next state
  search->search_state = asearch_end;
}
GEM_INLINE void approximate_search_neighborhood_search(approximate_search_t* const search,matches_t* const matches) {
  // Check max-differences allowed
  if (search->max_differences==0) {
    approximate_search_neighborhood_exact_search(search,matches);   // Exact Search
  } else {
    approximate_search_neighborhood_inexact_search(search,matches); // Basic brute force search
  }
  search->search_state = asearch_end; // Update search state
}
/*
 * Read Recovery
 */
GEM_INLINE void approximate_search_read_recovery(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  PROF_START(GP_AS_READ_RECOVERY);
  /*
   * If the number of wildcards (or errors required) is greater than the maximum number of
   * differences allowed, we try to recover as many matches as possible.
   * We extract feasible regions from the read and filter them trying to recover anything
   * out of bad quality reads.
   */
  // Parameters
  as_parameters_t* const actual_parameters = search->as_parameters;
  search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  // Select proper key
  const archive_t* archive = search->archive;
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
  region_profile_generate_adaptive(region_profile,search->archive->fm_index,key,key_length,
      parameters->allowed_enc,&parameters->rp_recovery,search->max_differences+1,false);
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) { region_profile_print(stderr,region_profile,false,false); }
  // Append the results to filter
  approximate_search_add_regions_to_filter(filtering_candidates,region_profile,search->mm_stack);
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
  // Process candidates
  filtering_candidates_process_candidates(filtering_candidates,
      search->archive,pattern,actual_parameters,verify_candidates,search->mm_stack);
  // Verify/Process candidates
  if (verify_candidates) {
    approximate_search_verify_candidates(search,matches);
  } else {
    // Process candidates (just prepare)
    filtering_candidates_process_candidates(filtering_candidates,search->archive,pattern,actual_parameters,false,search->mm_stack);
    search->search_state = asearch_verify_candidates;
  }
  PROF_STOP(GP_AS_READ_RECOVERY);
}
GEM_INLINE void approximate_search_exact_filtering_adaptive(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT);
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  // Region-Minimal Profile (Reduce the number of candidates per region and maximize number of regions)
  approximate_search_generate_region_profile(search,
      &parameters->rp_minimal,region_profile_adaptive_lightweight,0,false);
  // Check corner cases
  if (search->search_state==asearch_no_regions || search->search_state==asearch_exact_matches) {
    asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
  } else {
    // Region-Minimal Profile Stats
    approximate_search_generate_region_profile_minimal_stats(search);
    // Generate exact-candidates
    approximate_search_generate_candidates(search,matches,region_filter_adaptive_exact,0,false,verify_candidates);
    // Process candidates (just prepare to verification)
    filtering_candidates_process_candidates(search->filtering_candidates,
        search->archive,pattern,actual_parameters,verify_candidates,search->mm_stack);
    // Verify Candidates (if needed)
    if (verify_candidates) {
      approximate_search_verify_candidates(search,matches);
    } else {
      search->search_state = asearch_verify_candidates;
    }
  }
  PROF_STOP(GP_AS_FILTERING_EXACT);
}
GEM_INLINE void approximate_search_exact_filtering_boost(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT_BOOST);
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  // Region-Boost Profile
//  fm_index_t* const fm_index = search->archive->fm_index;
//  const uint64_t base_num_regions = search->region_profile.num_filtering_regions;
//  const uint64_t potential_num_regions = pattern->key_length/(fm_index->proper_length);
//  if (base_num_regions+1 >= potential_num_regions) {
//    asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
//    PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
//  }


  approximate_search_generate_region_profile(search,
      &parameters->rp_boost,region_profile_adaptive_extensive,0,false);
//  if (//base_num_regions <= search->region_profile.num_filtering_regions ||
//      search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
//    asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
//    PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
//    return;
//  }

  approximate_search_generate_region_profile_boost_stats(search); // Stats
  // Generate exact-candidates (Dynamic filtering incorporated)
  approximate_search_generate_candidates(search,matches,region_filter_adaptive_exact,0,true,true);
  // Process candidates (just prepare to verification)
  filtering_candidates_process_candidates(search->filtering_candidates,
      search->archive,pattern,actual_parameters,true,search->mm_stack);
  // Verify candidates
  approximate_search_verify_candidates(search,matches);
  search->search_state = asearch_exact_filtering_boost; // Correct State


  // Next State
  asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
  PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_BOOST_MAPPED,matches_is_mapped(matches)?1:0);
  PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_BOOST_MCS,search->max_complete_stratum);
  PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
}
GEM_INLINE void approximate_search_inexact_filtering(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_INEXACT);
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  // Region-Delimit Profile (Maximize number of unique-regions and try to isolate standard/repetitive regions)
  approximate_search_generate_region_profile(search,
      &parameters->rp_delimit,region_profile_adaptive_lightweight,0,true);
  if (search->region_profile.num_filtering_regions <= 1) {
    asearch_control_next_state(search,asearch_inexact_filtering,matches);
    PROF_STOP(GP_AS_FILTERING_INEXACT);
    return;
  }
  if (search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
    PROF_STOP(GP_AS_FILTERING_INEXACT);
    return;
  }
  approximate_search_generate_region_profile_delimit_stats(search); // Stats
  // Generate exact-candidates (Dynamic filtering incorporated)
  approximate_search_generate_candidates(search,matches,region_filter_adaptive_dynamic,0,true,true);
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
}
GEM_INLINE void approximate_search_local_alignments_filtering(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_LOCAL_ALIGNMENTS);
  // Local-Align discarded candidates
  filtering_candidates_local_alignment(
      search->filtering_candidates,search->archive,search->text_collection,
      &search->pattern,search->emulated_rc_search,search->as_parameters,
      matches,search->mm_stack);
  // Next State
  asearch_control_next_state(search,asearch_local_alignment,matches);
  PROF_ADD_COUNTER(GP_AS_FILTERING_LOCAL_ALIGNMENTS_MAPPED,matches_is_mapped(matches)?1:0);
  PROF_ADD_COUNTER(GP_AS_FILTERING_LOCAL_ALIGNMENTS_MCS,search->max_complete_stratum);
  PROF_STOP(GP_AS_FILTERING_LOCAL_ALIGNMENTS);
}
/*
 * [GEM-workflow 4.0] Adaptive mapping
 *
 *   Filtering-only approach indented to adjust the degree of filtering w.r.t
 *   the structure of the read. Thus, in general terms, a read with many regions
 *   will enable this approach to align the read up to more mismatches than a read
 *   with less number of regions.
 *   Fast-mapping (in all its kinds) tries to detect the proper degree of filtering
 *   to achieve a compromise between speed and depth of the search (max_mismatches)
 */
GEM_INLINE void approximate_search_adaptive_mapping_basic_cases(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  // Parameters
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // Check if all characters are wildcards
  if (key_length==num_wildcards) {
    search->hi_exact_matches = 0;
    search->lo_exact_matches = 0;
    approximate_search_update_mcs(search,key_length);
    search->search_state = asearch_end;
    return;
  }
  // Check if recovery is needed
  if (num_wildcards > 0 && num_wildcards >= search->max_differences) {
    search->search_state = asearch_read_recovery;
    return;
  }
  // Exact search
  if (search->max_differences==0) {
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
  return;
}
GEM_INLINE void approximate_search_adaptive_mapping(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  pattern_t* const pattern = &search->pattern;
  const bool verify_candidates = search->verify_candidates;
  // Process proper search-stage
  while (search->search_state != asearch_end) {
    switch (search->search_state) {
      case asearch_begin:
        // Search Start. Check basic cases
        approximate_search_adaptive_mapping_basic_cases(search,verify_candidates,matches);
        if (search->search_state==asearch_verify_candidates && !verify_candidates) return; // Return if no-filtering
        break;
      case asearch_read_recovery:
        // Read recovery
        approximate_search_read_recovery(search,verify_candidates,matches);
        if (!verify_candidates) return; // Return if no-filtering
        break;
      case asearch_exact_filtering_adaptive:
        // Exact-Filtering (Adaptive)
        approximate_search_exact_filtering_adaptive(search,verify_candidates,matches);
        if (!verify_candidates) return; // Return if no-filtering
        break;
      case asearch_verify_candidates:
        // Verify Candidates
        approximate_search_verify_candidates(search,matches);
        break;
      case asearch_candidates_verified:
        // Next State
        PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
        PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
        asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
        break;
      case asearch_exact_filtering_boost:
        // Exact-Filtering (Fixed)
        approximate_search_exact_filtering_boost(search,matches);
        break;
      case asearch_inexact_filtering:
        // Inexact-Filtering
        approximate_search_inexact_filtering(search,matches);
        break;
      case asearch_neighborhood:
        approximate_search_neighborhood_search(search,matches);
        break;
      case asearch_exact_matches: // Exact Matches
        // Add interval
        matches_add_interval_match(matches,search->lo_exact_matches,
            search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search);
        search->max_differences = actual_parameters->complete_strata_after_best_nominal; // Adjust max-differences
        approximate_search_update_mcs(search,1);
        search->search_state = asearch_end;
        break;
      case asearch_local_alignment: // Local alignments
        approximate_search_local_alignments_filtering(search,matches);
        break;
      case asearch_no_regions: // No regions found
        approximate_search_update_mcs(search,pattern->num_wildcards);
        search->search_state = asearch_end;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Done!!
  PROF_ADD_COUNTER(GP_AS_ADAPTIVE_MCS,search->max_complete_stratum);
}
/*
 * Complete Search
 */
GEM_INLINE void approximate_search_filtering_complete_search(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT);
  if (search->max_differences==0) {
    approximate_search_neighborhood_exact_search(search,matches); // Exact Search
  } else {
    // Parameters
    const as_parameters_t* const actual_parameters = search->as_parameters;
    search_parameters_t* const parameters = actual_parameters->search_parameters;
    fm_index_t* const fm_index = search->archive->fm_index;
    pattern_t* const pattern = &search->pattern;
    region_profile_t* const region_profile = &search->region_profile;
    // Select proper key
    const archive_t* archive = search->archive;
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
    region_profile_generate_adaptive_limited(region_profile,fm_index,key,key_length,
        parameters->allowed_enc,&parameters->rp_minimal,search->max_differences+1);
    if (region_profile_has_exact_matches(region_profile)) {
      const region_search_t* const first_region = region_profile->filtering_region;
      matches_add_interval_match(matches,first_region->lo,
          first_region->hi,pattern->key_length,0,search->emulated_rc_search);
      approximate_search_update_mcs(search,1); // Update MCS
    } else {
      approximate_search_generate_region_profile_minimal_stats(search); // Stats
      // Generate exact-candidates
      parameters->filtering_threshold = ALL; // No restriction
      approximate_search_generate_candidates(search,
          matches,region_filter_fixed,REGION_FILTER_DEGREE_ZERO,false,false);
      // Process candidates (just prepare to verification)
      filtering_candidates_process_candidates(search->filtering_candidates,
          search->archive,pattern,actual_parameters,true,search->mm_stack);
      // Verify candidates
      filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
      filtering_candidates_verify_candidates(
          filtering_candidates,search->archive,search->text_collection,
          pattern,actual_parameters,matches,search->mm_stack);
      // Align candidates
      filtering_candidates_align_candidates(filtering_candidates,search->archive->text,search->text_collection,
          pattern,search->emulated_rc_search,actual_parameters,false,matches,search->mm_stack);
      approximate_search_update_mcs(search,search->region_profile.errors_allowed + pattern->num_wildcards); // Update MCS
    }
  }
  PROF_STOP(GP_AS_FILTERING_EXACT);
}
/*
 * Approximate String Matching using the FM-index
 */
GEM_INLINE void approximate_search(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_MAIN);
  /*
   * Select mapping strategy
   */
  const search_parameters_t* const parameters = search->as_parameters->search_parameters;
  switch (parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
    case mapping_adaptive_filtering_match:
    case mapping_adaptive_filtering_complete:
      approximate_search_adaptive_mapping(search,matches); // Adaptive & incremental mapping
      break;
    case mapping_neighborhood_search:
      approximate_search_neighborhood_search(search,matches); // Brute-force mapping
      break;
    case mapping_filtering_complete:
      approximate_search_filtering_complete_search(search,matches); // Filtering Complete
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Set matches-MCS
  if (matches) matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->max_complete_stratum);
  PROF_STOP(GP_AS_MAIN);
}
/*
 * Approximate String Matching using the FM-index (Verification)
 */
GEM_INLINE void approximate_search_verify(approximate_search_t* const search,matches_t* const matches) {
  // Verify
  const uint64_t num_accepted_regions = filtering_candidates_verify_candidates(
      search->filtering_candidates,search->archive,search->text_collection,
      &search->pattern,search->as_parameters,matches,search->mm_stack);
  if (num_accepted_regions > 0) {
    // Realign
    filtering_candidates_align_candidates(search->filtering_candidates,search->archive->text,
        search->text_collection,&search->pattern,search->emulated_rc_search,
        search->as_parameters,true,matches,search->mm_stack);
  }
  // Update state
  if (search->search_state==asearch_verify_candidates) {
    search->search_state = asearch_candidates_verified;
  }
}
GEM_INLINE void approximate_search_verify_using_bpm_buffer(
    approximate_search_t* const search,
    matches_t* const matches,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end) {
  // Retrieve
  const uint64_t num_accepted_regions = filtering_candidates_bpm_buffer_retrieve(
      search->filtering_candidates,search->archive->text,search->text_collection,
      &search->pattern,bpm_gpu_buffer,candidate_offset_begin,candidate_offset_end,search->mm_stack);
  if (num_accepted_regions > 0) {
    // Realign
    PROF_START(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS);
    filtering_candidates_align_candidates(search->filtering_candidates,search->archive->text,
        search->text_collection,&search->pattern,search->emulated_rc_search,
        search->as_parameters,true,matches,search->mm_stack);
    PROF_STOP(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS);
  }
  // Update state
  if (search->search_state==asearch_verify_candidates) {
    search->search_state = asearch_candidates_verified;
  }
}
GEM_INLINE void approximate_search_hold_verification_candidates(approximate_search_t* const search) {
  filtering_candidates_set_all_regions_pending(search->filtering_candidates);
}
GEM_INLINE void approximate_search_release_verification_candidates(approximate_search_t* const search) {
  filtering_candidates_set_all_regions_unverified(search->filtering_candidates);
  if (search->search_state==asearch_candidates_verified) {
    search->search_state = asearch_verify_candidates;
  }
}
