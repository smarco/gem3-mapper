/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 * TODO
 *   Implement some type of re-seeding
 */

#include "approximate_search.h"

#include "quality_model.h"
#include "rank_mtable.h"
#include "neighborhood_search.h"

#include "mapper_profile.h"

/*
 * Debug
 */
#define DEBUG_REGION_PROFILE_PRINT  false
#define DEBUG_REGION_SCHEDULE_PRINT false

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
    search_actual_parameters_t* const search_actual_parameters) {
  // Index Structures & Parameters
  search->archive = archive;
  search->search_actual_parameters = search_actual_parameters;
  filtering_candidates_init(&search->filtering_candidates);
}
GEM_INLINE void approximate_search_configure(
    approximate_search_t* const search,text_collection_t* text_collection,
    interval_set_t* const interval_set,mm_stack_t* const mm_stack) {
  // Set Auxiliary Structures (external)
  search->text_collection = text_collection;  // Text-Collection
  search->interval_set = interval_set;        // Interval Set
  search->mm_stack = mm_stack;                // Set MM
}
GEM_INLINE void approximate_search_reset(approximate_search_t* const search) {
  // Reset Approximate Search State
  search->search_state = asearch_begin;
  search->verify_candidates = true;
  search->stop_before_neighborhood_search = false;
  search->max_differences = search->search_actual_parameters->max_search_error_nominal;
  search->max_complete_stratum = ALL;
  search->max_matches_reached = false;
  filtering_candidates_clear(&search->filtering_candidates);
}
GEM_INLINE void approximate_search_destroy(approximate_search_t* const search) {
  filtering_candidates_destroy(&search->filtering_candidates);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t approximate_search_get_num_potential_candidates(const approximate_search_t* const search) {
  const filtering_candidates_t* const filtering_candidates = &search->filtering_candidates;
  return filtering_candidates_get_num_candidate_regions(filtering_candidates);
}
/*
 * Pattern
 */
GEM_INLINE void approximate_search_pattern_prepare(approximate_search_t* const search,sequence_t* const sequence) {
  // MM
  mm_stack_t* const mm_stack = search->mm_stack;
  // Parameters
  const search_actual_parameters_t* const actual_parameters = search->search_actual_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  // Set quality search
  search->do_quality_search = (parameters->quality_format!=qualities_ignore) && sequence_has_qualities(sequence);
  // Allocate pattern memory
  const uint64_t read_length = sequence_get_length(sequence);
  pattern_t* const pattern = &search->pattern;
  pattern->key_length = read_length;
  pattern->key = mm_stack_calloc(mm_stack,read_length,uint8_t,false);
  // Build quality model & mask
  if (search->do_quality_search) {
    pattern->quality_mask =  mm_stack_calloc(mm_stack,read_length,uint8_t,false);
    quality_model(sequence,parameters->quality_model,
        parameters->quality_format,parameters->quality_threshold,pattern->quality_mask);
  } else {
    pattern->quality_mask = NULL;
  }
  /*
   * Check all characters in the key & encode key
   * Counts the number of wildcards(characters not allowed as replacements) & low_quality_bases
   */
  uint64_t i, num_wildcards=0, num_low_quality_bases=0;
  const char* const read = sequence_get_read(sequence);
  if (pattern->quality_mask == NULL) {
    for (i=0;i<read_length;++i) {
      const char character = read[i];
      if (!parameters->allowed_chars[(uint8_t)character]) ++num_wildcards;
      pattern->key[i] = dna_encode(character);
    }
  } else {
    for (i=0;i<read_length;++i) {
      const char character = read[i];
      if (!parameters->allowed_chars[(uint8_t)character]) {
        ++num_low_quality_bases; ++num_wildcards;
      } else if (pattern->quality_mask[i]!=qm_real) {
        ++num_low_quality_bases;
      }
      pattern->key[i] = dna_encode(character);
    }
  }
  pattern->num_wildcards = num_wildcards;
  pattern->num_low_quality_bases = num_low_quality_bases;
  // Calculate the effective number of differences
  const int64_t max_allowed_error = (int64_t)read_length - (int64_t) actual_parameters->min_matching_length_nominal;
  uint64_t effective_filtering_max_error;
  if (gem_expect_false(max_allowed_error<=0)) { // Constrained by min_matching_length_nominal
    effective_filtering_max_error = 0;
  } else {
    // Constrained by num_low_quality_bases
    effective_filtering_max_error = actual_parameters->max_filtering_error_nominal + pattern->num_low_quality_bases;
    if (effective_filtering_max_error > max_allowed_error) {
      effective_filtering_max_error = max_allowed_error;
    }
  }
  pattern->max_effective_filtering_error = effective_filtering_max_error;
  pattern->max_effective_bandwidth = actual_parameters->max_bandwidth_nominal;
  // Prepare region profile
  region_profile_new(&search->region_profile,read_length,mm_stack);
  // Prepare kmer-counting filter
//  kmer_counting_compile(&pattern->kmer_counting,(bool* const)parameters->allowed_enc,pattern->key,read_length,mm_stack);
  // Prepare BPM pattern
  bpm_pattern_compile(&pattern->bpm_pattern,pattern->key,read_length,effective_filtering_max_error,mm_stack);
  // Prepare SWG query-profile
//  if (parameters->alignment_model == alignment_model_gap_affine) {
//    swg_init_query_profile(&pattern->swg_query_profile,&parameters->swg_penalties,read_length,mm_stack);
//  }
}
GEM_INLINE void approximate_search_pattern_clear(approximate_search_t* const search) {
  search->pattern.key_length = 0; // Clear the pattern
}
GEM_INLINE bool approximate_search_pattern_is_null(approximate_search_t* const search) {
  return (search->pattern.key_length == 0);
}
/*
 * Debug
 */
GEM_INLINE void approximate_search_region_profile_schedule_filtering_degree_print(
    region_profile_t* const region_profile,const uint64_t max_differences,
    const uint64_t sensibility_error_length) {
  // Header
  gem_slog("[GEM]>Region.Filtering.Schedule\n");
  gem_slog("  => Max.Differences %lu\n",max_differences);
  gem_slog("  => Min.Sensibility.Region.Length %lu\n",sensibility_error_length);
  gem_slog("  => Regions\n");
  REGION_LOCATOR_ITERATE(region_profile,region,position) {
    gem_slog("    [%lu] Region-%s\t(%3lu,%3lu]\t\t(Length,Cand)=(%3lu,%4lu)",
        position,region->type==region_unique ? "unique" : "standard",
        region->end,region->start,region->start-region->end,region->hi-region->lo);
    if (region->min==0) {
      gem_slog("\tDegree=none\n");
    } else {
      gem_slog("\tDegree=%lu\n",region->min-1);
    }
  }
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
  const uint64_t delta = search->search_actual_parameters->complete_strata_after_best_nominal;
  /*
   * If delta parameter is set (and is below the maximum number of mismatches),
   * finds the minimum non zero stratum (mnzs) and adjusts
   * the maximum number of mismatches to (mnzs+delta)
   */
  if (delta < max_differences) {
    const int64_t fms = matches_counters_get_min(matches);
    if (fms>=0 && fms+delta < max_differences) {
      search->max_differences = fms+delta;
    }
  }
}
GEM_INLINE bool approximate_search_adaptive_mapping_fulfilled(
    approximate_search_t* const search,matches_t* const matches) {
  const search_parameters_t* const parameters = search->search_actual_parameters->search_parameters;
  // Determines when the search is done following the mapping criteria
  switch (parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast: // Anything is good
      return true;
    case mapping_adaptive_filtering_match: { // We need a match (or any evindence that we will get a match)
      if (search->search_state==asearch_no_regions) return false;
      if (search->search_state==asearch_exact_matches) return true;
      if (search->search_state==asearch_inexact_filtering) return true; // Stop before NS
      if (matches_is_mapped(matches)) return true;
      if (search->max_complete_stratum >= search->max_differences) {
        // Not mapped but meeting the search requirements // TODO Tricky case => Maybe straight to NS
        return true;
      } else {
        return false;
      }
    }
    case mapping_adaptive_filtering_complete:
      return search->max_complete_stratum >= search->max_differences;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}
/*
 * Region Profile
 */
GEM_INLINE void approximate_search_generate_region_profile(
    approximate_search_t* const search,const region_profile_type profile_type,
    const region_profile_model_t* const profile_model) {
  // Parameters
  const search_actual_parameters_t* const actual_parameters = search->search_actual_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  fm_index_t* const fm_index = search->archive->fm_index;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Compute the region profile
  switch (profile_type) {
    case region_profile_fixed:
      GEM_NOT_IMPLEMENTED(); // TODO
      break;
    case region_profile_adaptive:
      region_profile_generate_adaptive(region_profile,fm_index,pattern,
          parameters->allowed_enc,profile_model,search->max_differences+1,false);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) { region_profile_print(stderr,region_profile,false,false); }
  // Check Zero-Region & Exact-Matches
  if (region_profile->num_filtering_regions==0) {
    search->max_complete_stratum = pattern->num_wildcards;
    search->search_state = asearch_no_regions;
    return;
  } else if (region_profile_has_exact_matches(region_profile)) {
    const region_search_t* const first_region = region_profile->filtering_region;
    search->hi_exact_matches = first_region->hi;
    search->lo_exact_matches = first_region->lo;
    search->max_complete_stratum = 1;
    search->search_state = asearch_exact_matches;
    return;
  }
}
/*
 * Region Scheduling
 */
GEM_INLINE void approximate_search_region_profile_schedule_filtering_degree_fixed(
    region_profile_t* const region_profile,const uint64_t filtering_degree) {
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    region_locator_t* const loc = region_profile->loc+position;
    loc->id = position;
    loc->value = 0;
    region->min = filtering_degree;
  }
}
GEM_INLINE void approximate_search_region_schedule_filtering_degree(
    region_search_t* const region,const uint64_t num_standard_regions_left,
    const uint64_t num_unique_regions_left,const uint64_t max_differences,
    const uint64_t sensibility_error_length,const uint64_t errors_allowed) {
  if (errors_allowed > max_differences) {
    region->min = REGION_FILTER_NONE; // Search is fulfilled. Don't filter
    return;
  }
  if (region->type == region_standard) {
    region->min = REGION_FILTER_DEGREE_ZERO; // Cannot be filtered with errors
    return;
  }
  // Compute the scope of the search using zero-filter
  const uint64_t total_errors_zero_filter = errors_allowed + num_standard_regions_left + num_unique_regions_left + 1;
  if (total_errors_zero_filter >= max_differences+1) {
    region->min = REGION_FILTER_DEGREE_ZERO; // Search is fulfilled just by filtering-zero
  } else {
    // Compute number of errors left to be assigned to unique-regions
    const uint64_t region_length = region->start-region->end;
    const uint64_t pending_errors_at_unique_regions = (max_differences+1) - errors_allowed - num_standard_regions_left;
    if (num_unique_regions_left >= pending_errors_at_unique_regions  // Enough Regions to filter-zero
        || region_length < sensibility_error_length                  // Region is too small for other degree
        || max_differences==1) {                                     // Search doesn't require to allow more errors
      region->min=REGION_FILTER_DEGREE_ZERO;
    } else if (2*num_unique_regions_left >= pending_errors_at_unique_regions // Enough Regions to filter-one
        || region_length < 2*sensibility_error_length                // Region is too small for other degree
        || max_differences==2) {                                     // Search doesn't require to allow more errors
      region->min=REGION_FILTER_DEGREE_ONE;
    } else {                                                         // Maximum degree reached
      region->min=REGION_FILTER_DEGREE_ONE;
//      region->min=REGION_FILTER_DEGREE_TWO; // FIXME
    }
//    // Mandatory error condition // TODO
//    if (region->hi-region->lo==0 &&                        // Zero candidates regions (Need at least 1 error to match)
//        pending_errors_at_unique_regions > region->min &&  // We have errors pending to be assigned
//        region->min < REGION_FILTER_DEGREE_TWO) {          // We don't go beyond 2 errors threshold
//      ++(region->min);
//    }
  }
}
GEM_INLINE void approximate_search_region_profile_schedule_filtering_degree(
    region_profile_t* const region_profile,const uint64_t max_differences,
    const uint64_t sensibility_misms_length) {
  /*
   * PRE: (region_profile->num_filtering_regions <= max_mismatches)
   * Tries to assign the best possible filtering degree distribution among
   * all the regions to fulfill the requirements of a search up to @max_mismatches
   */
  const uint64_t num_regions = region_profile->num_filtering_regions;
  const uint64_t num_standard_regions = region_profile->num_standard_regions;
  const uint64_t num_unique_regions = num_regions - num_standard_regions;
  // Sort regions
  region_profile_sort_by_estimated_mappability(region_profile);
  // Try to schedule a distribution of the errors over the regions
  uint64_t num_standard_regions_left = num_standard_regions;
  uint64_t num_unique_regions_left = num_unique_regions;
  uint64_t errors_allowed = 0;
  REGION_LOCATOR_ITERATE(region_profile,region,position) {
    if (region->type == region_standard) {
      --num_standard_regions_left;
    } else {
      --num_unique_regions_left;
    }
    approximate_search_region_schedule_filtering_degree(
        region,num_standard_regions_left,num_unique_regions_left,
        max_differences,sensibility_misms_length,errors_allowed);
    errors_allowed += region->min;
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
    approximate_search_t* const search,
    const bool dynamic_scheduling,const uint64_t sensibility_error_length,
    const uint64_t filtering_threshold,matches_t* const matches) {
  PROF_START(GP_AS_GENERATE_CANDIDATES);
  /*
   * Filters all the regions up to the scheduled degree
   *  - Dynamic Scheduling: Assigns a filtering degree to each region as the search goes on
   *      (as opposed to a static scheduling giving in advanced)
   *  - Dynamic Filtering: Filters candidates per each queried region (thus, it can reduce
   *      the scope of the search by delta)
   *  - Locator guided: The region filtering is conducted by the locator order of region.
   */
  // Data
  const search_actual_parameters_t* const actual_parameters = search->search_actual_parameters;
  mm_stack_t* const mm_stack = search->mm_stack;
  // Approximate Search Data
  fm_index_t* const fm_index = search->archive->fm_index;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = &search->filtering_candidates;
  interval_set_t* const intervals_result = search->interval_set;
  // Pattern
  const pattern_t* const pattern = &search->pattern;
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
      approximate_search_region_profile_schedule_filtering_degree(region_profile,search->max_differences,sensibility_error_length);
      approximate_search_region_profile_schedule_filtering_degree_print(region_profile,search->max_differences,sensibility_error_length);
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
      approximate_search_region_schedule_filtering_degree(
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
     * Otherwise, ignore the region (intractable)
     */
    // Check candidates (Dynamic filtering)
    if (dynamic_scheduling && actual_parameters->complete_strata_after_best_nominal < search->max_differences) {
      PROF_PAUSE(GP_AS_GENERATE_CANDIDATES);
      PROF_START(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING);
      filtering_candidates_process_candidates(filtering_candidates,search->archive,pattern,mm_stack);
      filtering_candidates_verify_candidates(filtering_candidates,search->archive,search->text_collection,
          pattern,search->search_strand,actual_parameters,matches,mm_stack);
      approximate_search_adjust_max_differences_using_strata(search,matches);
      PROF_STOP(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING);
      PROF_CONTINUE(GP_AS_GENERATE_CANDIDATES);
    }
    // Check error-condition
    if (errors_allowed > search->max_differences) {
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
    bool filter_dynamic) {
  // Parameters
  const search_actual_parameters_t* const actual_parameters = search->search_actual_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Process the proper region-filter scheme
  const uint64_t filtering_threshold = parameters->filtering_threshold;
  const uint64_t sensibility_error_length =
      parameters->filtering_region_factor*fm_index_get_proper_length(search->archive->fm_index);
  if (filter_type == region_filter_fixed) {
    approximate_search_region_profile_schedule_filtering_degree_fixed(region_profile,filtering_degree);
    approximate_search_region_profile_generate_candidates(
        search,false,sensibility_error_length,filtering_threshold,matches);
  } else {
    if (!filter_dynamic) {
      approximate_search_region_profile_schedule_filtering_degree(
          region_profile,search->max_differences,sensibility_error_length);
    }
    approximate_search_region_profile_generate_candidates(
        search,filter_dynamic,sensibility_error_length,filtering_threshold,matches);
  }
  // Update MCS (maximum complete stratum) [Hint]
  search->max_complete_stratum = region_profile->errors_allowed + pattern->num_wildcards;
}
/*
 * Mapping workflows 4.0
 */
/*
 * Exact search
 */
GEM_INLINE void approximate_search_exact_search(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_EXACT_SEARCH);
  pattern_t* const pattern = &search->pattern;
  // FM-Index basic exact search
  fm_index_bsearch(search->archive->fm_index,pattern->key,
      pattern->key_length,&search->hi_exact_matches,&search->lo_exact_matches);
  // Add to matches
  matches_add_interval_match(matches,search->lo_exact_matches,
      search->hi_exact_matches,pattern->key_length,0,search->search_strand);
  // Update MCS
  search->max_complete_stratum = 1;
  search->search_state = asearch_exact_matches;
  PROF_STOP(GP_AS_EXACT_SEARCH);
}
/*
 * Inexact Search (neighbohood generation)
 */
GEM_INLINE void approximate_search_inexact_search(approximate_search_t* const search,matches_t* const matches) {
  // Parameters, pattern & interval-set
  const search_actual_parameters_t* const actual_parameters = search->search_actual_parameters;
  pattern_t* const pattern = &search->pattern;
  interval_set_t* const intervals_result = search->interval_set;
  // Basic search (Brute force mitigated by mrank_table)
  interval_set_clear(search->interval_set); // Clear
  neighborhood_search(search->archive->fm_index,pattern->key,pattern->key_length,
      actual_parameters->max_search_error_nominal,intervals_result,search->mm_stack);
  // Add results
  matches_add_interval_set(matches,intervals_result,pattern->key_length,search->search_strand);
  // Update MCS
  search->max_complete_stratum = actual_parameters->max_search_error_nominal + 1;
}
/*
 * Read Recovery
 */
GEM_INLINE void approximate_search_read_recovery(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  /*
   * If the number of wildcards (or errors required) is greater than the maximum number of
   * differences allowed, we try to recover as many matches as possible.
   * We extract feasible regions from the read and filter them trying to recover anything
   * out of bad quality reads.
   */
  // Parameters
  search_actual_parameters_t* const actual_parameters = search->search_actual_parameters;
  search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Compute the region profile
  region_profile_generate_adaptive(region_profile,search->archive->fm_index,pattern,
      parameters->allowed_enc,&parameters->rp_recovery,search->max_differences+1,false);
  // Append the results to filter
  filtering_candidates_t* const filtering_candidates = &search->filtering_candidates;
  approximate_search_add_regions_to_filter(filtering_candidates,region_profile,search->mm_stack);
  // Verify/Process candidates
  if (verify_candidates) {
    // Verify candidates
    filtering_candidates_process_candidates(filtering_candidates,search->archive,pattern,search->mm_stack);
    filtering_candidates_verify_candidates(filtering_candidates,search->archive,search->text_collection,
        pattern,search->search_strand,actual_parameters,matches,search->mm_stack);
    // Update MCS (maximum complete stratum)
    search->max_matches_reached = filtering_candidates->total_candidates_accepted >= parameters->max_search_matches;
    search->max_complete_stratum = (search->max_matches_reached) ? 0 : region_profile->errors_allowed + pattern->num_wildcards;
    search->search_state = asearch_candidates_verified;
  } else {
    // Process candidates (just prepare)
    filtering_candidates_process_candidates(filtering_candidates,search->archive,pattern,search->mm_stack);
    // Update MCS (maximum complete stratum)
    search->max_complete_stratum = region_profile->errors_allowed + pattern->num_wildcards;
    search->search_state = asearch_verify_candidates;
  }
}
/*
 * Brute-Force mapping
 */
GEM_INLINE void approximate_search_neighborhood_search(
    approximate_search_t* const search,matches_t* const matches) {
  if (search->max_differences==0) {
    approximate_search_exact_search(search,matches); // Exact Search
  } else {
    approximate_search_inexact_search(search,matches); // Basic brute force search
  }
  search->search_state = asearch_end; // Update search state
  return; // End
}
/*
 * [GEM-workflow 4.0] Basic cases
 */
GEM_INLINE bool approximate_search_adaptive_mapping_basic_cases(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  PROF_START(GP_AS_BASIC);
  // Parameters
  const search_actual_parameters_t* const actual_parameters = search->search_actual_parameters;
  pattern_t* const pattern = &search->pattern;
  // Check if all characters are wildcards
  if (search->pattern.key_length==search->pattern.num_wildcards) {
    search->max_complete_stratum = search->pattern.key_length;
    search->hi_exact_matches = 0;
    search->lo_exact_matches = 0;
    search->search_state = asearch_end;
    PROF_STOP(GP_AS_BASIC);
    return true;
  }
  // Check if recovery is needed
  if (pattern->num_wildcards >= actual_parameters->max_search_error_nominal) {
    PROF_START(GP_AS_READ_RECOVERY);
    approximate_search_read_recovery(search,verify_candidates,matches);
    PROF_STOP(GP_AS_READ_RECOVERY);
    PROF_STOP(GP_AS_BASIC);
    return true;
  }
  // Exact search
  if (actual_parameters->max_search_error_nominal==0) {
    approximate_search_exact_search(search,matches);
    search->search_state = asearch_end;
    PROF_STOP(GP_AS_BASIC);
    return true;
  }
  // Very short reads (Neighborhood search)
  if (pattern->key_length <= RANK_MTABLE_SEARCH_DEPTH ||
      pattern->key_length < search->archive->fm_index->proper_length) {
    PROF_START(GP_AS_SMALL_READS);
    approximate_search_inexact_search(search,matches);
    search->search_state = asearch_end;
    PROF_STOP(GP_AS_SMALL_READS);
    PROF_STOP(GP_AS_BASIC);
    return true;
  }
  PROF_STOP(GP_AS_BASIC);
  return false;
}
/*
 * [GEM-workflow 4.0] Adaptive mapping (Formerly known as fast-mapping)
 *
 *   Filtering-only approach indented to adjust the degree of filtering w.r.t
 *   the structure of the read. Thus, in general terms, a read with many regions
 *   will enable this approach to align the read up to more mismatches than a read
 *   with less number of regions.
 *   Fast-mapping (in all its kinds) tries to detect the proper degree of filtering
 *   to achieve a compromise between speed and depth of the search (max_mismatches)
 */
GEM_INLINE void approximate_search_adaptive_mapping(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const search_actual_parameters_t* const actual_parameters = search->search_actual_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  const bool verify_candidates = search->verify_candidates;
  // Callback (switch to proper search stage)
  approximate_search_adaptive_mapping_callback:
  switch (search->search_state) {
    case asearch_begin:
      // Search Start. Check basic cases
      if (approximate_search_adaptive_mapping_basic_cases(search,verify_candidates,matches)) break;
      // No Break
    case asearch_exact_filtering:
      PROF_START(GP_AS_FILTERING_EXACT);
      // Region-Minimal Profile (Reduce the number of candidates per region and maximize number of regions)
      approximate_search_generate_region_profile(search,region_profile_adaptive,&parameters->rp_minimal);
      // Check corner cases
      if (search->search_state==asearch_no_regions || search->search_state==asearch_exact_matches) {
        PROF_STOP(GP_AS_FILTERING_EXACT);
        if (!verify_candidates) {
          if (!approximate_search_adaptive_mapping_fulfilled(search,matches)) {
            search->search_state = asearch_inexact_filtering;
          }
          return; // Return if no-filtering
        }
        if (!approximate_search_adaptive_mapping_fulfilled(search,matches)) {
          search->search_state = asearch_inexact_filtering;
        }
        goto approximate_search_adaptive_mapping_callback;
      }
      approximate_search_generate_region_profile_minimal_stats(search); // Stats
      // Generate exact-candidates
      approximate_search_generate_candidates(
          search,matches,region_filter_fixed,REGION_FILTER_DEGREE_ZERO,verify_candidates);
      // Process candidates (just prepare to verification)
      filtering_candidates_process_candidates(&search->filtering_candidates,search->archive,pattern,search->mm_stack);
      // Set state
      if (!verify_candidates) {
        search->search_state = asearch_verify_candidates;
        PROF_STOP(GP_AS_FILTERING_EXACT);
        return; // Return if no-filtering
      }
      // No Break
    case asearch_verify_candidates: {
      filtering_candidates_t* const filtering_candidates = &search->filtering_candidates;
      filtering_candidates_verify_candidates(
          filtering_candidates,search->archive,search->text_collection,
          pattern,search->search_strand,actual_parameters,matches,search->mm_stack);
      // Update MCS (maximum complete stratum)
      search->max_matches_reached = filtering_candidates->total_candidates_accepted >= parameters->max_search_matches;
      search->max_complete_stratum = (search->max_matches_reached) ? 0 :
          search->region_profile.errors_allowed + pattern->num_wildcards;
      PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_STOP(GP_AS_FILTERING_EXACT);
      }
      // No Break
    case asearch_candidates_verified:
      if (approximate_search_adaptive_mapping_fulfilled(search,matches)) {
        search->search_state = asearch_end;
        break;
      }
      search->search_state = asearch_inexact_filtering;
      // No Break
    case asearch_inexact_filtering:
      PROF_START(GP_AS_FILTERING_INEXACT);
      // Region-Delimit Profile (Maximize number of unique-regions and try to isolate standard/repetitive regions)
      approximate_search_generate_region_profile(search,region_profile_adaptive,&parameters->rp_delimit);
      if (search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
        PROF_STOP(GP_AS_FILTERING_INEXACT);
        goto approximate_search_adaptive_mapping_callback; // TODO link with complete search
      }
      approximate_search_generate_region_profile_delimit_stats(search); // Stats
      // Generate exact-candidates (Dynamic filtering incorporated)
      approximate_search_generate_candidates(search,matches,region_filter_adaptive_dynamic,0,true);
      filtering_candidates_process_candidates(&search->filtering_candidates,search->archive,pattern,search->mm_stack);
      filtering_candidates_verify_candidates(&search->filtering_candidates,search->archive,search->text_collection,
          pattern,search->search_strand,actual_parameters,matches,search->mm_stack);
      approximate_search_adjust_max_differences_using_strata(search,matches);
      search->search_state = asearch_end;
      PROF_ADD_COUNTER(GP_AS_FILTERING_INEXACT_MAPPED,matches_is_mapped(matches)?1:0);
      PROF_STOP(GP_AS_FILTERING_INEXACT);
      break;
    /*
     * Corner Cases
     */
    case asearch_exact_matches:
      // Add interval
      matches_add_interval_match(matches,search->lo_exact_matches,
          search->hi_exact_matches,pattern->key_length,0,search->search_strand);
      search->max_differences = actual_parameters->complete_strata_after_best_nominal; // Adjust max-differences
      search->search_state = asearch_end;
      break;
    case asearch_no_regions:
      search->search_state = asearch_end;
      break;
    /*
     * EOW (End-of-workflow)
     */
    case asearch_end:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Done!!
  PROF_ADD_COUNTER(GP_AS_ADAPTIVE_MCS,search->max_complete_stratum);
}
/*
 * Approximate String Matching using the FM-index.
 */
GEM_INLINE void approximate_search(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_MAIN);
  /*
   * Select mapping strategy
   */
  const search_parameters_t* const parameters = search->search_actual_parameters->search_parameters;
  switch (parameters->mapping_mode) {
    case mapping_adaptive_filtering_fast:
    case mapping_adaptive_filtering_match:
    case mapping_adaptive_filtering_complete:
      approximate_search_adaptive_mapping(search,matches); // Adaptive & incremental mapping
      break;
    case mapping_neighborhood_search:
      approximate_search_neighborhood_search(search,matches); // Brute-force mapping
      break;
    case mapping_lab_testing: // Reserved for testing purposes
      GEM_NOT_IMPLEMENTED();
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Set MCS
  if (matches) {
    matches->max_complete_stratum = MIN(matches->max_complete_stratum,search->max_complete_stratum);
  }
  PROF_STOP(GP_AS_MAIN);
}
GEM_INLINE void approximate_search_verify_using_bpm_buffer(
    approximate_search_t* const search,
    matches_t* const matches,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t candidate_offset_begin,const uint64_t candidate_offset_end) {
  if (search->search_state==asearch_verify_candidates) {
    filtering_candidates_bpm_buffer_align(
        &search->filtering_candidates,search->archive,search->text_collection,
        &search->pattern,search->search_strand,search->search_actual_parameters,
        bpm_gpu_buffer,candidate_offset_begin,candidate_offset_end,matches,search->mm_stack);
    search->search_state = asearch_candidates_verified;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//  // Soft Probing the matches to lower max_mismatches using delta
//  if (search_params->delta < search_params->max_mismatches) {
//    PROF_START(SC_PROBING_DELTA);
//    fmi_probe_matches(fmi,search_params,region_profile,false,matches,mpool);
//    if (region_profile->num_filtering_regions > search_params->max_mismatches) {
//      PROF_INC_COUNTER(GSC_DELTA_PROBE_HIT); PROF_STOP(SC_PROBING_DELTA);
//      FMI_ZERO_FILTERING();
//      return true; // Well done !
//    }
//    PROF_STOP(SC_PROBING_DELTA);
//  }

//
///*
// * Filters all the regions in the profile up to zero mismatches
// * From the results tries to reduce the maximum number of mismatches
// * using delta.
// */
//GEM_INLINE void fmi_probe_matches(const _FMI_* const fmi,
//    fmi_search_parameters* const search_params,region_profile* const region_profile,
//    const bool sparse_probing,matches* const matches,vector_pool* const mpool) {
//  // We add all the exact search of the region to the filter query buffer
//  vector* const query_buffer = mpool->fbuf1;
//  vector_init(query_buffer,filter_query);
//
//  // Add candidates to filter
//  if (!sparse_probing) {
//    approximate_search_add_regions_to_filter(region_profile,query_buffer);
//  } else {
//    fmi_append_sampled_region_zero_degree_to_filter(region_profile,
//        search_params->internal_parameters.filtering_threshold,query_buffer);
//  }
//
//  // Filter all the exact regions
//  PROF_ADD_COUNTER(GSC_DELTA_PROBE_CAND,vector_get_used(query_buffer));
//  fmi_matches_filter__append_decoded(fmi,matches,search_params,
//      query_buffer,mpool->ibuf1,mpool->ibuf2,true,false,mpool);
//
////  // Adjust the maximum number of mismatches
////  approximate_search_adjust_max_differences_using_strata(search_params,matches); // TODO; Skip as done internally
//
//  // Clean matches buffers
//  vector_clean(mpool->ibuf1);
//  vector_clean(mpool->ibuf2);
//  vector_clean(mpool->fbuf1);
//}
//
//
//GEM_INLINE void fmi_append_sampled_region_zero_degree_to_filter(
//    region_profile* const region_profile,const uint64_t filtering_threshold,
//    vector* const vector_buffer) {
//  const filtering_region* const region = region_profile->filtering_region;
//  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
//  uint64_t i, total_candidates;
//  for (i=0,total_candidates=0; i<num_filtering_regions; ++i) {
//    total_candidates+=region[i].hi-region[i].lo;
//  }
//  // Append sampled candidates
//  const uint64_t sampling_gap = (total_candidates+filtering_threshold)/filtering_threshold;
//  uint64_t total_used = vector_get_used(vector_buffer);
//  for (i=0; i<num_filtering_regions; ++i) {
//    const uint64_t count = region[i].hi-region[i].lo;
//    if (__builtin_expect(count==0,false)) continue;
//    uint64_t num_sampled = (count+sampling_gap)/sampling_gap;
//    vector_reserve(vector_buffer,total_used+num_sampled);
//    filter_query* queries_to_filter = (filter_query*)vector_get_mem(vector_buffer)+total_used;
//    uint64_t pos;
//    for (num_sampled=0,pos=region[i].lo; pos<region[i].hi; pos+=sampling_gap) {
//      queries_to_filter->start = region[i].start;
//      queries_to_filter->end = region[i].end;
//      queries_to_filter->misms = 0;
//      queries_to_filter->end_reg_idx = pos;
//      ++queries_to_filter; ++num_sampled;
//    }
//    total_used+=num_sampled;
//  }
//  vector_set_used(vector_buffer,total_used);
//}

