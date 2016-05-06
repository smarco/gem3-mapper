/*
 * PROJECT: GEMMapper
 * FILE: region_profile_schedule.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/region_profile_schedule.h"

/*
 * Region Scheduling
 */
void region_schedule_filtering_adaptive(
    region_search_t* const region,
    const uint64_t num_standard_regions_left,
    const uint64_t num_unique_regions_left,
    const uint64_t max_complete_error,
    const uint64_t sensibility_error_length,
    const uint64_t errors_allowed) {
  if (errors_allowed > max_complete_error) {
    region->degree = REGION_FILTER_NONE; // Search is fulfilled. Don't filter
    return;
  }
  if (region->type == region_standard) {
    region->degree = REGION_FILTER_DEGREE_ZERO; // Cannot be filtered with errors
    return;
  }
  // Compute the scope of the search using zero-filter
  const uint64_t total_errors_zero_filter = errors_allowed + num_standard_regions_left + num_unique_regions_left + 1;
  if (total_errors_zero_filter >= max_complete_error+1) {
    region->degree = REGION_FILTER_DEGREE_ZERO; // Search is fulfilled just by filtering-zero
  } else {
    // Compute number of errors left to be assigned to unique-regions
    const uint64_t region_length = region->end-region->begin;
    const uint64_t pending_errors_at_unique_regions = (max_complete_error+1) - errors_allowed - num_standard_regions_left;
    if (num_unique_regions_left >= pending_errors_at_unique_regions  // Enough Regions to filter-zero
        || region_length < sensibility_error_length                  // Region is too small for other degree
        || max_complete_error==1) {                                  // Search doesn't require to allow more errors
      region->degree=REGION_FILTER_DEGREE_ZERO;
    } else if (2*num_unique_regions_left >= pending_errors_at_unique_regions // Enough Regions to filter-one
        || region_length < 2*sensibility_error_length                // Region is too small for other degree
        || max_complete_error==2) {                                  // Search doesn't require to allow more errors
      region->degree=REGION_FILTER_DEGREE_ONE;
    } else {                                                         // Maximum degree reached
      region->degree=REGION_FILTER_DEGREE_ONE;
//      region->min=REGION_FILTER_DEGREE_TWO; // FIXME At the moment NS(2) is too expensive
    }
//    // Mandatory error condition // TODO
//    if (region->hi-region->lo==0 &&                        // Zero candidates regions (Need at least 1 error to match)
//        pending_errors_at_unique_regions > region->min &&  // We have errors pending to be assigned
//        region->min < REGION_FILTER_DEGREE_TWO) {          // We don't go beyond 2 errors threshold
//      ++(region->min);
//    }
  }
}
/*
 * Region Profile Scheduling
 */
void region_profile_schedule_filtering_fixed(
    region_profile_t* const region_profile,
    const uint64_t regions_required,
    const uint64_t filtering_degree,
    const uint64_t filtering_threshold) {
  /*
   * Schedules the @regions_required regions with fewer number of candidates
   *   to be filtered up to zero mismatches.
   */
  const uint64_t num_regions = region_profile->num_filtering_regions;
  region_search_t* const filtering_region = region_profile->filtering_region;
  region_locator_t* const loc = region_profile->loc;
  // Sort by number of candidates
  region_profile_sort_by_candidates(region_profile);
  // Check the number of regions in the profile
  uint64_t i;
  if (filtering_degree==REGION_FILTER_DEGREE_ZERO && num_regions > regions_required) {
    /*
     * More regions than required (2 strategies)
     *   (*) - Filter only those required which have the less number of candidates
     *       - Filter all of them, but only check those
     *           candidates which more than 1 region (H-samples)  // TODO
     */
    // Filter up to 0-mismatches those regions with the less number of candidates
    for (i=0;i<regions_required;++i) {
      filtering_region[loc[i].id].degree = REGION_FILTER_DEGREE_ZERO;
    }
    for (;i<num_regions;++i) {
      filtering_region[loc[i].id].degree = REGION_FILTER_NONE;
    }
  } else {
    for (i=0;i<num_regions;++i) {
      if (loc[i].value <= filtering_threshold) {
        filtering_region[i].degree = filtering_degree;
      } else {
        filtering_region[i].degree = REGION_FILTER_NONE;
      }
    }
  }
}
void region_profile_schedule_filtering_adaptive(
    region_profile_t* const region_profile,
    const uint64_t max_complete_error,
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
    region_schedule_filtering_adaptive(
        region,num_standard_regions_left,num_unique_regions_left,
        max_complete_error,sensibility_misms_length,errors_allowed);
    errors_allowed += region->degree;
  }
}
/*
 * Display
 */
void region_profile_schedule_print(
    region_profile_t* const region_profile,
    const uint64_t max_differences,
    const uint64_t sensibility_error_length) {
  // Header
  gem_slog("[GEM]>Region.Filtering.Schedule\n");
  gem_slog("  => Max.Differences %"PRIu64"\n",max_differences);
  gem_slog("  => Min.Sensibility.Region.Length %"PRIu64"\n",sensibility_error_length);
  gem_slog("  => Regions\n");
  REGION_LOCATOR_ITERATE(region_profile,region,position) {
    gem_slog("    [%"PRIu64"] Region-%s\t(%3"PRIu64",%3"PRIu64"]"
    		"\t\t(Length,Cand)=(%3"PRIu64",%4"PRIu64")",
        position,region->type==region_unique ? "unique" : "standard",
        region->begin,region->end,region->end-region->begin,region->hi-region->lo);
    if (region->degree==0) {
      gem_slog("\tDegree=none\n");
    } else {
      gem_slog("\tDegree=%"PRIu64"\n",region->degree-1);
    }
  }
}
