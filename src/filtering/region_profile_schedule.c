/*
 * PROJECT: GEMMapper
 * FILE: region_profile_schedule.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/region_profile_schedule.h"

/*
 * Region Profile Scheduling
 */
void region_profile_schedule_filtering_exact(region_profile_t* const region_profile) {
  // Parameters
  const uint64_t num_regions = region_profile->num_filtering_regions;
  region_search_t* const filtering_region = region_profile->filtering_region;
  // Assign degree-zero to all regions
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    filtering_region[i].degree = REGION_FILTER_DEGREE_ZERO;
  }
}
