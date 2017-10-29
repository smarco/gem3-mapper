/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Region-Profile module provides functions to generate candidate
 *   bwt-positions from a region-profile.
 */

#include "filtering/region_profile/region_profile_schedule.h"

/*
 * Sorting regions
 */
typedef struct {
  int idx;
  int value;
} region_sort_t;
int region_sort_cmp_value(
    const region_sort_t* const a,
    const region_sort_t* const b) {
  return (int)a->value - (int)b->value;
}
void region_sort_sort(
    region_sort_t* const region_sort,
    const uint64_t num_regions) {
  qsort(region_sort,num_regions,sizeof(region_sort_t),
      (int (*)(const void *,const void *))region_sort_cmp_value);
}
/*
 * Schedule all exact
 */
bool region_profile_schedule_exact_matches(region_profile_t* const region_profile) {
  // Parameters
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  region_search_t* const filtering_region = region_profile->filtering_region;
  // Detect
  if (num_filtering_regions==1) {
    if (filtering_region->begin==0 && filtering_region->end==region_profile->pattern_length) {
      filtering_region->degree = REGION_FILTER_DEGREE_ZERO;
      region_profile->total_candidates = filtering_region->hi - filtering_region->lo;
      region_profile->max_region_length = region_profile->pattern_length;
      region_profile->avg_region_length = region_profile->pattern_length;
      region_profile->num_filtered_regions = 0;
      return true;
    }
  }
  return false;
}
/*
 * Schedule all exact
 */
void region_profile_schedule_exact(
    region_profile_t* const region_profile,
    const uint64_t candidates_threshold) {
  // Check exact matches
  if (region_profile_schedule_exact_matches(region_profile)) return;
  // Parameters
  const uint64_t num_regions = region_profile->num_filtering_regions;
  region_search_t* const filtering_region = region_profile->filtering_region;
  // Init region-profile metrics
  region_profile->total_candidates = 0;
  region_profile->num_filtered_regions = 0;
  region_profile->max_region_length = 0;
  region_profile->avg_region_length = 0;
  // Assign degree-zero to all regions
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    const uint64_t num_candidates = filtering_region[i].hi - filtering_region[i].lo;
    if (num_candidates <= candidates_threshold) {
      const uint64_t region_length = filtering_region[i].end - filtering_region[i].begin;
      // Schedule
      filtering_region[i].degree = REGION_FILTER_DEGREE_ZERO;
      ++(region_profile->num_filtered_regions);
      // Compute metrics
      region_profile->total_candidates += num_candidates;
      region_profile->max_region_length = MAX(region_profile->max_region_length,region_length);
      region_profile->avg_region_length += region_length;
    } else {
      // Discard region
      filtering_region[i].degree = REGION_FILTER_NONE;
    }
  }
  if (region_profile->num_filtered_regions > 0) {
    region_profile->avg_region_length /= region_profile->num_filtered_regions;
  }
}
void region_profile_schedule_exact_best(
    region_profile_t* const region_profile,
    const uint64_t num_regions,
    const uint64_t candidates_threshold) {
  // Check exact matches
  if (region_profile_schedule_exact_matches(region_profile)) return;
  // Allocate sorting vector
  mm_allocator_t* const mm_allocator = region_profile->mm_allocator;
  mm_allocator_push_state(mm_allocator);
  const uint64_t total_regions = region_profile->num_filtering_regions;
  region_sort_t* const region_sort = mm_allocator_calloc(mm_allocator,total_regions,region_sort_t,false);
  // Prepare sorting vetor
  region_search_t* const filtering_region = region_profile->filtering_region;
  uint64_t i;
  for (i=0;i<total_regions;++i) {
    region_sort[i].idx = i;
    region_sort[i].value = filtering_region[i].hi - filtering_region[i].lo; // Num candidates
  }
  // Sort vector
  region_sort_sort(region_sort,total_regions);
  // Init region-profile metrics
  region_profile->total_candidates = 0;
  region_profile->num_filtered_regions = 0;
  region_profile->max_region_length = 0;
  region_profile->avg_region_length = 0;
  // Schedule regions
  for (i=0;i<total_regions;++i) {
    region_search_t* const filtering_region_idx = filtering_region + region_sort[i].idx;
    const uint64_t num_candidates = filtering_region_idx->hi - filtering_region_idx->lo;
    if (i < num_regions && num_candidates <= candidates_threshold) {
      const uint64_t region_length = filtering_region_idx->end - filtering_region_idx->begin;
      // Schedule
      filtering_region_idx->degree = REGION_FILTER_DEGREE_ZERO;
      ++(region_profile->num_filtered_regions);
      // Compute metrics
      region_profile->total_candidates += num_candidates;
      region_profile->max_region_length = MAX(region_profile->max_region_length,region_length);
      region_profile->avg_region_length += region_length;
    } else {
      // Discard region
      filtering_region_idx->degree = REGION_FILTER_NONE;
    }
  }
  if (region_profile->num_filtered_regions > 0) {
    region_profile->avg_region_length /= region_profile->num_filtered_regions;
  }
  // Free
  mm_allocator_pop_state(mm_allocator);
}
/*
 * Approximate Scheduling
 */
void region_profile_schedule_approximate(
    region_profile_t* const region_profile,
    const uint64_t region_error,
    const uint64_t candidates_threshold) {
  // Check exact matches
  if (region_profile_schedule_exact_matches(region_profile)) return;
  // Parameters
  const uint64_t num_regions = region_profile->num_filtering_regions;
  region_search_t* const filtering_region = region_profile->filtering_region;
  // Init region-profile metrics
  region_profile->total_candidates = 0;
  region_profile->num_filtered_regions = 0;
  region_profile->max_region_length = 0;
  region_profile->avg_region_length = 0;
  // Assign degree-zero to all regions
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    const uint64_t num_candidates = filtering_region[i].hi - filtering_region[i].lo;
    if (num_candidates <= candidates_threshold) {
      const uint64_t region_length = filtering_region[i].end - filtering_region[i].begin;
      // Schedule
      filtering_region[i].degree = region_error+1;
      ++(region_profile->num_filtered_regions);
      // Compute metrics
      region_profile->total_candidates += filtering_region[i].hi - filtering_region[i].lo;
      region_profile->max_region_length = MAX(region_profile->max_region_length,region_length);
      region_profile->avg_region_length += region_length;
    } else {
      // Discard region
      filtering_region[i].degree = REGION_FILTER_NONE;
    }
  }
  if (region_profile->num_filtered_regions > 0) {
    region_profile->avg_region_length /= region_profile->num_filtered_regions;
  }
}












