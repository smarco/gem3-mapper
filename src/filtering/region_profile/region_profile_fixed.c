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
 *   Region-Profile module provides functions to generate an
 *   fixed/static key partition (with a fixed length size)
 */

#include "align/pattern/pattern.h"
#include "filtering/region_profile/region_profile_fixed.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "fm_index/fm_index_search.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Region Profile Partition
 */
void region_profile_partition_fixed(
    region_profile_t* const region_profile,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t region_length,
    const bool allow_gap_extending) {
  // Init & allocate region profile
  const uint64_t max_regions = region_profile->max_expected_regions;
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,max_regions);
  // Traverse the key & delimit regions of @region_length
  region_search_t* region_search = region_profile->filtering_region;
  uint64_t i, current_region_length = 0, num_filtering_regions = 0;
  region_search->begin = 0;
  for (i=0;i<key_length;++i) {
    // Cut-off
    if (num_filtering_regions >= max_regions) break;
    // Get next character & check allowed
    const uint8_t enc_char = key[i];
    if (enc_char==ENC_DNA_CHAR_N) { // Skip bad characters
      // Extend last region to cover the gap (if any)
      if (allow_gap_extending && current_region_length>0 && num_filtering_regions>0) {
        region_search_t* const last_region_search =
            region_profile->filtering_region + (num_filtering_regions-1);
        if (region_search->begin==last_region_search->end) {
          last_region_search->end = i-1;
        }
      }
      // Reset region
      region_search->begin = i+1;
      current_region_length = 0;
    }
    ++current_region_length; // Add character
    if (current_region_length == region_length) {
      // Close region
      region_search->end = i+1;
      ++region_search;
      ++num_filtering_regions;
      // Start new region
      region_search->begin = i+1;
      current_region_length = 0;
    }
  }
  // Extend last region to cover the gap (if possible)
  if (allow_gap_extending && current_region_length>0 && num_filtering_regions>0) {
    region_search_t* const last_region_search =
        region_profile->filtering_region + (num_filtering_regions-1);
    if (region_search->begin==last_region_search->end) {
      last_region_search->end = key_length;
    }
  }
  // Close profile
  region_profile->num_filtering_regions = num_filtering_regions;
}
/*
 * Fixed region-length region profile
 */
void region_profile_generate_fixed(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t region_length,
    const uint64_t region_step,
    const uint64_t region_error,
    const uint64_t max_candidates) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // Init & allocate region profile
  const uint64_t max_regions = DIV_CEIL(key_length,region_step);
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,max_regions);
  // Traverse the key & delimit regions of @region_length
  region_search_t* region_search = region_profile->filtering_region;
  uint64_t num_regions = 0;
  uint64_t region_sentinel = 0;
  while (region_sentinel+region_length <= key_length) {
    // Region
    region_search->begin = region_sentinel;
    region_search->end = region_sentinel+region_length;
    // Next
    region_sentinel += region_step;
    ++region_search;
    ++num_regions;
  }
  region_profile->num_filtering_regions = num_regions;
  // Query region profile
  region_profile_query_regions(region_profile,fm_index,key);
  // Schedule to filter all
  region_profile_schedule_approximate(region_profile,region_error,max_candidates);
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}
/*
 * Cheap k-mer selection
 */
void region_profile_generate_cks(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t region_length,
    const uint64_t num_regions,
    const uint64_t max_candidates) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  /*
   * Select those fixed-length regions with fewer
   * candidates from the static first partition
   */
  // Init & allocate region profile
  const uint64_t max_regions = DIV_CEIL(key_length,region_length);
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,max_regions);
  // Traverse the key & delimit regions of @region_length
  region_search_t* region_search = region_profile->filtering_region;
  uint64_t total_regions_generated = 0;
  uint64_t region_sentinel = 0;
  while (region_sentinel+region_length <= key_length) {
    // Region
    region_search->begin = region_sentinel;
    region_search->end = region_sentinel+region_length;
    // Next
    region_sentinel += region_length;
    ++region_search;
    ++total_regions_generated;
  }
  region_profile->num_filtering_regions = total_regions_generated;
  // Query region profile
  region_profile_query_regions(region_profile,fm_index,key);
  // Schedule to filter all
  region_profile_schedule_exact_best(region_profile,num_regions,max_candidates);
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}
/*
 * Factors region profile (divide the read in equal parts)
 */
void region_profile_generate_factors(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t num_regions,
    const uint64_t region_error,
    const uint64_t max_candidates) {
  PROFILE_START(GP_REGION_PROFILE,PROFILE_LEVEL);
  // Init & allocate region profile
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,num_regions);
  // Traverse the key & delimit regions of @region_length
  region_search_t* region_search = region_profile->filtering_region;
  uint64_t i, region_sentinel = 0, key_left = key_length;
  for (i=num_regions;i>0;--i) {
    // Begin
    region_search->begin = region_sentinel;
    // Update
    const uint64_t region_length = key_left/i;
    region_sentinel += region_length;
    key_left -= region_length;
    // End
    region_search->end = region_sentinel;
    // Next
    ++region_search;
  }
  region_profile->num_filtering_regions = num_regions;
  // Query region profile
  region_profile_query_regions(region_profile,fm_index,key);
  // Schedule to filter all
  region_profile_schedule_approximate(region_profile,region_error,max_candidates);
  PROFILE_STOP(GP_REGION_PROFILE,PROFILE_LEVEL);
}
/*
 * Display/Benchmark
 */
void region_profile_print_benchmark(
    FILE* const stream,
    const region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* key) {
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    /*
     * <chunk_length> <chunk_string> <lo> <hi> <steps_perfomed_on_search>
     * Eg. 10 ACGTACGTTG 12 134 10
     */
    // Chunk length
    const uint64_t chunk_length = region->end-region->begin;
    fprintf(stream,"%"PRIu64"\t",chunk_length);
    // Chunk string
    uint64_t i;
    for (i=region->begin;i<region->end;++i) {
      fprintf(stream,"%c",dna_decode(key[i]));
    }
    // Results
    uint64_t lo, hi, steps;
    fm_index_bsearch_debug(fm_index,key+region->begin,
        region->end-region->begin,&hi,&lo,&steps);
    fprintf(stream,"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n",lo,hi,steps);
  }
  fflush(stream);
}
