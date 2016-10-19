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
 *   Approximate-String-Matching (ASM) module to produce
 *   region-profiles (i.e. key-partition with few candidates)
 */

#include "approximate_search/approximate_search_stages.h"
#include "approximate_search/approximate_search_region_profile.h"
#include "filtering/region_profile/region_profile.h"
#include "filtering/region_profile/region_profile_adaptive.h"
#include "filtering/region_profile/region_profile_fixed.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "fm_index/fm_index_search.h"

/*
 * Debug
 */
#define DEBUG_REGION_PROFILE_PRINT          GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Benchmark
 */
#ifdef CUDA_BENCHMARK_GENERATE_REGION_PROFILE
FILE* benchmark_region_profile = NULL;
void approximate_search_region_profile_buffered_print_benchmark(approximate_search_t* const search);
#endif

/*
 * Region Profile Stats
 */
void approximate_search_region_profile_stats(region_profile_t* const region_profile) {
#ifdef GEM_PROFILE
  uint64_t total_candidates = 0;
  PROF_INC_COUNTER(GP_REGION_PROFILE);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_NUM_REGIONS,region_profile->num_filtering_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_REGION_LENGTH,region->end-region->begin);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_TOTAL_CANDIDATES,total_candidates);
#endif
}
/*
 * Region Profile Utils
 */
void approximate_search_region_profile_close_region(
    region_search_t* const filtering_region,
    const uint64_t filtering_threshold,
    uint64_t* const total_candidates,
    uint64_t* const num_regions_filtered,
    uint64_t* const avg_region_length,
    uint64_t* const max_region_length) {
  // Check number of candidates
  const uint64_t num_candidates = filtering_region->hi - filtering_region->lo;
  if (num_candidates <= filtering_threshold) {
    filtering_region->degree = REGION_FILTER_DEGREE_ZERO;
    *total_candidates += num_candidates;
    ++(*num_regions_filtered);
    // Accumulate Mappability
    const uint64_t region_length = filtering_region->end - filtering_region->begin;
    *max_region_length = MAX(*max_region_length,region_length);
    *avg_region_length += region_length;
  } else {
    filtering_region->degree = REGION_FILTER_NONE;
  }
}
void approximate_search_region_profile_static_close_profile(
    approximate_search_t* const search,
    const uint64_t num_filtering_regions,
    const uint64_t num_regions_filtered,
    const uint64_t total_candidates,
    const uint64_t max_region_length) {
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  // Check total number of filtering-regions
  const uint64_t key_length = search->pattern.key_length;
  const uint64_t min_num_regions = num_filtering_regions/2;
  if (total_candidates!=0 && num_regions_filtered > min_num_regions) {
    // Close region profile
    region_profile->total_candidates = total_candidates;
    region_profile->max_region_length = max_region_length;
    region_profile->avg_region_length /= region_profile->num_filtering_regions;
    // Set State
    search->processing_state = asearch_processing_state_region_profiled;
    // STATS
    approximate_search_region_profile_stats(region_profile);
  } else {
    // Set State
    region_profile->num_filtering_regions = 0;
    region_profile->total_candidates = 0;
    region_profile->max_region_length = key_length;
    region_profile->avg_region_length = key_length;
    search->processing_state = asearch_processing_state_no_regions;
  }
}
void approximate_search_region_profile_adaptive_close_profile(
    approximate_search_t* const search,
    const uint64_t num_regions,
    const uint64_t num_regions_filtered,
    uint64_t total_candidates,
    const uint64_t avg_region_length,
    const uint64_t max_region_length) {
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  // Check exact matches
  const uint64_t key_length = search->pattern.key_length;
  if (num_regions == 1 && num_regions_filtered == 0) {
    region_search_t* const first_region = region_profile->filtering_region;
    if (first_region->begin==0 && first_region->end==key_length) {
      region_profile->num_filtering_regions = 1;
      total_candidates += first_region->hi - first_region->lo;
      first_region->degree = REGION_FILTER_DEGREE_ZERO;
    } else {
      region_profile->num_filtering_regions = 0;
    }
  } else {
    region_profile->num_filtering_regions = num_regions;
  }
  // Check number of regions
  if (region_profile->num_filtering_regions == 0) {
    // Set State & close region profile
    search->processing_state = asearch_processing_state_no_regions;
    region_profile->total_candidates = 0;
    region_profile->max_region_length = key_length;
    region_profile->avg_region_length = key_length;
  } else {
    // Set State & close region profile
    search->processing_state = asearch_processing_state_region_profiled;
    region_profile->total_candidates = total_candidates;
    region_profile->max_region_length = max_region_length;
    region_profile->avg_region_length = avg_region_length/region_profile->num_filtering_regions;
  }
  // STATS
  approximate_search_region_profile_stats(region_profile);
}
/*
 * Region Profile Adaptive
 */
void approximate_search_region_profile_adaptive(
    approximate_search_t* const search,
    const region_profile_strategy_t strategy) {
  // Parameters
  const search_parameters_t* const parameters = search->search_parameters;
  fm_index_t* const fm_index = search->archive->fm_index;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
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
  // Compute the region profile
  switch (strategy) {
    case region_profile_adaptive:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->region_profile_model,ALL,false);
      break;
    case region_profile_adaptive_limited:
      region_profile_generate_adaptive_limited(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->region_profile_model,
          search->current_max_complete_error+1);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Compute K-mer frequency
  region_profile_compute_kmer_frequency(region_profile,
      fm_index,key,key_length,parameters->allowed_enc);
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) { region_profile_print(stderr,region_profile,false); }
  // Check Zero-Region
  if (region_profile->num_filtering_regions==0) {
    search->processing_state = asearch_processing_state_no_regions;
    return;
  }
  // STATS
  approximate_search_region_profile_stats(region_profile);
}
/*
 * Region Partition Fixed
 */
void approximate_search_region_profile_static_partition(approximate_search_t* const search) {
  // Parameters
  const search_parameters_t* const parameters = search->search_parameters;
  pattern_t* const pattern = &search->pattern;
  const uint8_t* key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const bool* const allowed_enc = parameters->allowed_enc;
  // Generate region profile partition
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t min_region_length = parameters->region_profile_model.region_length;
  region_profile_generate_fixed_partition(region_profile,key,key_length,allowed_enc,min_region_length);
  // Check region-partition result
  if (region_profile->num_filtering_regions <= 1) {
    search->processing_state = asearch_processing_state_no_regions;
  } else {
    search->processing_state = asearch_processing_state_region_partitioned;
  }
}
void approximate_search_region_profile_static_compute(approximate_search_t* const search) {
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  const uint64_t filtering_threshold = search_parameters->region_profile_model.region_th;
  // Traverse Region-Partition
  uint64_t num_regions_filtered = 0, total_candidates = 0;
  uint64_t avg_region_length = 0, max_region_length = 0;
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    // Compute region interval
    fm_index_bsearch(search->archive->fm_index,search->pattern.key+filtering_region->begin,
        filtering_region->end-filtering_region->begin,&filtering_region->hi,&filtering_region->lo);
    // Close region
    approximate_search_region_profile_close_region(
        filtering_region,filtering_threshold,&total_candidates,
        &num_regions_filtered,&avg_region_length,&max_region_length);
  }
  // Close profile
  approximate_search_region_profile_static_close_profile(
      search,num_filtering_regions,num_regions_filtered,
      total_candidates,max_region_length);
}
/*
 * Static Buffered Copy/Retrieve
 */
void approximate_search_region_profile_static_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  // Parameters
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  // Store Buffer Position
  search->gpu_buffer_fmi_search_offset = gpu_buffer_fmi_ssearch_get_num_queries(gpu_buffer_fmi_ssearch);
  search->gpu_buffer_fmi_search_total = num_filtering_regions;
  // Traverse Region-Profile
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    gpu_buffer_fmi_ssearch_add_query(gpu_buffer_fmi_ssearch,pattern,
        region_profile->filtering_region[i].begin,region_profile->filtering_region[i].end);
  }
  // BENCHMARK
#ifdef CUDA_BENCHMARK_GENERATE_REGION_PROFILE
  approximate_search_region_profile_buffered_print_benchmark(search);
#endif
}
void approximate_search_region_profile_static_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  const uint64_t filtering_threshold = search_parameters->region_profile_model.region_th;
  // Buffer offsets
  const uint64_t buffer_offset_begin = search->gpu_buffer_fmi_search_offset;
  // Traverse Region-Partition
  uint64_t num_regions_filtered = 0, total_candidates = 0;
  uint64_t avg_region_length = 0, max_region_length = 0;
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    // Fetch region search-interval
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    gpu_buffer_fmi_ssearch_get_result(gpu_buffer_fmi_ssearch,
        buffer_offset_begin+i,&filtering_region->hi,&filtering_region->lo);
    // DEBUG
#ifdef CUDA_CHECK_BUFFERED_REGION_PROFILE
    uint64_t hi, lo;
    fm_index_bsearch(
        search->archive->fm_index,search->pattern.key+filtering_region->begin,
        filtering_region->end-filtering_region->begin,&hi,&lo);
    gem_cond_error_msg(
             (hi-lo!=0 || filtering_region->hi-filtering_region->lo!=0) &&
             (filtering_region->hi!=hi || filtering_region->lo!=lo),
             "ASM.Region.Profile.Buffered. Check Region-Profile failed (hi::%lu!=%lu)(lo::%lu!=%lu)",
             filtering_region->hi,hi,filtering_region->lo,lo);
#endif
    // Close region
    approximate_search_region_profile_close_region(
        filtering_region,filtering_threshold,&total_candidates,
        &num_regions_filtered,&avg_region_length,&max_region_length);
  }
  // Compute kmer frequency
  region_profile_compute_kmer_frequency(
      region_profile,search->archive->fm_index,pattern->key,
      pattern->key_length,search_parameters->allowed_enc);
  // Close profile
  approximate_search_region_profile_static_close_profile(
      search,num_filtering_regions,num_regions_filtered,
      total_candidates,max_region_length);
  // DEBUG
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) { region_profile_print(stderr,region_profile,false); }
}
/*
 * Adaptive Buffered Copy/Retrieve
 */
void approximate_search_region_profile_adaptive_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  // Parameters
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Add query
  search->gpu_buffer_fmi_search_offset = gpu_buffer_fmi_asearch_add_query(
      gpu_buffer_fmi_asearch,pattern,region_profile->max_expected_regions);
}
void approximate_search_region_profile_adaptive_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t filtering_threshold = search_parameters->region_profile_model.region_th;
  // Retrieve buffer offsets
  uint64_t regions_offset, num_regions;
  gpu_buffer_fmi_asearch_get_result_total_regions(
      gpu_buffer_fmi_asearch,search->gpu_buffer_fmi_search_offset,
      &regions_offset,&num_regions);
  // Allocate region-profile
  region_profile_allocate_regions(region_profile,num_regions); // Allocate
  // Traverse Region-Partition
  uint64_t num_regions_filtered = 0, total_candidates = 0;
  uint64_t avg_region_length = 0, max_region_length = 0;
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    // Fetch region search-interval
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    gpu_buffer_fmi_asearch_get_result_region(
        gpu_buffer_fmi_asearch,regions_offset+i,
        &filtering_region->begin,&filtering_region->end,
        &filtering_region->hi,&filtering_region->lo);
    // Close region
    approximate_search_region_profile_close_region(
        filtering_region,filtering_threshold,&total_candidates,
        &num_regions_filtered,&avg_region_length,&max_region_length);
  }
  // Compute kmer frequency
  region_profile_compute_kmer_frequency(
      region_profile,search->archive->fm_index,pattern->key,
      pattern->key_length,search_parameters->allowed_enc);
  // Close profile
  approximate_search_region_profile_adaptive_close_profile(
      search,num_regions,num_regions_filtered,
      total_candidates,avg_region_length,max_region_length);
  // DEBUG
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) { region_profile_print(stderr,region_profile,false); }
}
/*
 * Benchmark
 */
#ifdef CUDA_BENCHMARK_GENERATE_REGION_PROFILE
void approximate_search_region_profile_buffered_print_benchmark(approximate_search_t* const search) {
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  // Prepare benchmark file
  if (benchmark_region_profile==NULL) {
    benchmark_region_profile = fopen("gem3.region.profile.benchmark","w+");
  }
  // Print Region-Profile benchmark
  region_profile_print_benchmark(benchmark_region_profile,
      region_profile,search->archive->fm_index,search->pattern.key);
}
#endif

