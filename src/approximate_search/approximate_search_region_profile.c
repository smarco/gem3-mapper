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
#include "filtering/region_profile/region_profile_fixed.h"
#include "filtering/region_profile/region_profile_adaptive.h"
#include "filtering/region_profile/region_profile_optimum.h"
#include "filtering/region_profile/region_profile_mem.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "filtering/region_profile/region_profile_split.h"
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
  PROF_INC_COUNTER(GP_REGION_PROFILES_COMPUTED);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_NUM_REGIONS_FILTERED,region_profile->num_filtered_regions);
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
void approximate_search_region_profile_static_close_profile(
    approximate_search_t* const search,
    const uint64_t filtering_threshold) {
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  region_profile_schedule_exact(region_profile,filtering_threshold);
  // Check total number of filtering-regions & Set State
  const uint64_t min_num_regions = region_profile->num_filtering_regions/2;
  if (region_profile->total_candidates!=0 &&
      region_profile->num_filtered_regions > min_num_regions) {
    search->processing_state = asearch_processing_state_region_profiled;
    approximate_search_region_profile_stats(region_profile); // STATS
  } else {
    search->processing_state = asearch_processing_state_no_regions;
  }
}
void approximate_search_region_profile_adaptive_close_profile(
    approximate_search_t* const search,
    const uint64_t filtering_threshold) {
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  region_profile_schedule_exact(region_profile,filtering_threshold);
  // Check number of regions & Set State
  if (region_profile->num_filtering_regions == 0) {
    search->processing_state = asearch_processing_state_no_regions;
  } else {
    search->processing_state = asearch_processing_state_region_profiled;
  }
  approximate_search_region_profile_stats(region_profile); // STATS
}
/*
 * Region Profile Adaptive
 */
void approximate_search_region_profile(approximate_search_t* const search) {
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
  switch (parameters->region_profile_model.strategy) {
    case region_profile_fixed:
      region_profile_generate_fixed(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.region_length,
          parameters->region_profile_model.region_step,
          parameters->region_profile_model.region_error,
          parameters->region_profile_model.max_candidates);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_factor:
      region_profile_generate_factors(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.num_regions,
          parameters->region_profile_model.region_error,
          parameters->region_profile_model.max_candidates);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_CKS:
      region_profile_generate_cks(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.region_length,
          parameters->region_profile_model.num_regions,
          parameters->region_profile_model.max_candidates);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_OPS:
      region_profile_generate_optimum_fixed(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.region_length,
          parameters->region_profile_model.num_regions,
          parameters->region_profile_model.max_candidates);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_adaptive:
      region_profile_generate_adaptive(
          region_profile,fm_index,key,key_length,
          &parameters->region_profile_model,ALL,false);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_adaptive_limited:
      region_profile_generate_adaptive_limited(
          region_profile,fm_index,key,key_length,
          &parameters->region_profile_model);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_OPP:
      // Optimum
      region_profile_generate_optimum_variable(region_profile,fm_index,
          key,key_length,parameters->region_profile_model.num_regions);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_MEM:
      // MEM
      region_profile_generate_mem(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.max_candidates);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_SMEM:
      region_profile_generate_smem(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.max_candidates);
      // region_profile_print(stderr,region_profile,false);
      break;
    case region_profile_test:
      // Fixed
      region_profile_generate_fixed(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.region_length,
          parameters->region_profile_model.region_step,
          parameters->region_profile_model.region_error,
          parameters->region_profile_model.max_candidates);
      region_profile_print_pretty(stderr,region_profile,"Fixed",false);
      // CKS
      region_profile_generate_cks(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.region_length,
          parameters->region_profile_model.num_regions,
          parameters->region_profile_model.max_candidates);
      region_profile_print_pretty(stderr,region_profile,"CKS",false);
      // OPS
      region_profile_generate_optimum_fixed(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.region_length,
          parameters->region_profile_model.num_regions,
          parameters->region_profile_model.max_candidates);
      region_profile_print_pretty(stderr,region_profile,"OPS",false);
      // Factor
      region_profile_generate_factors(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.num_regions,
          parameters->region_profile_model.region_error,
          parameters->region_profile_model.max_candidates);
      region_profile_print_pretty(stderr,region_profile,"Factors",false);
      // MEM
      region_profile_generate_mem(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.max_candidates);
      region_profile_print_pretty(stderr,region_profile,"MEM",false);
      // SMEM
      region_profile_generate_smem(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.max_candidates);
      region_profile_print_pretty(stderr,region_profile,"SMEM",false);
      // Adaptive
      region_profile_generate_adaptive(
          region_profile,fm_index,key,key_length,
          &parameters->region_profile_model,ALL,false);
      region_profile_print_pretty(stderr,region_profile,"Adaptive",false);
      // Optimum
      region_profile_generate_optimum_variable(
          region_profile,fm_index,key,key_length,
          parameters->region_profile_model.num_regions+2);
      region_profile_print_pretty(stderr,region_profile,"OPP",false);
//      // Splitters
//      region_profile_splitters(
//          region_profile,fm_index,key,key_length,allowed_enc,16,
//          parameters->region_profile_model.num_regions+2,region_profile->mm_allocator);
//      region_profile_print_pretty(stderr,region_profile,"Splitters",false);
      // Disable the rest
      region_profile->num_filtering_regions = 0;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Compute K-mer frequency
  region_profile_compute_kmer_frequency(region_profile,fm_index,key,key_length);
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
  // Generate region profile partition
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t region_length = parameters->region_profile_model.region_length;
  region_profile_partition_fixed(region_profile,key,key_length,region_length,true);
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
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    // Compute region interval
    fm_index_bsearch(search->archive->fm_index,search->pattern.key+filtering_region->begin,
        filtering_region->end-filtering_region->begin,&filtering_region->hi,&filtering_region->lo);
  }
  // Close profile
  approximate_search_region_profile_static_close_profile(search,filtering_threshold);
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
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    // Fetch region search-interval
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    gpu_buffer_fmi_ssearch_get_result(gpu_buffer_fmi_ssearch,
        buffer_offset_begin+i,&filtering_region->hi,&filtering_region->lo);
    // DEBUG
#ifdef GPU_CHECK_REGION_PROFILE
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
  }
  // Compute kmer frequency
  region_profile_compute_kmer_frequency(region_profile,
      search->archive->fm_index,pattern->key,pattern->key_length);
  // Close profile
  approximate_search_region_profile_static_close_profile(search,filtering_threshold);
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
  // Retrieve buffer offsets
  uint64_t regions_offset, num_filtering_regions;
  gpu_buffer_fmi_asearch_get_result_total_regions(
      gpu_buffer_fmi_asearch,search->gpu_buffer_fmi_search_offset,
      &regions_offset,&num_filtering_regions);
  // Allocate region-profile
  region_profile_allocate_regions(region_profile,num_filtering_regions); // Allocate
  // Traverse Region-Partition
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    // Fetch region search-interval
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    gpu_buffer_fmi_asearch_get_result_region(
        gpu_buffer_fmi_asearch,regions_offset+i,
        &filtering_region->begin,&filtering_region->end,
        &filtering_region->hi,&filtering_region->lo);
  }
  region_profile->num_filtering_regions = num_filtering_regions;
  // Compute kmer frequency
  region_profile_compute_kmer_frequency(region_profile,
      search->archive->fm_index,pattern->key,pattern->key_length);
  // Close profile
  const uint64_t filtering_threshold = search_parameters->region_profile_model.region_th;
  approximate_search_region_profile_adaptive_close_profile(search,filtering_threshold);
  // DEBUG
#ifdef GPU_CHECK_REGION_PROFILE
  // Init
  mm_allocator_push_state(region_profile->mm_allocator);
  region_profile_t region_profile_cpu;
  region_profile_init(&region_profile_cpu,region_profile->pattern_length);
  region_profile_clear(&region_profile_cpu);
  region_profile_inject_mm(&region_profile_cpu,region_profile->mm_allocator);
  // Compute CPU adaptive
  region_profile_generate_adaptive(
      &region_profile_cpu,search->archive->fm_index,pattern->key,
      pattern->key_length,&search_parameters->region_profile_model,ALL,false);
  // Compare region-profiles
  if (region_profile_cmp(region_profile,&region_profile_cpu)!=0) {
    region_profile_print_pretty(stderr,region_profile,"GPU",false);
    region_profile_print_pretty(stderr,&region_profile_cpu,"CPU",false);
    gem_error_msg("ASM.Region.Profile.Buffered. Check Region-Profile failed");
  }
  mm_allocator_pop_state(region_profile->mm_allocator);
#endif
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

