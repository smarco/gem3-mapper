/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_region_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_region_profile.h"
#include "approximate_search/approximate_search_filtering_stages.h"
#include "filtering/region_profile.h"
#include "filtering/region_profile_adaptive.h"
#include "filtering/region_profile_fixed.h"
#include "filtering/region_profile_schedule.h"
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
#endif

/*
 * Region Profile Stats
 */
void approximate_search_region_profile_fixed_stats(region_profile_t* const region_profile) {
#ifdef GEM_PROFILE
  uint64_t total_candidates = 0, num_elegible_regions = 0;
  PROF_INC_COUNTER(GP_REGION_PROFILE_FIXED);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    if (region->degree!=REGION_FILTER_DEGREE_ZERO) continue;
    ++num_elegible_regions;
    PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_REGION_LENGTH,region->end-region->begin);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS,num_elegible_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS_STANDARD,num_elegible_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_TOTAL_CANDIDATES,total_candidates);
#endif
}
void approximate_search_region_profile_lightweight_stats(region_profile_t* const region_profile) {
#ifdef GEM_PROFILE
  uint64_t total_candidates = 0;
  PROF_INC_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT);
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
#endif
}
void approximate_search_region_profile_heavyweight_stats(region_profile_t* const region_profile) {
#ifdef GEM_PROFILE
  uint64_t total_candidates = 0;
  PROF_INC_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_NUM_REGIONS_STANDARD,region_profile->num_standard_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_NUM_REGIONS_UNIQUE,
      region_profile->num_filtering_regions-region_profile->num_standard_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_REGION_LENGTH,region->end-region->begin);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
  PROF_ADD_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_TOTAL_CANDIDATES,total_candidates);
#endif
}
void approximate_search_region_profile_delimit_stats(region_profile_t* const region_profile) {
#ifdef GEM_PROFILE
  uint64_t total_candidates = 0;
  PROF_INC_COUNTER(GP_REGION_PROFILE_DELIMIT);
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
#endif
}
/*
 * Region Profile Adaptive
 */
void approximate_search_region_profile_adaptive(
    approximate_search_t* const search,
    const region_profile_strategy_t strategy,
    mm_stack_t* const mm_stack) {
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
    case region_profile_adaptive_lightweight:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_lightweight,ALL,false);
      break;
    case region_profile_adaptive_heavyweight:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_heavyweight,ALL,false);
      break;
    case region_profile_adaptive_limited:
      region_profile_generate_adaptive_limited(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_lightweight,search->max_complete_error+1);
      break;
    case region_profile_adaptive_delimit:
      region_profile_generate_adaptive(region_profile,fm_index,key,key_length,
          parameters->allowed_enc,&parameters->rp_delimit,UINT64_MAX,true);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Set Metrics
  approximate_search_metrics_set_max_region_length(&search->metrics,region_profile->max_region_length);
  approximate_search_metrics_set_num_zero_regions(&search->metrics,region_profile->num_zero_regions);
  approximate_search_metrics_set_mappability(&search->metrics,
      region_profile->mappability_p,region_profile->mappability_2p);
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) { region_profile_print(stderr,region_profile,false); }
  // Check Zero-Region
  if (region_profile->num_filtering_regions==0) {
    approximate_search_update_mcs(search,pattern->num_wildcards);
    search->processing_state = asearch_processing_state_no_regions;
    return;
  }
  // STATS
  switch (strategy) {
    case region_profile_adaptive_limited:
    case region_profile_adaptive_lightweight:
      approximate_search_region_profile_lightweight_stats(region_profile);
      break;
    case region_profile_adaptive_heavyweight:
      approximate_search_region_profile_heavyweight_stats(region_profile);
      break;
    case region_profile_adaptive_delimit:
      approximate_search_region_profile_delimit_stats(region_profile);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Region Partition Fixed
 */
void approximate_search_region_partition_fixed(approximate_search_t* const search) {
  // Parameters
  const search_parameters_t* const parameters = search->search_parameters;
  fm_index_t* const fm_index = search->archive->fm_index;
  pattern_t* const pattern = &search->pattern;
  const uint8_t* key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  // Generate region profile partition
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t proper_length = fm_index_get_proper_length(fm_index);
  region_profile_generate_fixed_partition(
      region_profile,key,key_length,parameters->allowed_enc,proper_length);
  // Check region-partition result
  if (region_profile->num_filtering_regions <= 1) {
    search->processing_state = asearch_processing_state_no_regions;
  } else {
    search->processing_state = asearch_processing_state_region_partitioned;
  }
}
/*
 * Buffered Copy/Retrieve
 */
void approximate_search_region_profile_buffered_print_benchmark(approximate_search_t* const search);
void approximate_search_region_profile_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  // Parameters
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  // Store Buffer Position
  search->gpu_buffer_fmi_search_offset = gpu_buffer_fmi_search_get_num_queries(gpu_buffer_fmi_search);
  search->gpu_buffer_fmi_search_total = num_filtering_regions;
  // Traverse Region-Profile
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    gpu_buffer_fmi_search_add_query(gpu_buffer_fmi_search,pattern,
        region_profile->filtering_region[i].begin,region_profile->filtering_region[i].end);
  }
  // BENCHMARK
#ifdef CUDA_BENCHMARK_GENERATE_REGION_PROFILE
  approximate_search_region_profile_buffered_print_benchmark(search);
#endif
}
void approximate_search_region_profile_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  const uint64_t filtering_threshold = search_parameters->rp_heavyweight.region_th;
  // Buffer offsets
  const uint64_t buffer_offset_begin = search->gpu_buffer_fmi_search_offset;
  // Traverse Region-Partition
  uint64_t num_regions_filtered = 0, num_zero_regions = 0, max_region_length = 0;
  uint64_t i, total_candidates = 0;
  for (i=0;i<num_filtering_regions;++i) {
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    gpu_buffer_fmi_search_get_result(gpu_buffer_fmi_search,
        buffer_offset_begin+i,&filtering_region->hi,&filtering_region->lo);
    // DEBUG
#ifdef CUDA_CHECK_BUFFERED_REGION_PROFILE
      uint64_t hi, lo;
      fm_index_bsearch(search->archive->fm_index,search->pattern.key+filtering_region->begin,
          filtering_region->end-filtering_region->begin,&hi,&lo);
      gem_cond_error_msg(
               (hi-lo!=0 || filtering_region->hi-filtering_region->lo!=0) &&
               (filtering_region->hi!=hi || filtering_region->lo!=lo),
               "ASM.Region.Profile.Buffered. Check Region-Profile failed (hi::%lu!=%lu)(lo::%lu!=%lu)",
               filtering_region->hi,hi,filtering_region->lo,lo);
#endif
    // Check number of candidates
    const uint64_t num_candidates = filtering_region->hi - filtering_region->lo;
    if (num_candidates <= filtering_threshold) {
      filtering_region->degree = REGION_FILTER_DEGREE_ZERO;
      total_candidates += num_candidates;
      ++num_regions_filtered;
      if (num_candidates==0) ++num_zero_regions;
    } else {
      filtering_region->degree = REGION_FILTER_NONE;
    }
    // Accumulate Mappability
    const uint64_t region_length = filtering_region->end - filtering_region->begin;
    max_region_length = MAX(max_region_length,region_length);
    if (num_candidates>0) region_profile->mappability_p += log2((double)num_candidates);
  }
  // Check total number of filtering-regions
  if (total_candidates!=0 && (num_regions_filtered-num_zero_regions) > num_filtering_regions/2) {
    // Close region profile
    region_profile->mappability_p /= (double)(2*num_filtering_regions);
    region_profile->mappability_2p = 0.0;
    region_profile->total_candidates = total_candidates;
    region_profile->max_region_length = max_region_length;
    region_profile->num_zero_regions = num_zero_regions;
    // Set State
    search->processing_state = asearch_processing_state_region_profiled;
    // Set Metrics
    approximate_search_metrics_set_max_region_length(&search->metrics,max_region_length);
    approximate_search_metrics_set_num_zero_regions(&search->metrics,num_zero_regions);
    approximate_search_metrics_set_mappability(&search->metrics,
        region_profile->mappability_p,region_profile->mappability_2p);
    // STATS
    approximate_search_region_profile_fixed_stats(region_profile);
  } else {
    // Set State
    region_profile->num_filtering_regions = 0;
    region_profile->total_candidates = 0;
    search->processing_state = asearch_processing_state_no_regions;
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) {
    region_profile_print(stderr,region_profile,false);
  }
}
void approximate_search_region_profile_buffered_recompute(approximate_search_t* const search) {
  // Re-Compute region profile
  search->processing_state = asearch_processing_state_begin;
  approximate_search_region_profile_adaptive(search,region_profile_adaptive_lightweight,search->mm_stack);
  if (search->processing_state==asearch_processing_state_no_regions) {
    approximate_search_update_mcs(search,search->pattern.num_wildcards); // Set MCS
    return;
  }
  // Schedule exact-candidates
  const search_parameters_t* const search_parameters = search->search_parameters;
  region_profile_schedule_filtering_fixed(&search->region_profile,ALL,
      REGION_FILTER_DEGREE_ZERO,search_parameters->filtering_threshold);
  // Set State
  search->processing_state = asearch_processing_state_region_profiled;
  // DEBUG
  gem_cond_debug_block(DEBUG_REGION_PROFILE_PRINT) {
    region_profile_print(stderr,&search->region_profile,false);
  }
}
void approximate_search_region_profile_buffered_print_benchmark(approximate_search_t* const search) {
#ifdef CUDA_BENCHMARK_GENERATE_REGION_PROFILE
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  // Prepare benchmark file
  if (benchmark_region_profile==NULL) {
    benchmark_region_profile = fopen("gem3.region.profile.benchmark","w+");
  }
  // Print Region-Profile benchmark
  region_profile_print_benchmark(benchmark_region_profile,
      region_profile,search->archive->fm_index,search->pattern.key);
#endif
}

