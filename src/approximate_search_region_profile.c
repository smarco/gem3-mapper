/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_region_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_region_profile.h"
#include "region_profile.h"
#include "region_profile_adaptive.h"
#include "region_profile_boost.h"
#include "region_profile_fixed.h"
#include "region_profile_schedule.h"
#include "fm_index_search.h"

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
  uint64_t total_candidates = 0;
  PROF_INC_COUNTER(GP_REGION_PROFILE_FIXED);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS,region_profile->num_filtering_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS_STANDARD,region_profile->num_standard_regions);
  PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS_UNIQUE,
      region_profile->num_filtering_regions-region_profile->num_standard_regions);
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_REGION_LENGTH,region->end-region->begin);
    PROF_ADD_COUNTER(GP_REGION_PROFILE_FIXED_REGION_CANDIDATES,(region->hi-region->lo));
    total_candidates += (region->hi-region->lo);
  }
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
void approximate_search_region_profile_boost_stats(region_profile_t* const region_profile) {
#ifdef GEM_PROFILE
  uint64_t total_candidates = 0;
  PROF_INC_COUNTER(GP_REGION_PROFILE_BOOST);
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
 * Region Partition Fixed
 */
void approximate_search_region_partition_fixed(approximate_search_t* const search) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  const archive_t* archive = search->archive;
  fm_index_t* const fm_index = search->archive->fm_index;
  pattern_t* const pattern = &search->pattern;
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
  // Generate region profile partition
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t proper_length = fm_index_get_proper_length(fm_index);
  region_profile_generate_fixed_partition(
      region_profile,key,key_length,parameters->allowed_enc,proper_length);
  /*
   * TODO Second profile
   */
  // Check region-partition result
  if (region_profile->num_filtering_regions <= 1) {
    search->processing_state = asearch_processing_state_no_regions;
  } else {
    search->processing_state = asearch_processing_state_region_partitioned;
  }
}
/*
 * Region Profile Adaptive
 */
void approximate_search_region_profile_adaptive(
    approximate_search_t* const search,
    const approximate_search_region_profile_strategy_t region_profile_strategy,
    mm_stack_t* const mm_stack) {
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
  switch (region_profile_strategy) {
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
    search->processing_state = asearch_processing_state_no_regions;
    return;
  } else if (region_profile_has_exact_matches(region_profile)) {
    const region_search_t* const first_region = region_profile->filtering_region;
    search->hi_exact_matches = first_region->hi;
    search->lo_exact_matches = first_region->lo;
    approximate_search_update_mcs(search,1);
    search->processing_state = asearch_processing_state_exact_matches;
    return;
  }
  // STATS
  switch (region_profile_strategy) {
    case region_profile_adaptive_lightweight:
      approximate_search_region_profile_lightweight_stats(region_profile);
      break;
    case region_profile_adaptive_boost:
      approximate_search_region_profile_boost_stats(region_profile);
      break;
    case region_profile_adaptive_delimit:
      approximate_search_region_profile_delimit_stats(region_profile);
      break;
    case region_profile_adaptive_limited:
    case region_profile_adaptive_recovery:
      approximate_search_region_profile_lightweight_stats(region_profile);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Buffered Copy/Retrieve
 */
void approximate_search_region_profile_buffered_print_benchmark(approximate_search_t* const search);
void approximate_search_region_profile_buffered_copy(
    approximate_search_t* const search,gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
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
    approximate_search_t* const search,gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  // Parameters
  search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  // Buffer offsets
  const uint64_t buffer_offset_begin = search->gpu_buffer_fmi_search_offset;
  // Traverse Region-Partition
  uint64_t i, num_regions_filtered = 0, total_candidates = 0;
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
    if (num_candidates <= search_parameters->gpu_filtering_threshold) {
      filtering_region->degree = REGION_FILTER_DEGREE_ZERO;
      total_candidates += num_candidates;
      ++num_regions_filtered;
    } else {
      filtering_region->degree = REGION_FILTER_NONE;
    }
  }
  region_profile->total_candidates = total_candidates; // Set total candidates
  // Set MCS & state
  const uint64_t num_wildcards = search->pattern.num_wildcards;
  approximate_search_update_mcs(search,num_regions_filtered + num_wildcards);
  search->processing_state = asearch_processing_state_region_profiled;
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  //  // Check total filtering-regions
  //  if (num_feasible_regions >= num_filtering_regions-2) {
  //    search->processing_state = asearch_processing_state_region_profiled;
  //  } else {
  //    search->processing_state = asearch_processing_state_no_regions;
  //    PROF_INC_COUNTER(GP_ASSW_REGION_PROFILE_UNSUCCESSFUL);
  //  }
  // STATS
  approximate_search_region_profile_fixed_stats(region_profile);
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

