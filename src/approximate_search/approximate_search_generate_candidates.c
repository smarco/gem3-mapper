/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_generate_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_generate_candidates.h"
#include "approximate_search/approximate_search_control.h"
#include "align/align.h"
#include "filtering/filtering_candidates_process.h"
#include "filtering/filtering_candidates_process_buffered.h"
#include "filtering/filtering_candidates_verify.h"
#include "filtering/filtering_candidates_align.h"
#include "filtering/region_profile.h"
#include "filtering/region_profile_adaptive.h"
#include "filtering/region_profile_schedule.h"

/*
 * Debug
 */
#define DEBUG_REGION_SCHEDULE_PRINT GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Benchmark
 */
#ifdef CUDA_BENCHMARK_GENERATE_DECODE_CANDIDATES
FILE* benchmark_decode_candidates = NULL;
#endif

/*
 * Generate Candidates from Region-Profile (Exact Regions)
 */
void approximate_search_generate_candidates_exact(
    approximate_search_t* const search,
    matches_t* const matches) {
  PROFILE_START(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = search->search_parameters;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  region_search_t* const region_search = region_profile->filtering_region;
  pattern_t* const pattern = &search->pattern;
  // Generate candidates for each region
  PROF_ADD_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS,region_profile->num_filtering_regions);
  region_profile->num_filtered_regions = 0;
  if (region_profile_has_exact_matches(region_profile)) {
    // Add exact matches
    filtering_candidates_add_region_interval(
        filtering_candidates,search_parameters,pattern,
        region_search->lo,region_search->hi,
        region_search->begin,region_search->end,
        (pattern->run_length) ? ALIGN_DISTANCE_INF : 0,
        &region_profile->candidates_limited);
    ++(region_profile->num_filtered_regions);
  } else {
    // Add all candidates from the region-profile
    uint64_t i;
    for (i=0;i<num_filtering_regions;++i) {
      PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_PROCESSED);
      // Generate exact-candidates for the region
      if (region_search[i].degree == REGION_FILTER_DEGREE_ZERO) {
        filtering_candidates_add_region_interval(
            filtering_candidates,search_parameters,pattern,
            region_search[i].lo,region_search[i].hi,
            region_search[i].begin,region_search[i].end,0,
            &region_profile->candidates_limited);
      }
      ++(region_profile->num_filtered_regions);
    }
  }
  PROFILE_STOP(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
}
/*
 * Buffered Copy/Retrieve
 */
void approximate_search_generate_candidates_buffered_print_benchmark(approximate_search_t* const search);
void approximate_search_generate_candidates_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  mm_stack_t* const mm_stack = search->mm_stack;
  // Store Buffer Positions
  const uint64_t num_filtering_positions = region_profile->total_candidates;
  search->gpu_buffer_fmi_decode_offset = gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode);
  search->gpu_buffer_fmi_decode_total = num_filtering_positions;
  filtering_position_buffered_t* gpu_filtering_positions =
      mm_stack_calloc(mm_stack,num_filtering_positions,filtering_position_buffered_t,false);
  search->gpu_filtering_positions = gpu_filtering_positions;
  search->gpu_num_filtering_positions = num_filtering_positions;
  // Copy all candidates positions
  uint64_t i;
  region_profile->num_filtered_regions = 0;
  for (i=0;i<num_filtering_regions;++i) {
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    if (filtering_region->degree==REGION_FILTER_DEGREE_ZERO) {
      uint64_t bwt_position;
      for (bwt_position=filtering_region->lo;bwt_position<filtering_region->hi;++bwt_position) {
        gpu_buffer_fmi_decode_add_query(gpu_buffer_fmi_decode,bwt_position);
        gpu_filtering_positions->source_region_begin = filtering_region->begin;
        gpu_filtering_positions->source_region_end = filtering_region->end;
        ++gpu_filtering_positions;
      }
      ++(region_profile->num_filtered_regions);
    }
  }
  PROF_ADD_COUNTER(GP_ASSW_DECODE_CANDIDATES_COPIED,region_profile->total_candidates);
  // BENCHMARK
  #ifdef CUDA_BENCHMARK_GENERATE_DECODE_CANDIDATES
  approximate_search_generate_candidates_buffered_print_benchmark(search);
  #endif
}
void approximate_search_generate_candidates_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  // Add all candidates positions
  const bool decode_text_enabled = gpu_buffer_fmi_decode->decode_text_enabled;
  const bool decode_sa_enabled = gpu_buffer_fmi_decode->decode_sa_enabled;
  uint64_t gpu_buffer_fmi_decode_offset = search->gpu_buffer_fmi_decode_offset;
  uint64_t i, gpu_filtering_positions_offset = 0;
  for (i=0;i<num_filtering_regions;++i) {
    region_search_t* const region_search = region_profile->filtering_region + i;
    if (region_search->degree==REGION_FILTER_DEGREE_ZERO) {
      // Retrieve/Decode all pending candidates
      const uint64_t pending_candidates = region_search->hi - region_search->lo;
      if (decode_text_enabled) {
        filtering_candidates_decode_text_filtering_positions_buffered(
            search->filtering_candidates,&search->pattern,
            region_search,search->gpu_filtering_positions+gpu_filtering_positions_offset,
            gpu_buffer_fmi_decode,gpu_buffer_fmi_decode_offset);
      } else if (decode_sa_enabled) {
        filtering_candidates_decode_sa_filtering_positions_buffered(
            search->filtering_candidates,&search->pattern,
            region_search,search->gpu_filtering_positions+gpu_filtering_positions_offset,
            gpu_buffer_fmi_decode,gpu_buffer_fmi_decode_offset);
      } else {
        filtering_candidates_decode_filtering_positions_buffered(
            search->filtering_candidates,&search->pattern,
            region_search,search->gpu_filtering_positions+gpu_filtering_positions_offset);
      }
      gpu_buffer_fmi_decode_offset += pending_candidates;
      gpu_filtering_positions_offset += pending_candidates;
    }
  }
  // Process all candidates
  filtering_candidates_process_candidates_buffered(search->filtering_candidates,&search->pattern,true);
  // Set state to verify-candidates
  search->processing_state = asearch_processing_state_candidates_processed;
}
/*
 * Display/Benchmark
 */
void approximate_search_generate_candidates_buffered_print_benchmark(approximate_search_t* const search) {
#ifdef CUDA_BENCHMARK_GENERATE_DECODE_CANDIDATES
  // Parameters
  region_profile_t* const region_profile = &search->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  // Prepare benchmark file
  if (benchmark_decode_candidates==NULL) {
    benchmark_decode_candidates = fopen("gem3.decode.candidates.benchmark","w+");
  }
  // Add all candidates positions
  uint64_t i;
  for (i=0;i<num_filtering_regions;++i) {
    region_search_t* const filtering_region = region_profile->filtering_region + i;
    if (filtering_region->degree==REGION_FILTER_DEGREE_ZERO) {
      uint64_t bwt_position;
      for (bwt_position=filtering_region->begin;bwt_position<filtering_region->end;++bwt_position) {
        fm_index_decode_print_benchmark(benchmark_decode_candidates,search->archive->fm_index,bwt_position);
      }
    }
  }
#endif
}

