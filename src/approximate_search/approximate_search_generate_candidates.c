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
 *   Approximate-String-Matching (ASM) module to generate filtering-candidates
 *   from a region-profile/key-partition
 */

#include "approximate_search/approximate_search_generate_candidates.h"
#include "align/alignment.h"
#include "filtering/candidates/filtering_candidates_buffered_process.h"
#include "filtering/candidates/filtering_candidates_accessors.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_align.h"
#include "filtering/region_profile/region_profile.h"
#include "filtering/region_profile/region_profile_adaptive.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "neighborhood_search/nsearch_levenshtein.h"

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
void approximate_search_generate_candidates_limit_exact_matches(
    approximate_search_t* const search) {
  // Check exact matches (limit the number of matches)
  region_profile_t* const region_profile = &search->region_profile;
  search_parameters_t* const search_parameters = search->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  region_search_t* const filtering_region = region_profile->filtering_region;
  const uint64_t total_candidates = filtering_region->hi - filtering_region->lo;
  if (total_candidates > select_parameters->max_searched_matches) {
    search->num_limited_exact_matches = total_candidates - select_parameters->max_searched_matches;
    filtering_region->hi = filtering_region->lo + select_parameters->max_searched_matches;
    region_profile->total_candidates = select_parameters->max_searched_matches;
  }
}
void approximate_search_generate_candidates(
    approximate_search_t* const search) {
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
    // Limit exact matches
    approximate_search_generate_candidates_limit_exact_matches(search);
    // Add exact matches
    filtering_candidates_add_positions_from_interval(
        filtering_candidates,search_parameters,pattern,
        region_search->lo,region_search->hi,
        region_search->begin,region_search->end,
        (pattern->run_length) ? ALIGN_DISTANCE_UNKNOWN : 0);
    ++(region_profile->num_filtered_regions);
  } else {
    // Add all candidates from the region-profile
    uint64_t i;
    for (i=0;i<num_filtering_regions;++i) {
      PROF_INC_COUNTER(GP_AS_GENERATE_CANDIDATES_PROCESSED);
      // Generate exact-candidates for the region
      const uint64_t degree = region_search[i].degree;
      if (degree == REGION_FILTER_NONE) continue;
      if (degree == REGION_FILTER_DEGREE_ZERO) {
        filtering_candidates_add_positions_from_interval(
            filtering_candidates,search_parameters,pattern,
            region_search[i].lo,region_search[i].hi,
            region_search[i].begin,region_search[i].end,0);
      } else {
        nsearch_levenshtein_base(
            search,pattern->key+region_search[i].begin,
            region_search[i].end-region_search[i].begin,
            region_search[i].begin,degree-1);
      }
      ++(region_profile->num_filtered_regions);
    }
  }
  PROFILE_STOP(GP_AS_GENERATE_CANDIDATES,PROFILE_LEVEL);
}
/*
 * Buffered Copy/Retrieve
 */
#ifdef CUDA_BENCHMARK_GENERATE_DECODE_CANDIDATES
void approximate_search_generate_candidates_buffered_print_benchmark(approximate_search_t* const search);
#endif
void approximate_search_generate_candidates_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  // Parameters
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_buffered_t* const filtering_candidates_buffered = &search->filtering_candidates_buffered;
  // Generate candidates
  approximate_search_generate_candidates(search);
  const uint64_t num_filtering_positions = filtering_candidates_get_num_positions(filtering_candidates);
  if (num_filtering_positions==0) {
    filtering_candidates_buffered->num_positions = 0;
    return;
  }
  // Allocate buffered positions
  filtering_position_t** const positions = filtering_candidates_get_positions(filtering_candidates);
  filtering_candidates_buffered_allocate_positions(filtering_candidates_buffered,num_filtering_positions);
  filtering_position_t** const positions_buffered = filtering_candidates_buffered->positions;
  // Copy candidates positions to GPU buffer
  uint64_t i;
  search->gpu_buffer_fmi_decode_offset = gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode);
  for (i=0;i<num_filtering_positions;++i) {
    positions_buffered[i] = positions[i];
    gpu_buffer_fmi_decode_add_query(gpu_buffer_fmi_decode,positions[i]->region_index_position);
  }
  filtering_candidates_buffered->num_positions = num_filtering_positions;
  // Clear positions in filtering candidates (stoted buffered)
  filtering_candidates_clear_positions(filtering_candidates,false);
  PROF_ADD_COUNTER(GP_ASSW_DECODE_CANDIDATES_COPIED,num_filtering_positions);
  // BENCHMARK
  #ifdef CUDA_BENCHMARK_GENERATE_DECODE_CANDIDATES
  approximate_search_generate_candidates_buffered_print_benchmark(search);
  #endif
}
void approximate_search_generate_candidates_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  // Parameters
  const uint64_t gpu_buffer_fmi_decode_offset = search->gpu_buffer_fmi_decode_offset;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_buffered_t* const filtering_candidates_buffered = &search->filtering_candidates_buffered;
  // Check filtering positions
  const uint64_t num_filtering_positions = filtering_candidates_buffered->num_positions;
  if (num_filtering_positions==0) {
    search->processing_state = asearch_processing_state_candidates_processed;
    return; // No filtering positions
  }
  // Retrieve decoded positions
  if (gpu_buffer_fmi_decode->decode_text_enabled) {
    filtering_candidates_buffered_decode_text_positions(
        filtering_candidates,filtering_candidates_buffered,
        &search->pattern,gpu_buffer_fmi_decode,gpu_buffer_fmi_decode_offset);
  } else if (gpu_buffer_fmi_decode->decode_sa_enabled) {
    filtering_candidates_buffered_decode_sampled_sa_positions(
        filtering_candidates,filtering_candidates_buffered,
        &search->pattern,gpu_buffer_fmi_decode,gpu_buffer_fmi_decode_offset);
  } else {
    filtering_candidates_buffered_decode_sa_positions(
        filtering_candidates,filtering_candidates_buffered,
        &search->pattern);
  }
  // Process all candidates
  filtering_candidates_buffered_process_candidates(search->filtering_candidates,&search->pattern,true);
  // Free buffered positions
  filtering_candidates_buffered_free_positions(&search->filtering_candidates_buffered);
  // Set state to verify-candidates
  search->processing_state = asearch_processing_state_candidates_processed;
}
/*
 * Display/Benchmark
 */
#ifdef CUDA_BENCHMARK_GENERATE_DECODE_CANDIDATES
void approximate_search_generate_candidates_buffered_print_benchmark(approximate_search_t* const search) {
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
}
#endif

