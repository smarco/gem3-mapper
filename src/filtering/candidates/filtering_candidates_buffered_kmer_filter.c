/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/candidates/filtering_candidates_buffered_kmer_filter.h"
#include "align/alignment.h"
#include "filtering/region/filtering_region_verify.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Constants
 */
#define FILTERING_CANDIDATES_KMER_FILTER_MIN_KEY_LENGTH 250000

/*
 * Kmer-Filter Buffered Add
 */
void filtering_candidates_buffered_kmer_filter_prepare_buffers(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    uint64_t* const gpu_buffer_kmer_filter_offset) {
  // Parameters
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  // Allocate buffered filtering regions
  filtering_candidates_buffered_allocate_regions(filtering_candidates_buffered,num_filtering_regions);
  filtering_candidates_buffered->num_regions = num_filtering_regions;
  // Set buffer offset
  *gpu_buffer_kmer_filter_offset = gpu_buffer_kmer_filter_get_num_candidates(gpu_buffer_kmer_filter);
}
void filtering_candidates_buffered_kmer_filter_add_pattern(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  // Add pattern
  pattern_tiled_t* const pattern_tiled = &pattern->pattern_tiled;
  gpu_buffer_kmer_filter_add_pattern(gpu_buffer_kmer_filter,&pattern_tiled->kmer_filter_nway);
}
void filtering_candidates_buffered_kmer_filter_add_filtering_regions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  // Parameters
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  filtering_region_t** const region_buffered = filtering_candidates_buffered->regions;
  pattern_tiled_t* const pattern_tiled = &pattern->pattern_tiled;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Add all candidates
  filtering_region_t** const regions_in = filtering_candidates_get_regions(filtering_candidates);
  uint64_t candidate_pos;
  for (candidate_pos=0;candidate_pos<num_filtering_regions;++candidate_pos) {
    // Filtering region
    filtering_region_t* const filtering_region = regions_in[candidate_pos];
    region_buffered[candidate_pos] = filtering_region;
    // Filter out exact-matches & key-trimmed regions
    if (filtering_region->alignment.distance_min_bound==0 ||
        filtering_region->key_trimmed ||
        pattern->key_length < FILTERING_CANDIDATES_KMER_FILTER_MIN_KEY_LENGTH ||
        !pattern_tiled->kmer_filter_nway.enabled) {
      continue; // Next
    }
    // Prepare Kmer-filter (candidate)
    const uint64_t text_length = filtering_region->text_end_position - filtering_region->text_begin_position;
    kmer_counting_prepare_tiling(&pattern_tiled->kmer_filter_nway,text_length,max_error);
    // Add Candidate
    gpu_buffer_kmer_filter_add_candidate(
        gpu_buffer_kmer_filter,&pattern_tiled->kmer_filter_nway,
        filtering_region->text_begin_position);
  }
  // Clear filtering regions (regions buffered)
  filtering_candidates_clear_regions(filtering_candidates,false);
}
void filtering_candidates_buffered_kmer_filter_add(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    uint64_t* const gpu_buffer_kmer_filter_offset) {
  // Check number of pending filtering-regions
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  if (num_filtering_regions==0) {
    filtering_candidates_buffered->num_regions = 0;
    return;
  }
  // Prepare Buffers
  filtering_candidates_buffered_kmer_filter_prepare_buffers(
      filtering_candidates,filtering_candidates_buffered,
      gpu_buffer_kmer_filter,gpu_buffer_kmer_filter_offset);
  // Add the pattern
  filtering_candidates_buffered_kmer_filter_add_pattern(
      filtering_candidates,pattern,gpu_buffer_kmer_filter);
  // Add all filtering regions
  filtering_candidates_buffered_kmer_filter_add_filtering_regions(
      filtering_candidates,filtering_candidates_buffered,
      pattern,gpu_buffer_kmer_filter);
}
/*
 * Kmer-Filter Computation (emulation/checker)
 */
uint64_t filtering_candidates_buffered_kmer_filter_compute_alignment(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  mm_allocator_t* const mm_allocator = filtering_candidates->mm_allocator;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Retrieve text
  const uint64_t candidate_position = filtering_region->text_begin_position;
  const uint64_t candidate_length =
      filtering_region->text_end_position - filtering_region->text_begin_position;
  archive_text_retrieve(
      archive_text,candidate_position,candidate_length,
      false,false,&filtering_region->text_trace,
      filtering_candidates->mm_allocator);
  text_trace_t* const text_trace = &filtering_region->text_trace;
  // Check kmer-filter compiled
  const uint64_t distance_min_bound = kmer_counting_min_bound_nway(
      &pattern->pattern_tiled.kmer_filter_nway,
      text_trace->text,text_trace->text_length,
      max_error,mm_allocator);
  // Return
  return distance_min_bound;
}
/*
 * Kmer-Filter Buffered Retrieve
 */
void filtering_candidates_buffered_kmer_filter_retrieve_region(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    uint64_t* const gpu_buffer_offset) {
  // Parameters
  alignment_t* const alignment = &filtering_region->alignment;
  pattern_tiled_t* const pattern_tiled = &pattern->pattern_tiled;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Detect exact-matches
  if (filtering_region->alignment.distance_min_bound == 0) {
    vector_insert(filtering_candidates->filtering_regions,filtering_region,filtering_region_t*);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
    return;
  }
  // Detect trimmed matches
  if (filtering_region->key_trimmed ||
      pattern->key_length < FILTERING_CANDIDATES_KMER_FILTER_MIN_KEY_LENGTH ||
      !pattern_tiled->kmer_filter_nway.enabled) {
    alignment->distance_min_bound = ALIGN_DISTANCE_UNKNOWN;
    vector_insert(filtering_candidates->filtering_regions,filtering_region,filtering_region_t*);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
    return;
  }
  // Retrieve min-bound
  uint64_t min_distance_bound;
  if (gpu_buffer_kmer_filter->kmer_filter_enabled) {
    min_distance_bound = gpu_buffer_kmer_filter_get_min_edit_bound(
        gpu_buffer_kmer_filter,*gpu_buffer_offset,
        pattern_tiled->kmer_filter_nway.num_tiles);
    *gpu_buffer_offset += pattern_tiled->kmer_filter_nway.num_tiles;
    // DEBUG
    #ifdef GPU_CHECK_KMER_FILTER
    const uint64_t min_distance_bound_check =
        filtering_candidates_buffered_kmer_filter_compute_alignment(
            filtering_candidates,filtering_region,pattern);
    if (min_distance_bound_check != min_distance_bound) {
      fprintf(stderr,"GPU.Kmer-Filter. Difference detected (CPU=%lu;GPU=%lu)\n",
          min_distance_bound_check,min_distance_bound);
    }
    #endif
  } else {
    min_distance_bound = 0; // Unknown
//        filtering_candidates_buffered_kmer_filter_compute_alignment(
//            filtering_candidates,filtering_region,pattern);
  }
  if (min_distance_bound <= max_error) {
    alignment->distance_min_bound = ALIGN_DISTANCE_UNKNOWN;
    vector_insert(filtering_candidates->filtering_regions,filtering_region,filtering_region_t*);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
  } else {
    // Discarded candidate
    filtering_region->status = filtering_region_verified_discarded;
    alignment->distance_min_bound = ALIGN_DISTANCE_INF;
    // Add to discarded candidates (buffered)
    const uint64_t num_discarded_regions = filtering_candidates_buffered->num_discarded_regions;
    filtering_candidates_buffered->discarded_regions[num_discarded_regions] = filtering_region;
    ++(filtering_candidates_buffered->num_discarded_regions);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
  }
}
void filtering_candidates_buffered_kmer_filter_retrieve(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter,
    const uint64_t gpu_buffer_kmer_filter_offset) {
  PROFILE_START(GP_FC_VERIFY_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Check filtering regions
  const uint64_t num_filtering_regions = filtering_candidates_buffered->num_regions;
  if (num_filtering_regions==0) return; // No filtering regions
  // Allocate discarded regions (buffered)
  filtering_candidates_buffered_allocate_discarded_regions(filtering_candidates_buffered,num_filtering_regions);
  filtering_candidates_buffered->num_discarded_regions = 0;
  // Traverse all filtering regions buffered
  uint64_t region_pos, gpu_buffer_offset = gpu_buffer_kmer_filter_offset;
  for (region_pos=0;region_pos<num_filtering_regions;++region_pos) {
    // Retrieve region
    filtering_region_t* const region_buffered = filtering_candidates_buffered->regions[region_pos];
    filtering_candidates_buffered_kmer_filter_retrieve_region(
        filtering_candidates_buffered,filtering_candidates,
        region_buffered,pattern,gpu_buffer_kmer_filter,&gpu_buffer_offset);
  }
  // Free (buffered regions and kmer-counting)
  filtering_candidates_buffered_free_regions(filtering_candidates_buffered);
  kmer_counting_destroy(&pattern->pattern_tiled.kmer_filter_nway,filtering_candidates->mm_allocator);
  // DEBUG
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions_BPM_buffer)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
}

