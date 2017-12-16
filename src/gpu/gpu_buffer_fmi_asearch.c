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
 *   GPU-adaptor module provides support data-structures and functions
 *   to offload the adaptive search of key partitions (adaptive region-profile)
 *   from a batch of searches to the GPU
 */

#include "gpu/gpu_buffer_fmi_asearch.h"
#include "fm_index/fm_index_search.h"

/*
 * Errors
 */
#define GEM_ERROR_GPU_FMI_ASEARCH_MAX_QUERIES "GPU.Adaptive.Search. Number of queries (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" queries)"
#define GEM_ERROR_GPU_FMI_ASEARCH_MAX_PATTERN_LENGTH "GPU.Adaptive.Search. Query pattern (%"PRIu64" bases) exceeds maximum query capacity (%"PRIu64" bases)"
#define GEM_ERROR_GPU_FMI_ASEARCH_MAX_REGIONS "GPU.Adaptive.Search. Query pattern (%"PRIu64" regions) exceeds maximum query capacity (%"PRIu64" regions)"

/*
 * Constants :: Buffer Hints
 */
#define GPU_ASEARCH_MIN_NUM_SAMPLES           1
#define GPU_ASEARCH_AVERAGE_QUERY_LENGTH      150
#define GPU_ASEARCH_BASES_PADDING             8

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
gpu_buffer_fmi_asearch_t* gpu_buffer_fmi_asearch_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const bool fmi_search_enabled,
    const uint32_t occ_min_threshold,
    const uint32_t extra_search_steps,
    const uint32_t alphabet_size) {
  gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch = mm_alloc(gpu_buffer_fmi_asearch_t);
  gpu_buffer_fmi_asearch->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_pos);
  gpu_buffer_fmi_asearch->fmi_search_enabled = fmi_search_enabled;
  // Init Region-Profile Parameters
  gpu_buffer_fmi_asearch->occ_min_threshold = occ_min_threshold;
  gpu_buffer_fmi_asearch->extra_search_steps = extra_search_steps;
  gpu_buffer_fmi_asearch->alphabet_size = alphabet_size;
  // Profile
  COUNTER_RESET(&gpu_buffer_fmi_asearch->query_length);
  TIMER_RESET(&gpu_buffer_fmi_asearch->timer);
  // Init gpu-buffer
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  const int64_t thread_id = gtid(); // Between [1,num_threads] (zero is master)
  gpu_alloc_buffer_(gpu_buffer_fmi_asearch->buffer, thread_id);
  gpu_fmi_asearch_init_buffer_(
      gpu_buffer_fmi_asearch->buffer,
      gpu_buffer_fmi_asearch_get_mean_query_length(gpu_buffer_fmi_asearch),
      REGION_MAX_REGIONS_FACTOR);
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  // Init Buffer state
  gpu_buffer_fmi_asearch->num_queries = 0;
  gpu_buffer_fmi_asearch->num_bases = 0;
  gpu_buffer_fmi_asearch->num_regions = 0;
  // Return
  return gpu_buffer_fmi_asearch;
}
void gpu_buffer_fmi_asearch_clear(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  // Init Buffer
  gpu_fmi_asearch_init_buffer_(
      gpu_buffer_fmi_asearch->buffer,
      gpu_buffer_fmi_asearch_get_mean_query_length(gpu_buffer_fmi_asearch),
      REGION_MAX_REGIONS_FACTOR);
  // Clear Buffer state
  gpu_buffer_fmi_asearch->num_queries = 0;
  gpu_buffer_fmi_asearch->num_bases = 0;
  gpu_buffer_fmi_asearch->num_regions = 0;
}
void gpu_buffer_fmi_asearch_delete(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  mm_free(gpu_buffer_fmi_asearch);
}
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_asearch_get_max_queries(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  return gpu_fmi_asearch_buffer_get_max_queries_(gpu_buffer_fmi_asearch->buffer);
}
uint64_t gpu_buffer_fmi_asearch_get_max_bases(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  return gpu_fmi_asearch_buffer_get_max_bases_(gpu_buffer_fmi_asearch->buffer);
}
uint64_t gpu_buffer_fmi_asearch_get_max_regions(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  return gpu_fmi_asearch_buffer_get_max_regions_(gpu_buffer_fmi_asearch->buffer);
}
uint64_t gpu_buffer_fmi_asearch_get_num_queries(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  return gpu_buffer_fmi_asearch->num_queries;
}
uint64_t gpu_buffer_fmi_asearch_get_num_bases(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  return gpu_buffer_fmi_asearch->num_bases;
}
uint64_t gpu_buffer_fmi_asearch_get_num_regions(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  return gpu_buffer_fmi_asearch->num_regions;
}
bool gpu_buffer_fmi_asearch_fits_in_buffer(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t total_queries,
    const uint64_t total_bases,
    const uint64_t total_regions) {
  // Get Limits
  uint64_t max_queries = gpu_buffer_fmi_asearch_get_max_queries(gpu_buffer_fmi_asearch);
  uint64_t max_bases = gpu_buffer_fmi_asearch_get_max_bases(gpu_buffer_fmi_asearch);
  uint64_t max_regions = gpu_buffer_fmi_asearch_get_max_regions(gpu_buffer_fmi_asearch);
  // Compute total bases + padding (worst case)
  const uint64_t total_bases_padded = total_bases + 2*GPU_ASEARCH_BASES_PADDING;
  // Check available space in buffer
  if (gpu_buffer_fmi_asearch->num_queries+total_queries > max_queries ||
      gpu_buffer_fmi_asearch->num_bases+total_bases_padded > max_bases ||
      gpu_buffer_fmi_asearch->num_regions+total_regions > max_regions) {
    // Check number of queries in the buffer
    if (gpu_buffer_fmi_asearch->num_queries > 0) {
      return false; // Leave it to the next fresh buffer
    }
    // Reallocate buffer
    gpu_fmi_asearch_init_and_realloc_buffer_(gpu_buffer_fmi_asearch->buffer,
        REGION_MAX_REGIONS_FACTOR,total_bases_padded,total_queries,total_regions);
    // Check reallocated buffer dimensions (error otherwise)
    max_queries = gpu_buffer_fmi_asearch_get_max_queries(gpu_buffer_fmi_asearch);
    gem_cond_fatal_error(total_queries > max_queries,GPU_FMI_ASEARCH_MAX_QUERIES,total_queries,max_queries);
    max_bases = gpu_buffer_fmi_asearch_get_max_bases(gpu_buffer_fmi_asearch);
    gem_cond_fatal_error(total_bases_padded > max_bases,GPU_FMI_ASEARCH_MAX_PATTERN_LENGTH,total_bases_padded,max_bases);
    max_regions = gpu_buffer_fmi_asearch_get_max_regions(gpu_buffer_fmi_asearch);
    gem_cond_fatal_error(total_regions > max_regions,GPU_FMI_ASEARCH_MAX_REGIONS,total_regions,max_regions);
    // Return OK (after reallocation)
    return true;
  }
  // Return OK (fits in buffer)
  return true;
}
uint32_t gpu_buffer_fmi_asearch_get_mean_query_length(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  if (COUNTER_GET_NUM_SAMPLES(&gpu_buffer_fmi_asearch->query_length) >= GPU_ASEARCH_MIN_NUM_SAMPLES) {
    return (uint64_t)ceil(COUNTER_GET_MEAN(&gpu_buffer_fmi_asearch->query_length));
  } else {
    return GPU_ASEARCH_AVERAGE_QUERY_LENGTH;
  }
}
/*
 * Accessors
 */
uint64_t gpu_buffer_fmi_asearch_add_query(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    pattern_t* const pattern,
    const uint64_t max_regions) {
  // Parameters
  const uint64_t key_length = pattern->key_length;
  const uint64_t query_offset = gpu_buffer_fmi_asearch->num_queries;
  gpu_fmi_search_query_info_t* const search_query_info =
      gpu_fmi_asearch_buffer_get_queries_info_(gpu_buffer_fmi_asearch->buffer) + query_offset;
  gpu_fmi_search_region_t* const search_query_region =
      gpu_fmi_asearch_buffer_get_regions_(gpu_buffer_fmi_asearch->buffer) + query_offset;
  // Add query
  search_query_info->init_offset = gpu_buffer_fmi_asearch->num_bases;
  search_query_info->query_size = key_length;
  search_query_region->init_offset = gpu_buffer_fmi_asearch->num_regions; // Set resulting regions offsets
  // Add pattern (bases)
  gpu_fmi_search_query_t* const search_query =
      gpu_fmi_asearch_buffer_get_queries_(gpu_buffer_fmi_asearch->buffer) + search_query_info->init_offset;
  const uint8_t* const key = pattern->key;
  memcpy(search_query,key,key_length); // Copy pattern
  COUNTER_ADD(&gpu_buffer_fmi_asearch->query_length,key_length); // Profile
  // Next
  ++(gpu_buffer_fmi_asearch->num_queries);
  gpu_buffer_fmi_asearch->num_regions += max_regions;
  gpu_buffer_fmi_asearch->num_bases += key_length;
  const uint64_t num_bases_mod = gpu_buffer_fmi_asearch->num_bases % GPU_ASEARCH_BASES_PADDING;
  if (num_bases_mod > 0) {
    const uint64_t bases_to_skip = GPU_ASEARCH_BASES_PADDING - num_bases_mod;
    gpu_buffer_fmi_asearch->num_bases += bases_to_skip;
  }
  // Return query-offset
  return query_offset;
}
void gpu_buffer_fmi_asearch_get_result_total_regions(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t query_offset,
    uint64_t* const regions_offset,
    uint64_t* const num_regions) {
  const gpu_fmi_search_region_t* const search_region =
      gpu_fmi_asearch_buffer_get_regions_(gpu_buffer_fmi_asearch->buffer) + query_offset;
  *regions_offset = search_region->init_offset;
  *num_regions = search_region->num_regions;
}
void gpu_buffer_fmi_asearch_get_result_region(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t region_offset,
    uint64_t* const region_begin,
    uint64_t* const region_end,
    uint64_t* const region_hi,
    uint64_t* const region_lo) {
  // Fetch region span
  gpu_fmi_search_region_info_t* const search_region_inf =
      gpu_fmi_asearch_buffer_get_regions_offsets_(gpu_buffer_fmi_asearch->buffer) + region_offset;
  *region_begin = search_region_inf->init_offset;
  *region_end = search_region_inf->end_offset;
  // Fetch region search-interval
  gpu_sa_search_inter_t* const search_inter =
      gpu_fmi_asearch_buffer_get_regions_intervals_(gpu_buffer_fmi_asearch->buffer) + region_offset;
  *region_hi = search_inter->hi;
  *region_lo = search_inter->low;
}
/*
 * Send/Receive
 */
void gpu_buffer_fmi_asearch_send(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_queries = gpu_buffer_fmi_asearch_get_max_queries(gpu_buffer_fmi_asearch);
  const uint64_t num_queries = gpu_buffer_fmi_asearch_get_num_queries(gpu_buffer_fmi_asearch);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_USAGE_QUERIES,(100*num_queries)/max_queries);
  TIMER_START(&gpu_buffer_fmi_asearch->timer);
#endif
  // Select computing device
  if (gpu_buffer_fmi_asearch->fmi_search_enabled) {
    if (gpu_buffer_fmi_asearch->num_queries > 0) {
      gpu_fmi_asearch_send_buffer_(
          gpu_buffer_fmi_asearch->buffer,gpu_buffer_fmi_asearch->num_queries,
          gpu_buffer_fmi_asearch->num_bases,gpu_buffer_fmi_asearch->num_regions,
          gpu_buffer_fmi_asearch->occ_min_threshold,gpu_buffer_fmi_asearch->extra_search_steps,
          gpu_buffer_fmi_asearch->alphabet_size);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_SEND);
}
void gpu_buffer_fmi_asearch_receive(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE);
  if (gpu_buffer_fmi_asearch->fmi_search_enabled) {
    if (gpu_buffer_fmi_asearch->num_queries > 0) {
      gpu_fmi_asearch_receive_buffer_(gpu_buffer_fmi_asearch->buffer);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE);
#ifdef GEM_PROFILE
  TIMER_STOP(&gpu_buffer_fmi_asearch->timer);
  COUNTER_ADD(&PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_DUTY_CYCLE)->time_ns,gpu_buffer_fmi_asearch->timer.accumulated);
#endif
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
gpu_buffer_fmi_asearch_t* gpu_buffer_fmi_asearch_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const bool fmi_search_enabled,
    const uint32_t occ_min_threshold,
    const uint32_t extra_search_steps,
    const uint32_t alphabet_size) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
void gpu_buffer_fmi_asearch_clear(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_asearch_delete(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_asearch_get_max_queries(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_fmi_asearch_get_max_bases(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_fmi_asearch_get_max_regions(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_fmi_asearch_get_num_queries(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_fmi_asearch_get_num_bases(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_fmi_asearch_get_num_regions(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
bool gpu_buffer_fmi_asearch_fits_in_buffer(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t total_queries,
    const uint64_t total_bases,
    const uint64_t total_regions) { GEM_CUDA_NOT_SUPPORTED(); return false; }
uint32_t gpu_buffer_fmi_asearch_get_mean_query_length(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
/*
 * Accessors
 */
uint64_t gpu_buffer_fmi_asearch_add_query(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    pattern_t* const pattern,
    const uint64_t max_regions) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
void gpu_buffer_fmi_asearch_get_result_total_regions(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t query_offset,
    uint64_t* const regions_offset,
    uint64_t* const num_regions) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_asearch_get_result_region(
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch,
    const uint64_t region_offset,
    uint64_t* const region_begin,
    uint64_t* const region_end,
    uint64_t* const region_hi,
    uint64_t* const region_lo) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Send/Receive
 */
void gpu_buffer_fmi_asearch_send(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_asearch_receive(gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) { GEM_CUDA_NOT_SUPPORTED(); }
#endif /* HAVE_CUDA */
