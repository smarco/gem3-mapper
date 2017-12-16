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
 *   to offload the static search of key partitions (fixed region-profile)
 *   from a batch of searches to the GPU
 */

#include "gpu/gpu_buffer_fmi_ssearch.h"
#include "fm_index/fm_index_search.h"

/*
 * Constants
 */
#define GPU_BUFFER_FMI_SEARCH_MAX_QUERY_LENGTH (((64*2)-8)/2)  /* 60 bases */

/*
 * Errors
 */
#define GEM_ERROR_GPU_FMI_SEARCH_MAX_QUERY_LENGTH "GPU.FMI.Search. Query pattern (%"PRIu64" bases) exceeds maximum query capacity (%"PRIu64" bases)"
#define GEM_ERROR_GPU_FMI_SEARCH_MAX_QUERIES "GPU.FMI.Search. Number of queries (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" queries)"

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
gpu_buffer_fmi_ssearch_t* gpu_buffer_fmi_ssearch_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const bool fmi_search_enabled) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch = mm_alloc(gpu_buffer_fmi_ssearch_t);
  gpu_buffer_fmi_ssearch->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_pos);
  gpu_buffer_fmi_ssearch->num_queries = 0;
  gpu_buffer_fmi_ssearch->fmi_search_enabled = fmi_search_enabled;
  TIMER_RESET(&gpu_buffer_fmi_ssearch->timer);
  // Init Buffer
  const int64_t thread_id = gtid(); // Between [1,num_threads] (zero is master)
  gpu_alloc_buffer_(gpu_buffer_fmi_ssearch->buffer, thread_id);
  gpu_fmi_ssearch_init_buffer_(gpu_buffer_fmi_ssearch->buffer);
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  // Return
  return gpu_buffer_fmi_ssearch;
}
void gpu_buffer_fmi_ssearch_clear(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  gpu_fmi_ssearch_init_buffer_(gpu_buffer_fmi_ssearch->buffer); // Init Buffer
  gpu_buffer_fmi_ssearch->num_queries = 0; // Clear
}
void gpu_buffer_fmi_ssearch_delete(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  mm_free(gpu_buffer_fmi_ssearch);
}
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_ssearch_get_max_queries(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  return gpu_fmi_ssearch_buffer_get_max_seeds_(gpu_buffer_fmi_ssearch->buffer);
}
uint64_t gpu_buffer_fmi_ssearch_get_num_queries(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  return gpu_buffer_fmi_ssearch->num_queries;
}
bool gpu_buffer_fmi_ssearch_fits_in_buffer(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    const uint64_t total_queries) {
  // Get Limits
  uint64_t max_queries = gpu_buffer_fmi_ssearch_get_max_queries(gpu_buffer_fmi_ssearch);
  // Check available space in buffer
  if (gpu_buffer_fmi_ssearch->num_queries+total_queries > max_queries) {
    // Check number of queries in the buffer
    if (gpu_buffer_fmi_ssearch->num_queries > 0) {
      return false; // Leave it to the next fresh buffer
    }
    // Reallocate buffer
    gpu_fmi_ssearch_init_and_realloc_buffer_(gpu_buffer_fmi_ssearch->buffer,total_queries);
    // Check reallocated buffer dimensions (error otherwise)
    max_queries = gpu_buffer_fmi_ssearch_get_max_queries(gpu_buffer_fmi_ssearch);
    gem_cond_fatal_error(total_queries > max_queries,GPU_FMI_SEARCH_MAX_QUERIES,total_queries,max_queries);
    // Return OK (after reallocation)
    return true;
  }
  // Return OK (fits in buffer)
  return true;
}
/*
 * Accessors
 */
void gpu_buffer_fmi_ssearch_add_query(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    pattern_t* const pattern,
    const uint64_t begin,
    const uint64_t end) {
  PROF_INC_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_NUM_QUERIES);
  // Check enabled
  if (!gpu_buffer_fmi_ssearch->fmi_search_enabled) return;
  // Check query length
  const uint64_t query_length = end - begin;
  gem_cond_fatal_error(query_length>=GPU_BUFFER_FMI_SEARCH_MAX_QUERY_LENGTH,
      GPU_FMI_SEARCH_MAX_QUERY_LENGTH,query_length,(uint64_t)GPU_BUFFER_FMI_SEARCH_MAX_QUERY_LENGTH);
  // Get next free query in buffer
  gpu_fmi_search_seed_t* const gpu_fmi_ssearch_seed =
      gpu_fmi_ssearch_buffer_get_seeds_(gpu_buffer_fmi_ssearch->buffer) + gpu_buffer_fmi_ssearch->num_queries;
  // Adapt pattern chunk to query encoding
  const int64_t begin_region = begin;
  const int64_t end_region = end-1;
  const uint8_t* const key = pattern->key;
  uint64_t lo = 0, hi = 0, offset = 0; // Init
  int64_t i;
  for (i=end_region;i>=begin_region;--i) {
    hi |= (((uint64_t)key[i]) << offset);
    offset += 2;
    if (offset >= 64) break;
  }
  for (offset=0;i>=begin_region;--i) {
    lo |= (((uint64_t)key[i]) << offset);
    offset += 2;
  }
  // Adapt query length (higher 8 bits of @lo)
  lo |= (query_length << 56);
  // Store & Increment queries used
  gpu_fmi_ssearch_seed->hi = hi;
  gpu_fmi_ssearch_seed->low = lo;
  // Next query
  ++(gpu_buffer_fmi_ssearch->num_queries);
}
void gpu_buffer_fmi_ssearch_get_result(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    const uint64_t buffer_pos,
    uint64_t* const hi,
    uint64_t* const lo) {
  const gpu_sa_search_inter_t* const gpu_fmi_ssearch_sa_interval =
      gpu_fmi_ssearch_buffer_get_sa_intervals_(gpu_buffer_fmi_ssearch->buffer) + buffer_pos;
  *hi = gpu_fmi_ssearch_sa_interval->hi;
  *lo = gpu_fmi_ssearch_sa_interval->low;
}
/*
 * Send/Receive
 */
void gpu_buffer_fmi_ssearch_send(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_queries = gpu_buffer_fmi_ssearch_get_max_queries(gpu_buffer_fmi_ssearch);
  const uint64_t used_queries = gpu_buffer_fmi_ssearch_get_num_queries(gpu_buffer_fmi_ssearch);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_USAGE_QUERIES,(100*used_queries)/max_queries);
  TIMER_START(&gpu_buffer_fmi_ssearch->timer);
#endif
  // Select computing device
  if (gpu_buffer_fmi_ssearch->fmi_search_enabled) {
    if (gpu_buffer_fmi_ssearch->num_queries > 0) {
      gpu_fmi_ssearch_send_buffer_(gpu_buffer_fmi_ssearch->buffer,gpu_buffer_fmi_ssearch->num_queries);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_SEND);
}
void gpu_buffer_fmi_ssearch_receive(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE);
  if (gpu_buffer_fmi_ssearch->fmi_search_enabled) {
    if (gpu_buffer_fmi_ssearch->num_queries > 0) {
      gpu_fmi_ssearch_receive_buffer_(gpu_buffer_fmi_ssearch->buffer);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE);
#ifdef GEM_PROFILE
  TIMER_STOP(&gpu_buffer_fmi_ssearch->timer);
  COUNTER_ADD(&PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_DUTY_CYCLE)->time_ns,gpu_buffer_fmi_ssearch->timer.accumulated);
#endif
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
gpu_buffer_fmi_ssearch_t* gpu_buffer_fmi_ssearch_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const bool fmi_search_enabled) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
void gpu_buffer_fmi_ssearch_clear(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_ssearch_delete(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_ssearch_get_max_queries(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_fmi_ssearch_get_num_queries(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
bool gpu_buffer_fmi_ssearch_fits_in_buffer(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    const uint64_t total_queries) { GEM_CUDA_NOT_SUPPORTED(); return false; }
/*
 * Accessors
 */
void gpu_buffer_fmi_ssearch_add_query(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    pattern_t* const pattern,
    const uint64_t begin,
    const uint64_t end) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_ssearch_get_result(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    const uint64_t buffer_pos,
    uint64_t* const hi,
    uint64_t* const lo) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Send/Receive
 */
void gpu_buffer_fmi_ssearch_send(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_ssearch_receive(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) { GEM_CUDA_NOT_SUPPORTED(); }

#endif /* HAVE_CUDA */
