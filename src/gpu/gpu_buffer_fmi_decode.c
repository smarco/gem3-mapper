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
 *   to offload the decoding of bwt-positions to text-positions
 *   from a batch of searches to the GPU
 */

#include "gpu/gpu_buffer_fmi_decode.h"

/*
 * Errors
 */
#define GEM_ERROR_GPU_FMI_DECODE_MAX_QUERIES "GPU.FMI.Decode. Number of queries (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" decodes)"

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
gpu_buffer_fmi_decode_t* gpu_buffer_fmi_decode_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const uint32_t sampling_rate,
    const bool decode_sa_enabled,
    const bool decode_text_enabled) {
  PROF_START(GP_GPU_BUFFER_FMI_DECODE_ALLOC);
  gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode = mm_alloc(gpu_buffer_fmi_decode_t);
  gpu_buffer_fmi_decode->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_pos);
  gpu_buffer_fmi_decode->num_queries = 0;
  gpu_buffer_fmi_decode->sampling_rate = sampling_rate;
  gpu_buffer_fmi_decode->decode_sa_enabled = decode_sa_enabled;
  gpu_buffer_fmi_decode->decode_text_enabled = decode_text_enabled;
  TIMER_RESET(&gpu_buffer_fmi_decode->timer);
  // Init Buffer
  const int64_t thread_id = gtid(); // Between [1,num_threads] (zero is master)
  gpu_alloc_buffer_(gpu_buffer_fmi_decode->buffer, thread_id);
  gpu_fmi_decode_init_buffer_(gpu_buffer_fmi_decode->buffer);
  PROF_STOP(GP_GPU_BUFFER_FMI_DECODE_ALLOC);
  // Return
  return gpu_buffer_fmi_decode;
}
void gpu_buffer_fmi_decode_clear(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  // Init Buffer
  gpu_fmi_decode_init_buffer_(gpu_buffer_fmi_decode->buffer);
  // Clear
  gpu_buffer_fmi_decode->num_queries = 0;
}
void gpu_buffer_fmi_decode_delete(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  mm_free(gpu_buffer_fmi_decode);
}
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  return gpu_fmi_decode_buffer_get_max_positions_(gpu_buffer_fmi_decode->buffer);
}
uint64_t gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  return gpu_buffer_fmi_decode->num_queries;
}
bool gpu_buffer_fmi_decode_fits_in_buffer(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t total_queries) {
  // Get Limits
  uint64_t max_queries = gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode);
  // Check available space in buffer
  if (gpu_buffer_fmi_decode->num_queries+total_queries > max_queries) {
    // Check number of queries in the buffer
    if (gpu_buffer_fmi_decode->num_queries > 0) {
      return false; // Leave it to the next fresh buffer
    }
    // Reallocate buffer
    gpu_fmi_decode_init_and_realloc_buffer_(gpu_buffer_fmi_decode->buffer,total_queries);
    // Check reallocated buffer dimensions (error otherwise)
    max_queries = gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode);
    gem_cond_fatal_error(total_queries > max_queries,GPU_FMI_DECODE_MAX_QUERIES,total_queries,max_queries);
    // Return OK (after reallocation)
    return true;
  }
  // Return OK (fits in buffer)
  return true;
}
/*
 * Accessors
 */
void gpu_buffer_fmi_decode_add_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t bwt_position) {
  PROF_INC_COUNTER(GP_GPU_BUFFER_FMI_DECODE_NUM_QUERIES);
  // Get next free query in buffer
  gpu_fmi_decode_init_pos_t* const gpu_fmi_decode_init_pos =
      gpu_fmi_decode_buffer_get_init_pos_(gpu_buffer_fmi_decode->buffer) + gpu_buffer_fmi_decode->num_queries;
  *gpu_fmi_decode_init_pos = bwt_position; // Set decode position
  // Increment number of queries
  ++(gpu_buffer_fmi_decode->num_queries);
}
void gpu_buffer_fmi_decode_get_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const position) {
  // Get query
  gpu_fmi_decode_init_pos_t* const gpu_fmi_decode_init_pos =
      gpu_fmi_decode_buffer_get_init_pos_(gpu_buffer_fmi_decode->buffer) + buffer_pos;
  *position = *gpu_fmi_decode_init_pos;
}
void gpu_buffer_fmi_decode_get_position_sa(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const bwt_sampled_position,
    uint64_t* const lf_steps) {
  // Get query
  const gpu_fmi_decode_end_pos_t* const gpu_fmi_decode_end_pos =
      gpu_fmi_decode_buffer_get_end_pos_(gpu_buffer_fmi_decode->buffer) + buffer_pos;
  *bwt_sampled_position = gpu_fmi_decode_end_pos->interval;
  *lf_steps = gpu_fmi_decode_end_pos->steps;
}
void gpu_buffer_fmi_decode_get_position_text(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const text_position) {
  // Get query
  const gpu_sa_decode_text_pos_t* const gpu_fmi_decode_end_pos =
      gpu_sa_decode_buffer_get_ref_pos_(gpu_buffer_fmi_decode->buffer) + buffer_pos;
  *text_position = *gpu_fmi_decode_end_pos;
}
/*
 * Send/Receive
 */
void gpu_buffer_fmi_decode_send(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROF_START(GP_GPU_BUFFER_FMI_DECODE_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_queries = gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode);
  const uint64_t used_queries = gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_FMI_DECODE_USAGE_CANDIDATES,(100*used_queries)/max_queries);
  TIMER_START(&gpu_buffer_fmi_decode->timer);
#endif
  if (gpu_buffer_fmi_decode->num_queries > 0) {
    if (gpu_buffer_fmi_decode->decode_sa_enabled || gpu_buffer_fmi_decode->decode_text_enabled) {
      gpu_fmi_decode_send_buffer_(gpu_buffer_fmi_decode->buffer,
          gpu_buffer_fmi_decode->num_queries,gpu_buffer_fmi_decode->sampling_rate);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_DECODE_SEND);
}
void gpu_buffer_fmi_decode_receive(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROF_START(GP_GPU_BUFFER_FMI_DECODE_RECEIVE);
  if (gpu_buffer_fmi_decode->num_queries > 0) {
    if (gpu_buffer_fmi_decode->decode_sa_enabled || gpu_buffer_fmi_decode->decode_text_enabled) {
      gpu_fmi_decode_receive_buffer_(gpu_buffer_fmi_decode->buffer);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_DECODE_RECEIVE);
#ifdef GEM_PROFILE
  TIMER_STOP(&gpu_buffer_fmi_decode->timer);
  COUNTER_ADD(&PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_DUTY_CYCLE)->time_ns,gpu_buffer_fmi_decode->timer.accumulated);
#endif
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
gpu_buffer_fmi_decode_t* gpu_buffer_fmi_decode_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const uint32_t sampling_rate,
    const bool decode_sa_enabled,
    const bool decode_text_enabled) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
void gpu_buffer_fmi_decode_clear(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_decode_delete(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_decode_get_max_queries(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
uint64_t gpu_buffer_fmi_decode_get_num_queries(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
bool gpu_buffer_fmi_decode_fits_in_buffer(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t total_queries) { GEM_CUDA_NOT_SUPPORTED(); return false; }
/*
 * Accessors
 */
void gpu_buffer_fmi_decode_add_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t bwt_position) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_decode_get_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const position) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_decode_get_position_sa(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const bwt_sampled_position,
    uint64_t* const lf_steps) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_decode_get_position_text(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const text_position) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Send/Receive
 */
void gpu_buffer_fmi_decode_send(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_decode_receive(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
#endif /* HAVE_CUDA */
