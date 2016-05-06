/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_decode.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu/gpu_buffer_fmi_decode.h"
#include "resources/gpu_modules/gpu_interface.h"


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
    fm_index_t* const fm_index,
    const bool gpu_decode_sa,
    const bool gpu_decode_text) {
  PROF_START(GP_GPU_BUFFER_FMI_DECODE_ALLOC);
  gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode = mm_alloc(gpu_buffer_fmi_decode_t);
  gpu_buffer_fmi_decode->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_pos);
  gpu_buffer_fmi_decode->num_queries = 0;
  gpu_buffer_fmi_decode->gpu_decode_sa = gpu_decode_sa;
  gpu_buffer_fmi_decode->gpu_decode_text = gpu_decode_text;
  gpu_buffer_fmi_decode->fm_index = fm_index;
  TIMER_RESET(&gpu_buffer_fmi_decode->timer);
  // Init Buffer
  gpu_alloc_buffer_(gpu_buffer_fmi_decode->buffer);
  gpu_fmi_decode_init_buffer_(gpu_buffer_fmi_decode->buffer);
  PROF_STOP(GP_GPU_BUFFER_FMI_DECODE_ALLOC);
  // Return
  return gpu_buffer_fmi_decode;
}
void gpu_buffer_fmi_decode_clear(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  // Init Buffer
  gpu_fmi_decode_init_buffer_(gpu_buffer_fmi_decode->buffer); // TODO-ALEJANDRO Always?
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
 * CPU emulated
 */
void gpu_buffer_fmi_decode_compute_cpu(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  // Parameters
  gpu_fmi_decode_init_pos_t* const gpu_fmi_decode_bwt_positions =
      gpu_fmi_decode_buffer_get_init_pos_(gpu_buffer_fmi_decode->buffer);
  gpu_fmi_decode_end_pos_t* const gpu_fmi_decode_sampled_bwt_positions =
      gpu_fmi_decode_buffer_get_end_pos_(gpu_buffer_fmi_decode->buffer);
  const uint64_t num_queries = gpu_buffer_fmi_decode->num_queries;
  // Traverse all queries
  uint64_t buffer_pos;
  for (buffer_pos=0;buffer_pos<num_queries;++buffer_pos) {
    // Get Query/Result
    gpu_fmi_decode_init_pos_t* const gpu_fmi_decode_bwt_position = gpu_fmi_decode_bwt_positions + buffer_pos;
    gpu_fmi_decode_end_pos_t* const gpu_fmi_decode_sampled_bwt_position = gpu_fmi_decode_sampled_bwt_positions + buffer_pos;
    // Retrieve sampled BWT-position
    fm_index_retrieve_bwt_sampled(gpu_buffer_fmi_decode->fm_index,*gpu_fmi_decode_bwt_position,
        &gpu_fmi_decode_sampled_bwt_position->interval,&gpu_fmi_decode_sampled_bwt_position->steps);
  }
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
  // Select computing device
  if (gpu_buffer_fmi_decode->gpu_decode_sa) {
    if (gpu_buffer_fmi_decode->num_queries > 0) {
      sampled_sa_t* const sampled_sa = gpu_buffer_fmi_decode->fm_index->sampled_sa;
      const uint32_t sampling_rate = sampled_sa_get_sa_sampling_rate(sampled_sa);
      gpu_fmi_decode_send_buffer_(gpu_buffer_fmi_decode->buffer,gpu_buffer_fmi_decode->num_queries,sampling_rate);
    }
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_DECODE_SEND);
}
void gpu_buffer_fmi_decode_receive(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROF_START(GP_GPU_BUFFER_FMI_DECODE_RECEIVE);
  if (gpu_buffer_fmi_decode->num_queries > 0) {
    // Select computing device
    if (gpu_buffer_fmi_decode->gpu_decode_sa) {
      gpu_fmi_decode_receive_buffer_(gpu_buffer_fmi_decode->buffer);
    } else {
      gpu_buffer_fmi_decode_compute_cpu(gpu_buffer_fmi_decode); // CPU emulated
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
    fm_index_t* const fm_index,
    const bool gpu_decode_sa,
    const bool gpu_decode_text) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
void gpu_buffer_fmi_decode_clear(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_decode_delete(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Computing Device
 */
void gpu_buffer_fmi_decode_set_device_cpu(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
void gpu_buffer_fmi_decode_set_device_gpu(
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
