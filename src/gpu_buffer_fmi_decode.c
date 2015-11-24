/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_decode.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu_buffer_fmi_decode.h"
#include "../resources/gpu_modules/gpu_interface.h"

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
GEM_INLINE gpu_buffer_fmi_decode_t* gpu_buffer_fmi_decode_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,fm_index_t* const fm_index) {
  PROF_START(GP_GPU_BUFFER_FMI_DECODE_ALLOC);
  gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode = mm_alloc(gpu_buffer_fmi_decode_t);
  gpu_buffer_fmi_decode->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_pos);
  gpu_buffer_fmi_decode->num_queries = 0;
  gpu_buffer_fmi_decode->compute_cpu = false;
  gpu_buffer_fmi_decode->fm_index = fm_index;
  // Init Buffer
  gpu_alloc_buffer_(gpu_buffer_fmi_decode->buffer);
  gpu_fmi_decode_init_buffer_(gpu_buffer_fmi_decode->buffer);
  PROF_STOP(GP_GPU_BUFFER_FMI_DECODE_ALLOC);
  // Return
  return gpu_buffer_fmi_decode;
}
GEM_INLINE void gpu_buffer_fmi_decode_clear(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  // Init Buffer
  gpu_fmi_decode_init_buffer_(gpu_buffer_fmi_decode->buffer); // TODO-ALEJANDRO Always?
  // Clear
  gpu_buffer_fmi_decode->num_queries = 0;
}
GEM_INLINE void gpu_buffer_fmi_decode_delete(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  mm_free(gpu_buffer_fmi_decode);
}
/*
 * Computing Device
 */
GEM_INLINE void gpu_buffer_fmi_decode_set_device_cpu(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  gpu_buffer_fmi_decode->compute_cpu = true;
}
GEM_INLINE void gpu_buffer_fmi_decode_set_device_gpu(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  gpu_buffer_fmi_decode->compute_cpu = false;
}
/*
 * Occupancy & Limits
 */
GEM_INLINE uint64_t gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  return gpu_fmi_decode_buffer_get_max_positions_(gpu_buffer_fmi_decode->buffer);
}
GEM_INLINE uint64_t gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  return gpu_buffer_fmi_decode->num_queries;
}
/*
 * Accessors
 */
GEM_INLINE void gpu_buffer_fmi_decode_add_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t bwt_position) {
  PROF_INC_COUNTER(GP_GPU_BUFFER_FMI_DECODE_NUM_QUERIES);
  // Check query length & buffer occupancy
  const uint64_t max_buffer_queries = gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode);
  gem_cond_fatal_error(gpu_buffer_fmi_decode->num_queries >= max_buffer_queries,
      GPU_FMI_DECODE_MAX_QUERIES,(uint64_t)gpu_buffer_fmi_decode->num_queries,max_buffer_queries);
  // Get next free query in buffer
  gpu_fmi_decode_init_pos_t* const gpu_fmi_decode_init_pos =
      gpu_fmi_decode_buffer_get_init_pos_(gpu_buffer_fmi_decode->buffer) + gpu_buffer_fmi_decode->num_queries;
  *gpu_fmi_decode_init_pos = bwt_position; // Set decode position
  // Increment number of queries
  ++(gpu_buffer_fmi_decode->num_queries);
}
GEM_INLINE void gpu_buffer_fmi_decode_get_result(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t buffer_pos,
    uint64_t* const bwt_sampled_position,uint64_t* const lf_steps) {
  // Get query
  const gpu_fmi_decode_end_pos_t* const gpu_fmi_decode_end_pos =
      gpu_fmi_decode_buffer_get_end_pos_(gpu_buffer_fmi_decode->buffer) + buffer_pos;
  *bwt_sampled_position = gpu_fmi_decode_end_pos->interval;
  *lf_steps = gpu_fmi_decode_end_pos->steps;
}
/*
 * CPU emulated
 */
GEM_INLINE void gpu_buffer_fmi_decode_compute_cpu(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
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
GEM_INLINE void gpu_buffer_fmi_decode_send(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROF_START(GP_GPU_BUFFER_FMI_DECODE_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_queries = gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode);
  const uint64_t used_queries = gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_FMI_DECODE_USAGE_CANDIDATES,(100*used_queries)/max_queries);
#endif
  // Select computing device
  if (!gpu_buffer_fmi_decode->compute_cpu) {
    sampled_sa_t* const sampled_sa = gpu_buffer_fmi_decode->fm_index->sampled_sa;
    const uint32_t sampling_rate = sampled_sa_get_sa_sampling_rate(sampled_sa);
    gpu_fmi_decode_send_buffer_(gpu_buffer_fmi_decode->buffer,gpu_buffer_fmi_decode->num_queries,sampling_rate);
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_DECODE_SEND);
}
GEM_INLINE void gpu_buffer_fmi_decode_receive(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROF_START(GP_GPU_BUFFER_FMI_DECODE_RECEIVE);
  // Select computing device
  if (!gpu_buffer_fmi_decode->compute_cpu) {
    gpu_fmi_decode_receive_buffer_(gpu_buffer_fmi_decode->buffer);
  } else {
    // CPU emulated
    gpu_buffer_fmi_decode_compute_cpu(gpu_buffer_fmi_decode);
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_DECODE_RECEIVE);
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
GEM_INLINE gpu_buffer_fmi_decode_t* gpu_buffer_fmi_decode_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,fm_index_t* const fm_index) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
GEM_INLINE void gpu_buffer_fmi_decode_clear(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_decode_delete(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Computing Device
 */
GEM_INLINE void gpu_buffer_fmi_decode_set_device_cpu(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_decode_set_device_gpu(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
GEM_INLINE uint64_t gpu_buffer_fmi_decode_get_max_queries(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
GEM_INLINE uint64_t gpu_buffer_fmi_decode_get_num_queries(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
/*
 * Accessors
 */
GEM_INLINE void gpu_buffer_fmi_decode_add_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t bwt_position) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_decode_get_result(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t buffer_pos,
    uint64_t* const bwt_sampled_position,uint64_t* const lf_steps) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Send/Receive
 */
GEM_INLINE void gpu_buffer_fmi_decode_send(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_decode_receive(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) { GEM_CUDA_NOT_SUPPORTED(); }
#endif /* HAVE_CUDA */
