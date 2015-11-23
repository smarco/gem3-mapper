/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_bsearch.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu_buffer_fmi_bsearch.h"
#include "../resources/gpu_modules/gpu_interface.h"
#include "fm_index_search.h"

/*
 * Constants
 */
#define GPU_BUFFER_FMI_SEARCH_MAX_QUERY_LENGTH (((64*2)-8)/2)  /* 60 bases */

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
GEM_INLINE gpu_buffer_fmi_search_t* gpu_buffer_fmi_search_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,fm_index_t* const fm_index) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search = mm_alloc(gpu_buffer_fmi_search_t);
  gpu_buffer_fmi_search->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_pos);
  gpu_buffer_fmi_search->num_queries = 0;
  gpu_buffer_fmi_search->compute_cpu = false;
  gpu_buffer_fmi_search->fm_index = fm_index;
  // Init Buffer
  gpu_alloc_buffer_(gpu_buffer_fmi_search->buffer);
  gpu_fmi_search_init_buffer_(gpu_buffer_fmi_search->buffer);
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  // Return
  return gpu_buffer_fmi_search;
}
GEM_INLINE void gpu_buffer_fmi_search_clear(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  gpu_fmi_search_init_buffer_(gpu_buffer_fmi_search->buffer); // Init Buffer
  gpu_buffer_fmi_search->num_queries = 0; // Clear
}
GEM_INLINE void gpu_buffer_fmi_search_delete(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  mm_free(gpu_buffer_fmi_search);
}
/*
 * Computing Device
 */
GEM_INLINE void gpu_buffer_fmi_search_set_device_cpu(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  gpu_buffer_fmi_search->compute_cpu = true;
}
GEM_INLINE void gpu_buffer_fmi_search_set_device_gpu(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  gpu_buffer_fmi_search->compute_cpu = false;
}
/*
 * Occupancy & Limits
 */
GEM_INLINE uint64_t gpu_buffer_fmi_search_get_max_queries(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  return gpu_fmi_search_buffer_get_max_seeds_(gpu_buffer_fmi_search->buffer);
}
GEM_INLINE uint64_t gpu_buffer_fmi_search_get_num_queries(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  return gpu_buffer_fmi_search->num_queries;
}
/*
 * Accessors
 */
GEM_INLINE void gpu_buffer_fmi_search_add_query(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,
    pattern_t* const pattern,const uint64_t begin,const uint64_t end) {
  PROF_INC_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_NUM_QUERIES);
  // Check query length & buffer occupancy
  const uint64_t max_buffer_queries = gpu_buffer_fmi_search_get_max_queries(gpu_buffer_fmi_search);
  gem_cond_fatal_error(gpu_buffer_fmi_search->num_queries >= max_buffer_queries,
      GPU_FMI_SEARCH_MAX_QUERIES,(uint64_t)gpu_buffer_fmi_search->num_queries,max_buffer_queries);
  const uint64_t query_length = end - begin;
  gem_cond_fatal_error(query_length>=GPU_BUFFER_FMI_SEARCH_MAX_QUERY_LENGTH,
      GPU_FMI_SEARCH_MAX_QUERY_LENGTH,query_length,(uint64_t)GPU_BUFFER_FMI_SEARCH_MAX_QUERY_LENGTH);
  // Get next free query in buffer
  gpu_fmi_search_seed_t* const gpu_fmi_search_seed =
      gpu_fmi_search_buffer_get_seeds_(gpu_buffer_fmi_search->buffer) + gpu_buffer_fmi_search->num_queries;
  // Select computing device
  if (gpu_buffer_fmi_search->compute_cpu) {
    gpu_fmi_search_seed->hi = (uint64_t)pattern; // Yes, I know what I'm doing here
    gpu_fmi_search_seed->low = (begin << 32) | end;
  } else {
    // Init
    uint64_t lo = 0, hi = 0, offset = 0;
    // Adapt pattern chunk to query encoding
    const uint8_t* const key = pattern->key;
    int64_t i;
    for (i=end-1;i>=begin;--i) {
      lo |= (key[i] << offset);
      offset += 2;
      if (offset >= 64) break;
    }
    for (offset=0;i>=begin;--i) {
      hi |= (key[i] << offset);
      offset += 2;
    }
    // Adapt query length (higher 8 bits of hi)
    hi |= (query_length << 56);
    // Store & Increment queries used
    gpu_fmi_search_seed->hi = hi;
    gpu_fmi_search_seed->low = lo;
  }
  //
  ++(gpu_buffer_fmi_search->num_queries);
}
GEM_INLINE void gpu_buffer_fmi_search_get_result(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,
    const uint64_t buffer_pos,uint64_t* const hi,uint64_t* const lo) {
  const gpu_fmi_search_sa_inter_t* const gpu_fmi_search_sa_interval =
      gpu_fmi_search_buffer_get_sa_intervals_(gpu_buffer_fmi_search->buffer) + buffer_pos;
  *hi = gpu_fmi_search_sa_interval->hi;
  *lo = gpu_fmi_search_sa_interval->low;
}
/*
 * CPU emulated
 */
GEM_INLINE void gpu_buffer_fmi_search_compute_cpu(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  // Parameters
  gpu_fmi_search_seed_t* const gpu_fmi_search_seeds =
      gpu_fmi_search_buffer_get_seeds_(gpu_buffer_fmi_search->buffer);
  gpu_fmi_search_sa_inter_t* const gpu_fmi_search_sa_intervals =
      gpu_fmi_search_buffer_get_sa_intervals_(gpu_buffer_fmi_search->buffer);
  const uint64_t num_queries = gpu_buffer_fmi_search->num_queries;
  // Traverse all queries
  uint64_t buffer_pos;
  for (buffer_pos=0;buffer_pos<num_queries;++buffer_pos) {
    // Get Query
    gpu_fmi_search_seed_t* const gpu_fmi_search_seed = gpu_fmi_search_seeds + buffer_pos;
    pattern_t* const pattern = (pattern_t*)gpu_fmi_search_seed->hi;
    const uint64_t begin = gpu_fmi_search_seed->low >> 32;
    const uint64_t end = gpu_fmi_search_seed->low & 0x00000000FFFFFFFF;
    // Search Query
    uint64_t hi, lo;
    fm_index_bsearch(gpu_buffer_fmi_search->fm_index,pattern->key+begin,end-begin,&hi,&lo);
    // Set Result
    gpu_fmi_search_sa_inter_t* const gpu_fmi_search_sa_interval = gpu_fmi_search_sa_intervals + buffer_pos;
    gpu_fmi_search_sa_interval->low = lo;
    gpu_fmi_search_sa_interval->hi = hi;
  }
}
/*
 * Send/Receive
 */
GEM_INLINE void gpu_buffer_fmi_search_send(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_SEND);
#ifdef GEM_PROFILE
  const uint64_t max_queries = gpu_buffer_fmi_search_get_max_queries(gpu_buffer_fmi_search);
  const uint64_t used_queries = gpu_buffer_fmi_search_get_num_queries(gpu_buffer_fmi_search);
  PROF_ADD_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_USAGE_CANDIDATES,(100*used_queries)/max_queries);
#endif
  // Select computing device
  if (!gpu_buffer_fmi_search->compute_cpu) {
    gpu_fmi_search_send_buffer_(gpu_buffer_fmi_search->buffer,gpu_buffer_fmi_search->num_queries);
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_SEND);
}
GEM_INLINE void gpu_buffer_fmi_search_receive(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE);
  // Select computing device
  if (!gpu_buffer_fmi_search->compute_cpu) {
    gpu_fmi_search_receive_buffer_(gpu_buffer_fmi_search->buffer);
  } else {
    // CPU emulated
    gpu_buffer_fmi_search_compute_cpu(gpu_buffer_fmi_search);
  }
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE);
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
GEM_INLINE gpu_buffer_fmi_search_t* gpu_buffer_fmi_search_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,fm_index_t* const fm_index) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
GEM_INLINE void gpu_buffer_fmi_search_clear(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_search_delete(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Computing Device
 */
GEM_INLINE void gpu_buffer_fmi_search_set_device_cpu(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_search_set_device_gpu(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Occupancy & Limits
 */
GEM_INLINE uint64_t gpu_buffer_fmi_search_get_max_queries(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
GEM_INLINE uint64_t gpu_buffer_fmi_search_get_num_queries(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); return 0; }
/*
 * Accessors
 */
GEM_INLINE void gpu_buffer_fmi_search_add_query(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,
    pattern_t* const pattern,const uint64_t begin,const uint64_t end) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_search_get_result(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,
    const uint64_t buffer_pos,uint64_t* const hi,uint64_t* const lo) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Send/Receive
 */
GEM_INLINE void gpu_buffer_fmi_search_send(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_search_receive(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); }
#endif /* HAVE_CUDA */
