/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_bsearch.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu_buffer_fmi_bsearch.h"
#include "../resources/gpu_modules/gpu_interface.h"

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
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_pos) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search = mm_alloc(gpu_buffer_fmi_search_t);
  gpu_buffer_fmi_search->buffer = gpu_buffer_collection_get_buffer(gpu_buffer_collection,buffer_pos);
  gpu_buffer_fmi_search->num_queries = 0;
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
GEM_INLINE void gpu_buffer_fmi_search_device(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,const device_t device) {
  gpu_buffer_fmi_search->device = device;
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
      GPU_FMI_SEARCH_MAX_QUERIES,gpu_buffer_fmi_search->num_queries,max_buffer_queries);
  const uint64_t query_length = end - begin;
  gem_cond_fatal_error(query_length>=GPU_BUFFER_FMI_SEARCH_MAX_QUERY_LENGTH,
      GPU_FMI_SEARCH_MAX_QUERY_LENGTH,query_length,GPU_BUFFER_FMI_SEARCH_MAX_QUERY_LENGTH);
  // Get next free query in buffer
  gpu_fmi_search_seed_t* const gpu_fmi_search_seed =
      gpu_fmi_search_buffer_get_seeds_(gpu_buffer_fmi_search->buffer) + gpu_buffer_fmi_search->num_queries;
  // Init
  uint64_t lo = 0, hi = 0;
  // Adapt pattern chunk to query encoding
  const uint8_t* const key = pattern->key;
  uint64_t i, offset = 0;
  for (i=begin;i<end;++i) {
    lo |= (key[i] << offset);
    offset += 2;
    if (offset >= 64) break;
  }
  for (;i<end;++i) {
    hi |= (key[i] << offset);
    offset += 2;
  }
  // Adapt query length (higher 8 bits of hi)
  hi |= (query_length << 56);
  // Store & Increment queries used
  gpu_fmi_search_seed->hi = hi;
  gpu_fmi_search_seed->low = lo;
  ++(gpu_buffer_fmi_search->num_queries);
}
GEM_INLINE void gpu_buffer_fmi_search_get_result(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,
    const uint64_t buffer_pos,uint64_t* const hi,uint64_t* const lo) {
  const gpu_fmi_search_sa_inter_t* const gpu_fmi_search_sa_inter =
      gpu_fmi_search_buffer_get_sa_intervals_(gpu_buffer_fmi_search->buffer) + buffer_pos;
  *hi = gpu_fmi_search_sa_inter->hi;
  *lo = gpu_fmi_search_sa_inter->low;
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
  gpu_fmi_search_send_buffer_(gpu_buffer_fmi_search->buffer,gpu_buffer_fmi_search->num_queries);
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_SEND);
}
GEM_INLINE void gpu_buffer_fmi_search_receive(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE);
  gpu_fmi_search_receive_buffer_(gpu_buffer_fmi_search->buffer);
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
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_pos) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
GEM_INLINE void gpu_buffer_fmi_search_clear(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_search_delete(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) { GEM_CUDA_NOT_SUPPORTED(); }
GEM_INLINE void gpu_buffer_fmi_search_device(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,const device_t device) {
  GEM_CUDA_NOT_SUPPORTED();
}
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
