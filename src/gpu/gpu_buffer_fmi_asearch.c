/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_asearch.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu/gpu_buffer_fmi_asearch.h"
#include "fm_index/fm_index_search.h"
#include "resources/gpu_modules/gpu_interface.h"

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
  // Init gpu-buffer
  PROF_START(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  gpu_alloc_buffer_(gpu_buffer_fmi_asearch->buffer);
  gpu_fmi_asearch_init_buffer_(
      gpu_buffer_fmi_asearch->buffer,
      gpu_buffer_fmi_asearch_get_mean_query_length(gpu_buffer_fmi_asearch),
      REGION_MAX_REGIONS_FACTOR);
  PROF_STOP(GP_GPU_BUFFER_FMI_SEARCH_ALLOC);
  // Init Buffer state
  gpu_buffer_fmi_asearch->num_queries = 0;
  gpu_buffer_fmi_asearch->num_bases = 0;
  gpu_buffer_fmi_asearch->num_regions = 0;
  // Profile
  COUNTER_RESET(&gpu_buffer_fmi_asearch->query_length);
  TIMER_RESET(&gpu_buffer_fmi_asearch->timer);
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
  // Check available space in buffer
  if (gpu_buffer_fmi_asearch->num_queries+total_queries > max_queries ||
      gpu_buffer_fmi_asearch->num_bases+total_bases > max_bases ||
      gpu_buffer_fmi_asearch->num_regions+total_regions > max_regions) {
    // Check number of queries in the buffer
    if (gpu_buffer_fmi_asearch->num_queries > 0) {
      return false; // Leave it to the next fresh buffer
    }
    // Reallocate buffer
    gpu_fmi_asearch_init_and_realloc_buffer_(gpu_buffer_fmi_asearch->buffer,
        REGION_MAX_REGIONS_FACTOR,total_bases,total_queries,total_regions);
    // Check reallocated buffer dimensions (error otherwise)
    max_queries = gpu_buffer_fmi_asearch_get_max_queries(gpu_buffer_fmi_asearch);
    gem_cond_fatal_error(total_queries > max_queries,GPU_FMI_ASEARCH_MAX_QUERIES,total_queries,max_queries);
    max_bases = gpu_buffer_fmi_asearch_get_max_bases(gpu_buffer_fmi_asearch);
    gem_cond_fatal_error(total_bases > max_bases,GPU_FMI_ASEARCH_MAX_PATTERN_LENGTH,total_bases,max_bases);
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
  // Add query
  const uint64_t num_queries = gpu_buffer_fmi_asearch->num_queries;
  gpu_fmi_search_query_info_t* const search_query_info =
      gpu_fmi_asearch_buffer_get_queries_info_(gpu_buffer_fmi_asearch->buffer) + num_queries;
  search_query_info->init_offset = gpu_buffer_fmi_asearch->num_bases;
  search_query_info->query_size = key_length;
  ++(gpu_buffer_fmi_asearch->num_queries); // Next
  // Add pattern (bases)
  gpu_fmi_search_query_t* const search_query =
      gpu_fmi_asearch_buffer_get_queries_(gpu_buffer_fmi_asearch->buffer) +
      gpu_buffer_fmi_asearch->num_bases;
  const uint8_t* const key = pattern->key;
  memcpy(search_query,key,key_length); // Copy pattern
  gpu_buffer_fmi_asearch->num_bases += key_length; // Next
  // Set resulting regions offsets
  const uint64_t num_regions = gpu_buffer_fmi_asearch->num_regions;
  gpu_fmi_search_region_t* const search_region =
      gpu_fmi_asearch_buffer_get_regions_(gpu_buffer_fmi_asearch->buffer) + num_queries;
  search_region->init_offset = num_regions;
  gpu_buffer_fmi_asearch->num_regions += max_regions; // Next
  // Return query-offset
  return num_queries;
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
    const uint64_t max_regions) { GEM_CUDA_NOT_SUPPORTED(); }
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
