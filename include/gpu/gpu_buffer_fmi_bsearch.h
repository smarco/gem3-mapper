/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_bsearch.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_BUFFER_FMI_BSEARCH_H_
#define GPU_BUFFER_FMI_BSEARCH_H_

#include "utils/essentials.h"
#include "system/profiler_timer.h"
#include "data_structures/pattern.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU FM-Index Backward Search Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;             // GPU Generic Buffer
  /* Buffer state */
  uint32_t num_queries;     // Buffer state
  /* CPU Computation */
  bool compute_cpu;         // Computing Using CPU (disable GPU)
  fm_index_t* fm_index;     // FM-Index
  /* Profile */
  gem_timer_t timer;
} gpu_buffer_fmi_search_t;

/*
 * Setup
 */
gpu_buffer_fmi_search_t* gpu_buffer_fmi_search_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    fm_index_t* const fm_index);
void gpu_buffer_fmi_search_clear(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);
void gpu_buffer_fmi_search_delete(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);

/*
 * Computing Device
 */
void gpu_buffer_fmi_search_set_device_cpu(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);
void gpu_buffer_fmi_search_set_device_gpu(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_search_get_max_queries(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);
uint64_t gpu_buffer_fmi_search_get_num_queries(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);
bool gpu_buffer_fmi_search_fits_in_buffer(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,
    const uint64_t total_queries);

/*
 * Accessors
 */
void gpu_buffer_fmi_search_add_query(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,
    pattern_t* const pattern,
    const uint64_t begin,
    const uint64_t end);
void gpu_buffer_fmi_search_get_result(
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search,
    const uint64_t buffer_pos,
    uint64_t* const hi,
    uint64_t* const lo);

/*
 * Send/Receive
 */
void gpu_buffer_fmi_search_send(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);
void gpu_buffer_fmi_search_receive(gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);

#endif /* GPU_BUFFER_FMI_BSEARCH_H_ */
