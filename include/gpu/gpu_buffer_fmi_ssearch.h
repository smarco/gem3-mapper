/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_ssearch.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_BUFFER_FMI_SSEARCH_H_
#define GPU_BUFFER_FMI_SSEARCH_H_

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
  bool fmi_search_enabled;  // Enabled GPU fmi-search
  uint32_t num_queries;     // Buffer state
  /* Profile */
  gem_timer_t timer;
} gpu_buffer_fmi_ssearch_t;

/*
 * Setup
 */
gpu_buffer_fmi_ssearch_t* gpu_buffer_fmi_ssearch_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    const bool fmi_search_enabled);
void gpu_buffer_fmi_ssearch_clear(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void gpu_buffer_fmi_ssearch_delete(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_ssearch_get_max_queries(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
uint64_t gpu_buffer_fmi_ssearch_get_num_queries(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
bool gpu_buffer_fmi_ssearch_fits_in_buffer(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    const uint64_t total_queries);

/*
 * Accessors
 */
void gpu_buffer_fmi_ssearch_add_query(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    pattern_t* const pattern,
    const uint64_t begin,
    const uint64_t end);
void gpu_buffer_fmi_ssearch_get_result(
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch,
    const uint64_t buffer_pos,
    uint64_t* const hi,
    uint64_t* const lo);

/*
 * Send/Receive
 */
void gpu_buffer_fmi_ssearch_send(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void gpu_buffer_fmi_ssearch_receive(gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);

#endif /* GPU_BUFFER_FMI_SSEARCH_H_ */
