/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_decode.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_BUFFER_FMI_DECODE_H_
#define GPU_BUFFER_FMI_DECODE_H_

#include "utils/essentials.h"
#include "system/profiler_timer.h"
#include "gpu/gpu_buffer_collection.h"

/*
 * GPU FM-Index Backward Search Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  void* buffer;           // GPU Generic Buffer
  /* Buffer state */
  uint32_t num_queries;   // Buffer state
  /* GPU Computation */
  bool gpu_decode_sa;     // GPU Computing decode-SA
  bool gpu_decode_text;   // GPU Computing decode-Text
  fm_index_t* fm_index;   // FM-Index
  /* Profile */
  gem_timer_t timer;
} gpu_buffer_fmi_decode_t;

/*
 * Setup
 */
gpu_buffer_fmi_decode_t* gpu_buffer_fmi_decode_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_pos,
    fm_index_t* const fm_index,
    const bool gpu_decode_sa,
    const bool gpu_decode_text);
void gpu_buffer_fmi_decode_clear(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void gpu_buffer_fmi_decode_delete(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
uint64_t gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
bool gpu_buffer_fmi_decode_fits_in_buffer(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t total_queries);

/*
 * Accessors
 */
void gpu_buffer_fmi_decode_add_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t bwt_position);
void gpu_buffer_fmi_decode_get_position_sa(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const bwt_sampled_position,
    uint64_t* const lf_steps);
void gpu_buffer_fmi_decode_get_position_text(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_pos,
    uint64_t* const text_position);

/*
 * Send/Receive
 */
void gpu_buffer_fmi_decode_send(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void gpu_buffer_fmi_decode_receive(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

#endif /* GPU_BUFFER_FMI_DECODE_H_ */
