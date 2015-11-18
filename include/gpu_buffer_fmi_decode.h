/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_fmi_decode.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_BUFFER_FMI_DECODE_H_
#define GPU_BUFFER_FMI_DECODE_H_

#include "essentials.h"
#include "profiler_timer.h"
#include "gpu_buffer_collection.h"

/*
 * GPU FM-Index Backward Search Buffer
 */
typedef struct {
  /* GPU Generic Buffer*/
  device_t device;          // Computing Device (GPU/CPU)
  void* buffer;             // GPU Generic Buffer
  /* Buffer state */
  uint32_t num_queries;     // Buffer state
  /* Profile */
  gem_timer_t timer;
} gpu_buffer_fmi_decode_t;

/*
 * Setup
 */
gpu_buffer_fmi_decode_t* gpu_buffer_fmi_decode_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_pos);
void gpu_buffer_fmi_decode_clear(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void gpu_buffer_fmi_decode_delete(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void gpu_buffer_fmi_decode_device(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const device_t device);

/*
 * Occupancy & Limits
 */
uint64_t gpu_buffer_fmi_decode_get_max_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
uint64_t gpu_buffer_fmi_decode_get_num_queries(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

/*
 * Accessors
 */
void gpu_buffer_fmi_decode_add_query(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t bwt_position);
void gpu_buffer_fmi_decode_get_result(
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t buffer_pos,
    uint64_t* const bwt_sampled_position,uint64_t* const lf_steps);

/*
 * Send/Receive
 */
void gpu_buffer_fmi_decode_send(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void gpu_buffer_fmi_decode_receive(gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

/*
 * Errors
 */
#define GEM_ERROR_GPU_FMI_DECODE_MAX_QUERIES "GPU-FMI-Decode. Number of queries (%"PRIu64") exceeds maximum buffer capacity (%"PRIu64" decodes)"

#endif /* GPU_BUFFER_FMI_DECODE_H_ */
