/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_BUFFER_COLLECTION_H_
#define GPU_BUFFER_COLLECTION_H_

#include "essentials.h"
#include "archive.h"
#include "profiler_timer.h"
#include "gpu_config.h"

typedef struct {
  void* gpu_buffers_dto;              // GPU-Buffer Initializer DTO
  void** internal_buffers;            // Internal Buffers
  uint64_t num_buffers;               // Total number of buffers allocated
} gpu_buffer_collection_t;

/*
 * Setup
 */
gpu_buffer_collection_t* gpu_buffer_collection_new(
    char* const gpu_index_name,const uint64_t num_buffers,
    const uint64_t buffer_size,const bool verbose);
void gpu_buffer_collection_delete(gpu_buffer_collection_t* const gpu_buffer_collection);

/*
 * Accessors
 */
void* gpu_buffer_collection_get_buffer(
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_no);

#endif /* GPU_BUFFER_COLLECTION_H_ */
