/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_SA_H_
#define GPU_SA_H_

#include "gpu_commons.h"
#include "gpu_devices.h"

typedef struct {
  uint64_t        sampligRate;
  uint64_t        numEntries;
  gpu_sa_entry_t  *h_sa;
  gpu_sa_entry_t  **d_sa;
  memory_stats_t  hostAllocStats;
  memory_alloc_t  *memorySpace;
} gpu_sa_buffer_t;


/* Get information functions */
gpu_error_t gpu_sa_index_get_size(const gpu_sa_buffer_t* const sa, size_t* const bytesPerSA);

/* Functions to initialize the index data on the DEVICE */
gpu_error_t gpu_sa_index_init_dto(gpu_sa_buffer_t* const sa);
gpu_error_t gpu_sa_index_init(gpu_sa_buffer_t* const sa, const uint64_t saNumEntries,
                              const uint32_t sampligRate, const uint32_t numSupportedDevices);
gpu_error_t gpu_sa_index_allocate(gpu_sa_buffer_t* const sa);

/* Data transfer functions */
gpu_error_t gpu_sa_index_transfer_CPU_to_GPUs(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices);

/* Stream index functions  */
gpu_error_t gpu_sa_index_read_specs(int fp, gpu_sa_buffer_t* const sa);
gpu_error_t gpu_sa_index_read(int fp, gpu_sa_buffer_t* const sa);
gpu_error_t gpu_sa_index_write_specs(int fp, const gpu_sa_buffer_t* const sa);
gpu_error_t gpu_sa_index_write(int fp, const gpu_sa_buffer_t* const sa);

/* Data load functions */
gpu_error_t gpu_sa_index_load_specs_MFASTA_FULL(const char* const indexRaw, gpu_sa_buffer_t* const sa);

/* Data transform functions  */
gpu_error_t gpu_sa_index_transform_ASCII(const char* const textBWT, gpu_sa_buffer_t* const sa);
gpu_error_t gpu_sa_index_transform_GEM_FULL(const gpu_gem_sa_dto_t* const gpu_gem_sa_dto, gpu_sa_buffer_t* const sa);
gpu_error_t gpu_sa_index_transform_MFASTA_FULL(const char* const indexRaw, gpu_sa_buffer_t* const sa);

/* Functions to release the index data from the DEVICE & HOST */
gpu_error_t gpu_sa_index_free_host(gpu_sa_buffer_t* const sa);
gpu_error_t gpu_sa_index_free_unused_host(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices);
gpu_error_t gpu_sa_index_free_device(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices);
gpu_error_t gpu_sa_index_free_metainfo(gpu_sa_buffer_t* const sa);

#endif /* GPU_SA_H_ */
