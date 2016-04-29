/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_INDEX_H_
#define GPU_FMI_INDEX_H_

#include "gpu_commons.h"
#include "gpu_devices.h"
#include "gpu_fmi_structure.h"

/*****************************
Internal Objects (General)
*****************************/

typedef struct {
  uint64_t        bwtSize;
  uint64_t        numEntries;
  gpu_fmi_entry_t *h_fmi;
  gpu_fmi_entry_t **d_fmi;
  memory_stats_t  hostAllocStats;
  memory_alloc_t  *memorySpace;
} gpu_fmi_buffer_t;


/* Get information functions */
gpu_error_t gpu_fmi_index_get_size(const gpu_fmi_buffer_t* const __restrict__ fmi, size_t* const __restrict__ bytesPerFMI);

/* Functions to initialize the index data on the DEVICE */
gpu_error_t gpu_fmi_index_init_dto(gpu_fmi_buffer_t* const __restrict__ fmi);
gpu_error_t gpu_fmi_index_init(gpu_fmi_buffer_t* const __restrict__ fmi, const uint64_t bwtSize, const uint32_t numSupportedDevices);
gpu_error_t gpu_fmi_index_allocate(gpu_fmi_buffer_t* const __restrict__ fmi);

/* Data transfer functions */
gpu_error_t gpu_fmi_index_transfer_CPU_to_GPUs(gpu_fmi_buffer_t* const __restrict__ fmi, gpu_device_info_t** const __restrict__ devices);

/* Stream index functions */
gpu_error_t gpu_fmi_index_read_specs(FILE* fp, gpu_fmi_buffer_t* const __restrict__ fmi);
gpu_error_t gpu_fmi_index_read(FILE* fp, gpu_fmi_buffer_t* const __restrict__ fmi);
gpu_error_t gpu_fmi_index_write(FILE* fp, const gpu_fmi_buffer_t* const __restrict__ fmi);

/* Data transform functions */
gpu_error_t gpu_fmi_index_transform_ASCII(const char* const __restrict__ textBWT, gpu_fmi_buffer_t* const __restrict__ fmi);
gpu_error_t gpu_fmi_index_transform_GEM_FULL(const gpu_gem_fmi_dto_t* const __restrict__ gpu_gem_fmi_dto, gpu_fmi_buffer_t* const __restrict__ fmi);
gpu_error_t gpu_fmi_index_transform_MFASTA_FULL(const char* const __restrict__ indexRaw, gpu_fmi_buffer_t* const __restrict__ fmi);

/* Data load functions */
gpu_error_t gpu_fmi_index_load_specs_MFASTA_FULL(const char* const __restrict__ fn, gpu_fmi_buffer_t* const __restrict__ fmi);
gpu_error_t gpu_fmi_index_load_MFASTA_FULL(const char* const __restrict__ fn, gpu_fmi_buffer_t* const __restrict__ fmi, char **h_BWT);

/* Functions to release the index data from the DEVICE & HOST */
gpu_error_t gpu_fmi_index_free_host(gpu_fmi_buffer_t* const __restrict__ fmi);
gpu_error_t gpu_fmi_index_free_unused_host(gpu_fmi_buffer_t* const __restrict__ fmi, gpu_device_info_t** const __restrict__ devices);
gpu_error_t gpu_fmi_index_free_device(gpu_fmi_buffer_t* const __restrict__ fmi, gpu_device_info_t** const __restrict__ devices);
gpu_error_t gpu_fmi_index_free_metainfo(gpu_fmi_buffer_t* const __restrict__ fmi);

/* Local functions to transform the data */
gpu_error_t gpu_fmi_index_build_PEQ(const gpu_fmi_buffer_t* const __restrict__ fmi, const char* const __restrict__ h_ascii_BWT,
                                    gpu_index_bitmap_entry_t* const __restrict__ h_bitmap_BWT);
gpu_error_t gpu_fmi_index_build_COUNTERS(const gpu_fmi_buffer_t* const __restrict__ fmi, gpu_index_counter_entry_t* const __restrict__ h_counters_FMI,
                                         const char* const __restrict__ h_ascii_BWT);
gpu_error_t gpu_fmi_index_build_FMI(gpu_fmi_buffer_t* const __restrict__ fmi, gpu_index_bitmap_entry_t* const __restrict__ h_bitmap_BWT,
                                    const gpu_index_counter_entry_t* const __restrict__ h_counters_FMI);

#endif /* GPU_FMI_INDEX_H_ */
