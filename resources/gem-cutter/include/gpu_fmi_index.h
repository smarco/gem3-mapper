/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
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
#include "gpu_fmi_table.h"

/*****************************
Internal Objects (General)
*****************************/

typedef struct {
  gpu_fmi_table_t table;
  uint64_t        bwtSize;
  uint64_t        numEntries;
  gpu_fmi_entry_t *h_fmi;
  gpu_fmi_entry_t **d_fmi;
  memory_stats_t  hostAllocStats;
  memory_alloc_t  *memorySpace;
} gpu_fmi_buffer_t;


/* Get information functions */
gpu_error_t gpu_fmi_index_get_size(const gpu_fmi_buffer_t* const fmi, size_t* const bytesPerFMI);

/* Functions to initialize the index data on the DEVICE */
gpu_error_t gpu_fmi_index_init_dto(gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_fmi_index_init(gpu_fmi_buffer_t* const fmi, const uint64_t bwtSize, const uint32_t numSupportedDevices);
gpu_error_t gpu_fmi_index_allocate(gpu_fmi_buffer_t* const fmi);

/* Data transfer functions */
gpu_error_t gpu_fmi_index_transfer_CPU_to_GPUs(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices);

/* Stream index functions */
gpu_error_t gpu_fmi_index_read_specs(int fp, gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_fmi_index_read(int fp, gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_fmi_index_write_specs(int fp, const gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_fmi_index_write(int fp, const gpu_fmi_buffer_t* const fmi);

/* Data transform functions */
gpu_error_t gpu_fmi_index_transform_ASCII(const char* const textBWT, gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_fmi_index_transform_GEM_FULL(const gpu_gem_fmi_dto_t* const gpu_gem_fmi_dto, gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_fmi_index_transform_MFASTA_FULL(const char* const indexRaw, gpu_fmi_buffer_t* const fmi);

/* Data load functions */
gpu_error_t gpu_fmi_index_load_specs_MFASTA_FULL(const char* const fn, gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_fmi_index_load_MFASTA_FULL(const char* const fn, gpu_fmi_buffer_t* const fmi, char **h_BWT);

/* Functions to release the index data from the DEVICE & HOST */
gpu_error_t gpu_fmi_index_free_host(gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_fmi_index_free_unused_host(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices);
gpu_error_t gpu_fmi_index_free_device(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices);
gpu_error_t gpu_fmi_index_free_metainfo(gpu_fmi_buffer_t* const fmi);

/* Local functions to transform the data */
gpu_error_t gpu_fmi_index_build_PEQ(const gpu_fmi_buffer_t* const fmi, const char* const h_ascii_BWT,
                                    gpu_index_bitmap_entry_t* const h_bitmap_BWT);
gpu_error_t gpu_fmi_index_build_COUNTERS(const gpu_fmi_buffer_t* const fmi, gpu_index_counter_entry_t* const h_counters_FMI,
                                         const char* const h_ascii_BWT);
gpu_error_t gpu_fmi_index_build_FMI(gpu_fmi_buffer_t* const fmi, gpu_index_bitmap_entry_t* const h_bitmap_BWT,
                                    const gpu_index_counter_entry_t* const h_counters_FMI);

#endif /* GPU_FMI_INDEX_H_ */
