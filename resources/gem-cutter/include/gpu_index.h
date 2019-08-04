/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_INDEX_H_
#define GPU_INDEX_H_

#include "gpu_commons.h"
#include "gpu_devices.h"
#include "gpu_index_modules.h"


/* Get information functions */
gpu_error_t gpu_index_get_size(const gpu_index_buffer_t* const index, size_t *bytesPerIndex, const gpu_module_t activeModules);

/* Functions to initialize the index data on the DEVICE */
gpu_error_t gpu_index_init_dto(gpu_index_buffer_t *index, const gpu_module_t activeModules);
gpu_error_t gpu_index_init(gpu_index_buffer_t** const index, const gpu_index_dto_t* const rawIndex, const uint32_t numSupportedDevices, const gpu_module_t activeModules);
gpu_error_t gpu_index_load(gpu_index_buffer_t* index, const gpu_index_dto_t * const rawIndex,const gpu_module_t activeModules);
gpu_error_t gpu_index_set_specs(gpu_index_buffer_t* const index, const gpu_index_dto_t* const indexRaw,const gpu_index_coding_t indexCoding, const gpu_module_t activeModules);
gpu_error_t gpu_index_allocate(gpu_index_buffer_t* index, const gpu_module_t activeModules);

/* Data transfer functions */
gpu_error_t gpu_index_transfer_CPU_to_GPUs(gpu_index_buffer_t* const index, gpu_device_info_t** const devices, const gpu_module_t activeModules);

/* Stream index functions  */
gpu_error_t gpu_index_read_specs(int fp, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_index_read(int fp, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_index_write_specs(int fp, const gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_index_write(int fp, const gpu_index_buffer_t* const index, const gpu_module_t activeModules);

/* Stream index functions  */
gpu_error_t gpu_index_transform(gpu_index_buffer_t* const index, const gpu_index_dto_t* const indexRaw, const gpu_index_coding_t indexCoding, const gpu_module_t activeModules);
gpu_error_t gpu_index_transform_ASCII(const gpu_index_dto_t* const textRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_index_transform_GEM_FULL(const gpu_index_dto_t* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_index_transform_MFASTA_FULL(const gpu_index_dto_t* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_index_load_specs_MFASTA_FULL(const gpu_index_dto_t* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules);

/* Functions to release the index data from the DEVICE & HOST */
gpu_error_t gpu_index_free_host(gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_index_free_unused_host(gpu_index_buffer_t* index, gpu_device_info_t** const devices, const gpu_module_t activeModules);
gpu_error_t gpu_index_free_device(gpu_index_buffer_t* index, gpu_device_info_t** const devices, const gpu_module_t activeModules);
gpu_error_t gpu_index_free_metainfo(gpu_index_buffer_t* index, const gpu_module_t activeModules);
gpu_error_t gpu_index_free(gpu_index_buffer_t **index, gpu_device_info_t** const devices, const gpu_module_t activeModules);

#endif /* GPU_INDEX_H_ */
