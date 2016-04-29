/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_REFERENCE_H_
#define GPU_REFERENCE_H_

#include "gpu_commons.h"
#include "gpu_devices.h"

/* Defines related to Reference representation */
#define GPU_REFERENCE_CHAR_LENGTH     2
#define GPU_REFERENCE_CHARS_PER_UINT1 (GPU_UINT32_LENGTH / GPU_REFERENCE_CHAR_LENGTH)
#define GPU_REFERENCE_CHARS_PER_UINT2 (GPU_REFERENCE_CHARS_PER_UINT1 * 2)
#define GPU_REFERENCE_CHARS_PER_UINT4 (GPU_REFERENCE_CHARS_PER_UINT1 * 4)
#define GPU_REFERENCE_CHARS_PER_ENTRY GPU_REFERENCE_CHARS_PER_UINT2
#define GPU_REFERENCE_BYTES_PER_ENTRY GPU_UINT64_SIZE
#define GPU_REFERENCE_MASK_BASE       (GPU_UINT64_ONES >> (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH))
#define GPU_REFERENCE_END_PADDING     625

/*****************************
Internal Objects (General)
*****************************/

typedef struct {
  uint64_t        size;
  uint64_t        numEntries;
  uint64_t        *h_reference;
  uint64_t        **d_reference;
  memory_stats_t  hostAllocStats;
  memory_alloc_t  *memorySpace;
  gpu_module_t    activeModules;
} gpu_reference_buffer_t;


/* Get information functions */
gpu_error_t gpu_reference_get_size(gpu_reference_buffer_t* const __restrict__ reference, size_t *bytesPerReference);

/* String basic functions */
uint64_t    gpu_char_to_bin_ASCII(const unsigned char base);
char        gpu_complement_base(const char character);

/* Transform reference functions */
gpu_error_t gpu_reference_transform(gpu_reference_buffer_t* const __restrict__ ref, const char* const __restrict__ referenceRaw, const gpu_ref_coding_t refCoding, const gpu_module_t activeModules);

/* Stream reference functions  */
gpu_error_t gpu_reference_read_specs(FILE* fp, gpu_reference_buffer_t* const __restrict__ reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_read(FILE* fp, gpu_reference_buffer_t* const __restrict__ reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_write(FILE* fp, const gpu_reference_buffer_t* const __restrict__ reference, const gpu_module_t activeModules);

/* Initialize reference functions */
gpu_error_t gpu_reference_init_dto(gpu_reference_buffer_t* const __restrict__ ref);
gpu_error_t gpu_reference_set_specs(gpu_reference_buffer_t* const __restrict__ ref, const char* const __restrict__ referenceRaw, const gpu_ref_coding_t refCoding, const gpu_module_t activeModules);
gpu_error_t gpu_reference_init(gpu_reference_buffer_t **reference, const gpu_reference_dto_t* const __restrict__ referenceRaw, const uint32_t numSupportedDevices, const gpu_module_t activeModules);
gpu_error_t gpu_reference_load(gpu_reference_buffer_t *reference, const gpu_reference_dto_t* const __restrict__ referenceRaw,const gpu_module_t activeModules);
gpu_error_t gpu_reference_allocate(gpu_reference_buffer_t *reference, const gpu_module_t activeModules);

/* Data transfer functions */
gpu_error_t gpu_reference_transfer_CPU_to_GPUs(gpu_reference_buffer_t* const __restrict__ reference, gpu_device_info_t** const __restrict__ devices,const gpu_module_t activeModules);

/* Free reference functions */
gpu_error_t gpu_reference_free_host(gpu_reference_buffer_t* const __restrict__ reference);
gpu_error_t gpu_reference_free_unused_host(gpu_reference_buffer_t* const __restrict__ reference, gpu_device_info_t** const __restrict__ devices, const gpu_module_t activeModules);
gpu_error_t gpu_reference_free_device(gpu_reference_buffer_t* const __restrict__ reference, gpu_device_info_t** const __restrict__ devices);
gpu_error_t gpu_reference_free(gpu_reference_buffer_t **reference, gpu_device_info_t** const __restrict__ devices, const gpu_module_t activeModules);


/* LOCAL functions */
gpu_error_t gpu_reference_transform_ASCII(const char* const __restrict__ referenceASCII, gpu_reference_buffer_t* const __restrict__ reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_transform_GEM(const gpu_gem_ref_dto_t* const __restrict__ gem_reference, gpu_reference_buffer_t* const __restrict__ reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_transform_GEM_FULL(const gpu_gem_ref_dto_t* const __restrict__ gem_reference, gpu_reference_buffer_t* const __restrict__ reference, const gpu_module_t activeModules);

#endif /* GPU_REFERENCE_H_ */
