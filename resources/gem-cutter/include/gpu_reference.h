/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_REFERENCE_H_
#define GPU_REFERENCE_H_

#include "gpu_commons.h"
#include "gpu_devices.h"

/* Defines of global reference representation */
#define GPU_REFERENCE_PLAIN__CHAR_LENGTH       2
#define GPU_REFERENCE_PLAIN__ENTRY_LENGTH      GPU_UINT64_LENGTH
#define GPU_REFERENCE_PLAIN__ENTRY_SIZE        GPU_UINT64_SIZE
#define GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY  (GPU_REFERENCE_PLAIN__ENTRY_LENGTH / GPU_REFERENCE_PLAIN__CHAR_LENGTH)
#define GPU_REFERENCE_PLAIN__MASK_BASE        (GPU_UINT64_ONES >> (GPU_REFERENCE_PLAIN__ENTRY_LENGTH - GPU_REFERENCE_PLAIN__CHAR_LENGTH))

/* Defines of masked reference representation */
#define GPU_REFERENCE_MASKED__CHAR_LENGTH      1
#define GPU_REFERENCE_MASKED__ENTRY_LENGTH     GPU_UINT64_LENGTH
#define GPU_REFERENCE_MASKED__ENTRY_SIZE       GPU_UINT64_SIZE
#define GPU_REFERENCE_MASKED__CHARS_PER_ENTRY (GPU_REFERENCE_MASKED__ENTRY_LENGTH / GPU_REFERENCE_MASKED__CHAR_LENGTH)
#define GPU_REFERENCE_MASKED__MASK_BASE       (GPU_UINT64_ONES >> (GPU_REFERENCE_MASKED__ENTRY_LENGTH - GPU_REFERENCE_MASKED__CHAR_LENGTH))

/* Defines of global reference representation */
#define GPU_REFERENCE_CHAR_LENGTH      GPU_REFERENCE_PLAIN__CHAR_LENGTH
#define GPU_REFERENCE_ENTRY_LENGTH     GPU_REFERENCE_PLAIN__ENTRY_LENGTH
#define GPU_REFERENCE_ENTRY_SIZE       GPU_REFERENCE_PLAIN__ENTRY_SIZE
#define GPU_REFERENCE_CHARS_PER_ENTRY  GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY
#define GPU_REFERENCE_MASK_BASE        GPU_REFERENCE_PLAIN__MASK_BASE
#define GPU_REFERENCE_UINT64_MASK_BASE (GPU_UINT64_ONES >> (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH))
#define GPU_REFERENCE_UINT32_MASK_BASE (GPU_UINT32_ONES >> (GPU_UINT32_LENGTH - GPU_REFERENCE_CHAR_LENGTH))
#define GPU_REFERENCE_END_PADDING      625

/*****************************
Internal Objects (General)
*****************************/

typedef struct {
  /* Data types for a Global Reference */
  uint64_t        size;
  /* Data types for a Plain Reference */
  uint64_t        numEntriesPlain;
  uint64_t        *h_reference_plain;
  uint64_t        **d_reference_plain;
  /* Data types for a Masked Reference */
  uint64_t        numEntriesMasked;
  uint64_t		  *h_reference_masked;
  uint64_t        **d_reference_masked;
  /* Memory allocation configuration */
  memory_stats_t  hostAllocStats;
  memory_alloc_t  *memorySpace;
  gpu_module_t    activeModules;
} gpu_reference_buffer_t;


/* Get information functions */
size_t		gpu_reference_get_size(gpu_reference_buffer_t* const reference, size_t* const referenceSize, const gpu_module_t activeModules);

/* String basic functions */
uint64_t    gpu_char_to_bin_ASCII(const unsigned char base);
char        gpu_complement_base(const char character);

/* Transform reference functions */
gpu_error_t gpu_reference_transform(gpu_reference_buffer_t* const ref, const char* const referenceRaw, const gpu_ref_coding_t refCoding, const gpu_module_t activeModules);

/* Stream reference functions  */
gpu_error_t gpu_reference_read_specs(int fp, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_read(int fp, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_write_specs(int fp, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_write(int fp, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);

/* Initialize reference functions */
gpu_error_t gpu_reference_init_dto(gpu_reference_buffer_t* const ref);
gpu_error_t gpu_reference_set_specs(gpu_reference_buffer_t* const ref, const char* const referenceRaw, const gpu_ref_coding_t refCoding, const gpu_module_t activeModules);
gpu_error_t gpu_reference_init(gpu_reference_buffer_t **reference, const gpu_reference_dto_t* const referenceRaw, const uint32_t numSupportedDevices, const gpu_module_t activeModules);
gpu_error_t gpu_reference_load(gpu_reference_buffer_t *reference, const gpu_reference_dto_t* const referenceRaw,const gpu_module_t activeModules);
gpu_error_t gpu_reference_allocate(gpu_reference_buffer_t *reference, const gpu_module_t activeModules);

/* Data transfer functions */
gpu_error_t gpu_reference_transfer_CPU_to_GPUs(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices,const gpu_module_t activeModules);

/* Free reference functions */
gpu_error_t gpu_reference_free_host(gpu_reference_buffer_t* const reference);
gpu_error_t gpu_reference_free_unused_host(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices, const gpu_module_t activeModules);
gpu_error_t gpu_reference_free_device(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices);
gpu_error_t gpu_reference_free(gpu_reference_buffer_t **reference, gpu_device_info_t** const devices, const gpu_module_t activeModules);


/* LOCAL functions */
gpu_error_t gpu_reference_transform_ASCII(const char* const referenceASCII, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_transform_GEM(const gpu_gem_ref_dto_t* const gem_reference, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_reference_transform_GEM_FULL(const gpu_gem_ref_dto_t* const gem_reference, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);

#endif /* GPU_REFERENCE_H_ */
