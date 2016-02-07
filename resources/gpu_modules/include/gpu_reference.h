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
  memory_alloc_t  *memorySpace;
  gpu_module_t    activeModules;
} gpu_reference_buffer_t;

/* String basic functions */
uint64_t    gpu_char_to_bin_ASCII(const unsigned char base);
char        gpu_complement_base(const char character);

/* Transform reference functions */
gpu_error_t gpu_transform_reference(const char* const referenceRaw, gpu_reference_buffer_t* const ref, const gpu_ref_coding_t refCoding);

/* Stream reference functions  */
gpu_error_t gpu_read_reference(FILE* fp, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_write_reference(FILE* fp, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);


/* Initialize reference functions */
gpu_error_t gpu_init_reference_dto(gpu_reference_buffer_t* const ref);
gpu_error_t gpu_init_reference(gpu_reference_buffer_t **reference, const char* const referenceRaw,
                               const uint64_t refSize, const gpu_ref_coding_t refCoding,
                               const uint32_t numSupportedDevices, const gpu_module_t activeModules);
gpu_error_t gpu_transfer_reference_CPU_to_GPUs(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices);

/* Free reference functions */
gpu_error_t gpu_free_reference_host(gpu_reference_buffer_t* const reference);
gpu_error_t gpu_free_unused_reference_host(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices);
gpu_error_t gpu_free_reference_device(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices);
gpu_error_t gpu_free_reference(gpu_reference_buffer_t **reference, gpu_device_info_t** const devices);

#endif /* GPU_REFERENCE_H_ */
