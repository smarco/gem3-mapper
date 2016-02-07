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
  memory_alloc_t  *memorySpace;
} gpu_fmi_buffer_t;


/* Functions to initialize the index data on the DEVICE */
gpu_error_t gpu_init_fmi_index_dto(gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_init_fmi_index(gpu_fmi_buffer_t* const fmi, const char* const indexRaw,
                               const uint64_t bwtSize, const gpu_index_coding_t indexCoding,
                               const uint32_t numSupportedDevices);

/* Data transfer functions */
gpu_error_t gpu_transfer_fmi_index_CPU_to_GPUs(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices);

/* Stream index functions  */
gpu_error_t gpu_read_fmi_index(FILE* fp, gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_write_fmi_index(FILE* fp, const gpu_fmi_buffer_t* const fmi);

/* Data transform functions  */
gpu_error_t gpu_transform_fmi_index_ASCII(const char* const textBWT, gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_transform_fmi_index_GEM_FULL(const gpu_gem_fmi_dto_t* const gpu_gem_fmi_dto, gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_transform_fmi_index_MFASTA_FULL(const char* const indexRaw, gpu_fmi_buffer_t* const fmi);

/* Functions to release the index data from the DEVICE & HOST */
gpu_error_t gpu_free_fmi_index_host(gpu_fmi_buffer_t* const fmi);
gpu_error_t gpu_free_unused_fmi_index_host(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices);
gpu_error_t gpu_free_fmi_index_device(gpu_fmi_buffer_t* const fmi, gpu_device_info_t** const devices);
gpu_error_t gpu_free_fmi_index_metainfo(gpu_fmi_buffer_t* const fmi);

/* Local functions to transform the data */
gpu_error_t gpu_fmi_index_build_PEQ(const gpu_fmi_buffer_t* const fmi, const char* const h_ascii_BWT,
                                    gpu_index_bitmap_entry_t* const h_bitmap_BWT);
gpu_error_t gpu_fmi_index_build_COUNTERS(const gpu_fmi_buffer_t* const fmi, gpu_index_counter_entry_t* const h_counters_FMI,
                                         const char* const h_ascii_BWT);
gpu_error_t gpu_fmi_index_build_FMI(gpu_fmi_buffer_t* const fmi, gpu_index_bitmap_entry_t* const h_bitmap_BWT,
                                    const gpu_index_counter_entry_t* const h_counters_FMI);

#endif /* GPU_FMI_INDEX_H_ */
