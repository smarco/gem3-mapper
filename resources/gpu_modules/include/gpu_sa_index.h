#ifndef GPU_SA_H_
#define GPU_SA_H_

#include "gpu_commons.h"
#include "gpu_devices.h"

typedef struct {
  uint64_t        sampligRate;
  uint64_t        numEntries;
  gpu_sa_entry_t  *h_sa;
  gpu_sa_entry_t  **d_sa;
  memory_alloc_t  *memorySpace;
} gpu_sa_buffer_t;


/* Functions to initialize the index data on the DEVICE */
gpu_error_t gpu_init_sa_index_dto(gpu_sa_buffer_t* const sa);
gpu_error_t gpu_init_sa_index(gpu_sa_buffer_t* const sa, const char* const indexRaw, const uint64_t textSize,
                              const uint32_t sampligRate, const gpu_index_coding_t indexCoding,
                              const uint32_t numSupportedDevices);

/* Data transfer functions */
gpu_error_t gpu_transfer_sa_index_CPU_to_GPUs(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices);

/* Stream index functions  */
gpu_error_t gpu_write_sa_index(FILE* fp, const gpu_sa_buffer_t* const sa);
gpu_error_t gpu_read_sa_index(FILE* fp, gpu_sa_buffer_t* const sa);

/* Data transform functions  */
gpu_error_t gpu_transform_sa_index_ASCII(const char* const textBWT, gpu_sa_buffer_t* const sa);
gpu_error_t gpu_transform_sa_index_GEM_FULL(const gpu_gem_sa_dto_t* const gpu_gem_sa_dto, gpu_sa_buffer_t* const sa);
gpu_error_t gpu_transform_sa_index_MFASTA_FULL(const char* const indexRaw, gpu_sa_buffer_t* const sa);

/* Functions to release the index data from the DEVICE & HOST */
gpu_error_t gpu_free_sa_index_host(gpu_sa_buffer_t* const sa);
gpu_error_t gpu_free_unused_sa_index_host(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices);
gpu_error_t gpu_free_sa_index_device(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices);
gpu_error_t gpu_free_sa_index_metainfo(gpu_sa_buffer_t* const sa);

#endif /* GPU_SA_H_ */
