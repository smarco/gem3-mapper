#ifndef GPU_INDEX_H_
#define GPU_INDEX_H_

#include "gpu_commons.h"
#include "gpu_devices.h"
#include "gpu_index_modules.h"

/* Functions to initialize the index data on the DEVICE */
gpu_error_t gpu_init_index_dto(gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_init_index(gpu_index_buffer_t **index, const char* const indexRaw,
                           const uint64_t bwtSize, const uint32_t samplingRate, const gpu_index_coding_t indexCoding,
                           const uint32_t numSupportedDevices, const gpu_module_t activeModules);

/* Data transfer functions */
gpu_error_t gpu_transfer_index_CPU_to_GPUs(gpu_index_buffer_t* const index, gpu_device_info_t** const devices, const gpu_module_t activeModules);

/* Stream index functions  */
gpu_error_t gpu_write_index(FILE* fp, const gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_read_index(FILE* fp, gpu_index_buffer_t* const index, const gpu_module_t activeModules);

/* Stream index functions  */
gpu_error_t gpu_transform_index_ASCII(const char* const textRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_transform_index_GEM_FULL(const char* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_transform_index_MFASTA_FULL(const char* const indexRaw, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_transform_index(const char* const indexRaw, gpu_index_buffer_t* const index, const gpu_index_coding_t indexCoding, const gpu_module_t activeModules);
/* Functions to release the index data from the DEVICE & HOST */
gpu_error_t gpu_free_index_host(gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_free_unused_index_host(gpu_index_buffer_t* index, gpu_device_info_t** const devices, const gpu_module_t activeModules);
gpu_error_t gpu_free_index_device(gpu_index_buffer_t* index, gpu_device_info_t** const devices, const gpu_module_t activeModules);
gpu_error_t gpu_free_index_metainfo(gpu_index_buffer_t* index, const gpu_module_t activeModules);
gpu_error_t gpu_free_index(gpu_index_buffer_t **index, gpu_device_info_t** const devices, const gpu_module_t activeModules);

#endif /* GPU_INDEX_H_ */
