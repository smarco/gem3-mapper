/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Common headers and data structures for BPM on GPU library
 */

#include "gpu_commons.h"

/* Include the required objects */
#include "gpu_devices.h"
#include "gpu_reference.h"
#include "gpu_index.h"

/* Include all the supported modules */
#include "gpu_fmi.h"
#include "gpu_bpm.h"

#ifndef GPU_BUFFER_H_
#define GPU_BUFFER_H_

typedef union{
	gpu_bpm_buffer_t 		bpm;
	gpu_fmi_search_buffer_t search;
	gpu_fmi_decode_buffer_t decode;
} gpu_data_buffer_t;

typedef struct {
	gpu_module_t			typeBuffer;
	uint32_t				numBuffers;
	uint32_t 				idBuffer;
	cudaStream_t			idStream;
	uint32_t				idSupportedDevice;
	gpu_device_info_t		**device;
	gpu_reference_buffer_t	*reference;
	gpu_index_buffer_t		*index;
	size_t					sizeBuffer;
	void 					*h_rawData;
	void 					*d_rawData;
	gpu_data_buffer_t		data;
} gpu_buffer_t;


/* Primitives to get information from buffers */
gpu_error_t gpu_get_min_memory_per_module(size_t *minimumMemorySize, gpu_reference_buffer_t *reference,
													 gpu_index_buffer_t *index,  const uint32_t numBuffers,
													 const gpu_module_t activeModules);
gpu_error_t gpu_module_memory_manager(const uint32_t idDevice, const uint32_t idSupDevice, const uint32_t numBuffers,
												 const gpu_data_location_t userReferenceAllocOption, gpu_module_t activeModules,
												 bool *dataFits, bool *lReference, bool *lIndex, size_t *recMemSize,
												 size_t *reqMemSize, gpu_reference_buffer_t *reference, gpu_index_buffer_t *index);
gpu_error_t gpu_configure_modules(gpu_device_info_t ***devices, const gpu_dev_arch_t selectedArchitectures,
						  					 const gpu_data_location_t userAllocOption, const gpu_module_t activeModules,
						  					 const uint32_t numBuffers, gpu_reference_buffer_t *reference, gpu_index_buffer_t *index);

/* Primitives to schedule and manage the buffers */
gpu_error_t gpu_get_min_memory_size_per_buffer(size_t *bytesPerBuffer); //gpu_get_min_memory_size_per_device
gpu_error_t gpu_configure_buffer(gpu_buffer_t *mBuff, const uint32_t idBuffer, const uint32_t idSupportedDevice, const size_t bytesPerBuffer,
											const uint32_t numBuffers, gpu_device_info_t **device, gpu_reference_buffer_t *reference, gpu_index_buffer_t *index);
gpu_error_t gpu_schedule_buffers(gpu_buffer_t ***gpuBuffer, const uint32_t numBuffers, gpu_device_info_t **device,
										   gpu_reference_buffer_t *reference, gpu_index_buffer_t *index,
										   float maxMbPerBuffer, const gpu_module_t activeModules);

/* Functions to free all the buffer resources (HOST & DEVICE) */
gpu_error_t gpu_free_buffer(gpu_buffer_t *mBuff);

#endif /* GPU_BUFFER_H_ */
