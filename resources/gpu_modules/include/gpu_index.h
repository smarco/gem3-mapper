#ifndef GPU_INDEX_H_
#define GPU_INDEX_H_

#include "gpu_commons.h"
#include "gpu_devices.h"

/* Defines related to BWT representation */
#define GPU_FMI_BITS_LOAD_PER_THREAD	128			// accesses of 16 Bytes / thread

#define GPU_FMI_ENTRY_LENGTH			((GPU_FMI_COUNTERS_PER_ENTRY * GPU_UINT64_LENGTH) \
										+ (GPU_FMI_ENTRY_SIZE * GPU_FMI_BWT_CHAR_LENGTH)) 		//   512 bits
#define GPU_FMI_THREADS_PER_ENTRY		(GPU_FMI_ENTRY_LENGTH / GPU_FMI_BITS_LOAD_PER_THREAD)  	//   4 threads
#define GPU_FMI_ENTRIES_PER_WARP		(GPU_WARP_SIZE / GPU_FMI_THREADS_PER_ENTRY)
#define GPU_FMI_ENTRIES_PER_BLOCK		(GPU_MAX_THREADS_PER_SM / GPU_FMI_THREADS_PER_ENTRY)
#define GPU_FMI_SLICES_PER_THREAD		(GPU_FMI_BITS_LOAD_PER_THREAD / GPU_UINT32_LENGTH)		// slices of 16 Bytes / thread


/*****************************
Internal Objects (General)
*****************************/

typedef struct {
	uint64_t        bwtSize;
	uint64_t        numEntries;
	gpu_fmi_entry_t *h_fmi;
	gpu_fmi_entry_t **d_fmi;
	memory_alloc_t	*memorySpace;
} gpu_index_buffer_t;


/*************************************
Specific types for the Devices (GPUs)
**************************************/

// Just to cast and exchange memory between threads for the shared mem
typedef union {
    uint4 s;
    uint32_t v[GPU_FMI_SLICES_PER_THREAD];
} gpu_fmi_exch_bmp_mem_t;

typedef union {
	uint4 v[GPU_FMI_THREADS_PER_ENTRY];
} gpu_fmi_device_entry_t;

/* Functions to initialize the index data on the DEVICE*/
gpu_error_t gpu_load_index_PROFILE(const char *fn, gpu_index_buffer_t *index);
gpu_error_t gpu_transfer_index_CPU_to_GPUs(gpu_index_buffer_t *index, gpu_device_info_t **devices);
gpu_error_t gpu_init_index(gpu_index_buffer_t **index, const char *indexRaw,
									  const uint64_t bwtSize, const gpu_index_coding_t indexCoding,
									  const uint32_t numSupportedDevices);
/* Functions to release the index data from the DEVICE & HOST*/
gpu_error_t gpu_free_index_host(gpu_index_buffer_t *index);
gpu_error_t gpu_free_unused_index_host(gpu_index_buffer_t *index, gpu_device_info_t **devices);
gpu_error_t gpu_free_index_device(gpu_index_buffer_t *index, gpu_device_info_t **devices);
gpu_error_t gpu_free_index(gpu_index_buffer_t **index, gpu_device_info_t **devices);


#endif /* GPU_INDEX_H_ */
