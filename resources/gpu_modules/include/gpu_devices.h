#ifndef GPU_DEVICES_H_
#define GPU_DEVICES_H_

#include "gpu_commons.h"

/*************************************
GPU Interface Objects
**************************************/

/* Defines related to GPU Architecture */
#define	GPU_WARP_SIZE					32
/* Defines related to GPU Architecture */
#define	GPU_THREADS_PER_BLOCK_FERMI		256
#define	GPU_THREADS_PER_BLOCK_KEPLER	128
#define	GPU_THREADS_PER_BLOCK_MAXWELL	64
#define	GPU_THREADS_PER_BLOCK_NEWGEN	64

typedef enum
{
	GPU_HOST_MAPPED,
	GPU_DEVICE_MAPPED,
	GPU_NONE_MAPPED
} memory_alloc_t;

typedef struct {
	/* System specifications */
	uint32_t		numDevices;
	uint32_t		numSupportedDevices;
	/* Device specifications */
	uint32_t		idDevice;
	uint32_t		idSupportedDevice;
	gpu_dev_arch_t	architecture;
	uint32_t		cudaCores;
	uint32_t		coreClockRate;   		// Mhz
	uint32_t		memoryBusWidth;  		// Bits
	uint32_t		memoryClockRate; 		// Mhz
	/* Device performance metrics */
	uint32_t		absolutePerformance;	// MOps/s
	float			relativePerformance;	// Ratio
	uint32_t		absoluteBandwidth; 		// MB/s
	float			relativeBandwidth;		// Ratio
	/* System performance metrics */
	uint32_t		allSystemPerformance;	// MOps/s
	uint32_t		allSystemBandwidth;		// MB/s
} gpu_device_info_t;

/* Primitives to get information for the scheduler */
size_t 			gpu_get_device_free_memory(uint32_t idDevice);
gpu_dev_arch_t 	gpu_get_device_architecture(uint32_t idDevice);
uint32_t 		gpu_get_SM_cuda_cores(gpu_dev_arch_t architecture);
uint32_t 		gpu_get_device_cuda_cores(uint32_t idDevice);
uint32_t 		gpu_get_num_devices();
uint32_t 		gpu_get_threads_per_block(gpu_dev_arch_t architecture);


/* Primitives to manage device driver options */
gpu_error_t 	gpu_set_devices_local_memory(gpu_device_info_t **devices, enum cudaFuncCache cacheConfig);
gpu_error_t 	gpu_fast_driver_awake();

/* Primitives to schedule and manage the devices */
gpu_error_t 	gpu_screen_status_device(const uint32_t idDevice, const bool deviceArchSupported,
										 const bool dataFitsMemoryDevice, const bool localReference, const bool localIndex,
										 const size_t recomendedMemorySize, const size_t requiredMemorySize);

/* Primitives to initialize device options */
gpu_error_t 	gpu_init_device(gpu_device_info_t **devices, uint32_t idDevice, uint32_t idSupportedDevice,
								const gpu_dev_arch_t selectedArchitectures);
gpu_error_t 	gpu_characterize_devices(gpu_device_info_t **devices, uint32_t numSupportedDevices);
void 			gpu_kernel_thread_configuration(const gpu_device_info_t *device, const uint32_t numThreads,
												dim3 *blocksPerGrid, dim3 *threadsPerBlock);

/* Functions to free all the buffer resources (HOST & DEVICE) */
gpu_error_t 	gpu_free_devices_list_host(gpu_device_info_t ***devices);
gpu_error_t 	gpu_free_devices_info(gpu_device_info_t **devices);

/* Collective device functions */
gpu_error_t 	gpu_reset_all_devices(gpu_device_info_t **devices);
gpu_error_t 	gpu_synchronize_all_devices(gpu_device_info_t **devices);

#endif /* GPU_DEVICES_H_ */
