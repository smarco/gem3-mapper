/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: gpu_scheduler.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Common headers and data structures for BPM on GPU library
 */

#ifndef GPU_SCHEDULER_H_
#define GPU_SCHEDULER_H_

#include "gpu_commons.h"
#include "gpu_bpm.h"
#include "gpu_fmi.h"

/* Defines related to GPU Architecture */
#define	GPU_WARP_SIZE					32
#ifndef GPU_MAX_THREADS_PER_BLOCK
	#define	GPU_MAX_THREADS_PER_BLOCK	128
#endif
#ifndef GPU_MAX_THREADS_PER_SM
	#define	GPU_MAX_THREADS_PER_SM		128
#endif


/*************************************
GPU Interface Objects
**************************************/

typedef enum
{
	GPU_HOST_MAPPED,
	GPU_DEVICE_MAPPED
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
	gpu_device_info_t		*device;
	gpu_reference_buffer_t	*reference;
	gpu_index_buffer_t		*index;
	size_t					sizeBuffer;
	void 					*h_rawData;
	void 					*d_rawData;
	gpu_data_buffer_t		data;
} gpu_buffer_t;

#endif /* GPU_SCHEDULER_H_ */
