/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_DEVICES_H_
#define GPU_DEVICES_H_

#include "gpu_commons.h"

/*************************************
GPU Interface Objects
**************************************/
/* Defines related to GPU configurations*/
#define GPU_DEVICE_STREAM_CONFIG      GPU_STREAM_BUFFER_MAPPED
#define GPU_DEVICE_STREAM_DEFAULT     0

/* Defines related to GPU Architecture */
#define GPU_WARP_SIZE                 32
/* Defines related to GPU Architecture */
#define GPU_THREADS_PER_BLOCK_FERMI   256
#define GPU_THREADS_PER_BLOCK_KEPLER  128
#define GPU_THREADS_PER_BLOCK_MAXWELL 64
#define GPU_THREADS_PER_BLOCK_PASCAL  64
#define GPU_THREADS_PER_BLOCK_VOLTA   64
#define GPU_THREADS_PER_BLOCK_NEWGEN  64

typedef enum
{
  GPU_HOST_MAPPED,
  GPU_DEVICE_MAPPED,
  GPU_NONE_MAPPED
} memory_alloc_t;

typedef enum
{
  /* Types of page-lock allocations */
  GPU_PAGE_LOCKED_PORTABLE       = GPU_UINT32_ONE_MASK << 0,
  GPU_PAGE_LOCKED_MAPPED         = GPU_UINT32_ONE_MASK << 1,
  GPU_PAGE_LOCKED_WRITECOMBINED  = GPU_UINT32_ONE_MASK << 2,
  /* Types of host allocations */
  GPU_PAGE_UNLOCKED              = GPU_UINT32_ONE_MASK << 3,
  GPU_PAGE_LOCKED                = GPU_PAGE_LOCKED_PORTABLE | GPU_PAGE_LOCKED_MAPPED | GPU_PAGE_LOCKED_WRITECOMBINED,
  /* Types for non-allocated pages */
  GPU_PAGE_UNALLOCATED           = 0,
  /* Types for assigned pages (contains data) */
  GPU_PAGE_ASSIGNED              = GPU_UINT32_ONE_MASK << 31,
  /*  */
  GPU_PAGE_ASSIGNED_AND_UNLOCKED = GPU_PAGE_ASSIGNED | GPU_PAGE_UNLOCKED,
  GPU_PAGE_ASSIGNED_AND_LOCKED   = GPU_PAGE_ASSIGNED | GPU_PAGE_LOCKED
} memory_stats_t;

typedef enum
{
  GPU_STREAM_BUFFER_MAPPED,
  GPU_STREAM_THREAD_MAPPED,
  GPU_STREAM_APPLICATION_MAPPED
} stream_config_t;

typedef struct {
  /* System specifications */
  uint32_t        numDevices;
  uint32_t        numSupportedDevices;
  /* Device specifications */
  uint32_t        idDevice;
  uint32_t        idSupportedDevice;
  gpu_dev_arch_t  architecture;
  uint32_t        cudaCores;
  float           coreClockRate;        // Ghz
  uint32_t        memoryBusWidth;       // Bits
  float           memoryClockRate;      // Ghz
  /* Device performance metrics */
  float           absolutePerformance;  // GOps/s
  float           relativePerformance;  // Ratio
  float           absoluteBandwidth;    // GB/s
  float           relativeBandwidth;    // Ratio
  /* System performance metrics */
  float           allSystemPerformance; // GOps/s
  float           allSystemBandwidth;   // GB/s
} gpu_device_info_t;

/* Primitives to get information for the scheduler */
size_t          gpu_device_get_free_memory(const uint32_t idDevice);
gpu_dev_arch_t  gpu_device_get_architecture(const uint32_t idDevice);
uint32_t        gpu_device_get_SM_cuda_cores(const gpu_dev_arch_t architecture);
uint32_t        gpu_device_get_cuda_cores(const uint32_t idDevice);
uint32_t        gpu_device_get_num_all();
uint32_t        gpu_device_get_threads_per_block(const gpu_dev_arch_t architecture);
uint32_t        gpu_device_get_stream_configuration(const stream_config_t streamConfig, const uint64_t idThread, const uint32_t idBuffer);

/* Primitives to manage device driver options */
gpu_error_t     gpu_device_set_local_memory_all(gpu_device_info_t **devices, enum cudaFuncCache cacheConfig);
gpu_error_t     gpu_device_fast_driver_awake();

/* Primitives to schedule and manage the devices */
gpu_error_t     gpu_device_setup_system(gpu_device_info_t **devices);
gpu_error_t     gpu_device_screen_status(const uint32_t idDevice, const bool deviceArchSupported, const size_t recomendedMemorySize, const size_t requiredMemorySize);

/* Primitives to initialize device options */
gpu_error_t     gpu_device_init(gpu_device_info_t **devices, uint32_t idDevice, uint32_t idSupportedDevice, const gpu_dev_arch_t selectedArchitectures);
gpu_error_t     gpu_device_characterize_all(gpu_device_info_t **devices, uint32_t numSupportedDevices);
void            gpu_device_kernel_thread_configuration(const gpu_device_info_t *device, const uint32_t numThreads, dim3* const blocksPerGrid, dim3* const threadsPerBlock);

/* Functions to free all the buffer resources (HOST & DEVICE) */
gpu_error_t     gpu_device_free_list(gpu_device_info_t ***devices);
gpu_error_t     gpu_device_free_info_all(gpu_device_info_t **devices);

/* Collective device functions */
gpu_error_t     gpu_device_reset_all(gpu_device_info_t **devices);
gpu_error_t     gpu_device_synchronize(gpu_device_info_t **devices, uint32_t idSupDevice, cudaStream_t idStream);
gpu_error_t     gpu_device_synchronize_all(gpu_device_info_t **devices);

#endif /* GPU_DEVICES_H_ */
