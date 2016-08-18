/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_DEVICES_C_
#define GPU_DEVICES_C_

#include "../include/gpu_devices.h"

/************************************************************
Primitives to get devices properties
************************************************************/

uint32_t gpu_device_get_threads_per_block(const gpu_dev_arch_t architecture)
{
  if(architecture & GPU_ARCH_FERMI)   return(GPU_THREADS_PER_BLOCK_FERMI);
  if(architecture & GPU_ARCH_KEPLER)  return(GPU_THREADS_PER_BLOCK_KEPLER);
  if(architecture & GPU_ARCH_MAXWELL) return(GPU_THREADS_PER_BLOCK_MAXWELL);
  return(GPU_THREADS_PER_BLOCK_NEWGEN);
}

gpu_dev_arch_t gpu_device_get_architecture(const uint32_t idDevice)
{
  struct cudaDeviceProp devProp;
  CUDA_ERROR(cudaGetDeviceProperties(&devProp, idDevice));
                                                                              /* {Y:Mayor X:Minor} */
  if (devProp.major <= 1) return(GPU_ARCH_TESLA);                             /* CC 1.X            */
  if (devProp.major == 2 && devProp.minor == 0) return(GPU_ARCH_FERMI_1G);    /* CC 2.0            */
  if (devProp.major == 2 && devProp.minor >  0) return(GPU_ARCH_FERMI_2G);    /* CC 2.1            */
  if (devProp.major == 3 && devProp.minor <  5) return(GPU_ARCH_KEPLER_1G);   /* CC 3.0, 3.2       */
  if (devProp.major == 3 && devProp.minor >= 5) return(GPU_ARCH_KEPLER_2G);   /* CC 3.5            */
  if (devProp.major == 5 && devProp.minor == 0) return(GPU_ARCH_MAXWELL_1G);  /* CC 5.0            */
  if (devProp.major == 5 && devProp.minor >= 5) return(GPU_ARCH_MAXWELL_2G);  /* CC 5.2            */
  if (devProp.major == 6 && devProp.minor == 0) return(GPU_ARCH_PASCAL_1G);   /* CC 6.0            */
  if (devProp.major == 6 && devProp.minor >= 5) return(GPU_ARCH_PASCAL_2G);   /* CC 6.2            */
  return(GPU_ARCH_NEWGEN);                                                    /* CC Y.X            */
}

uint32_t gpu_device_get_SM_cuda_cores(const gpu_dev_arch_t architecture)
{
  switch (architecture) {
    case GPU_ARCH_TESLA:      return(8);
    case GPU_ARCH_FERMI_1G:   return(32);
    case GPU_ARCH_FERMI_2G:   return(48);
    case GPU_ARCH_KEPLER_1G:  return(192);
    case GPU_ARCH_KEPLER_2G:  return(192);
    case GPU_ARCH_MAXWELL_1G: return(128);
    case GPU_ARCH_MAXWELL_2G: return(128);
    case GPU_ARCH_PASCAL_1G:  return(128);
    case GPU_ARCH_PASCAL_2G:  return(128);
    default:                  return(128);
  }
}

uint32_t gpu_device_get_cuda_cores(const uint32_t idDevice)
{
  uint32_t coresPerSM;
  gpu_dev_arch_t architecture;

  struct cudaDeviceProp devProp;
  CUDA_ERROR(cudaGetDeviceProperties(&devProp, idDevice));

  architecture = gpu_device_get_architecture(idDevice);
  coresPerSM = gpu_device_get_SM_cuda_cores(architecture);

  return(devProp.multiProcessorCount * coresPerSM);
}

size_t gpu_device_get_free_memory(const uint32_t idDevice)
{
  size_t free, total;
  CUDA_ERROR(cudaSetDevice(idDevice));
  CUDA_ERROR(cudaMemGetInfo(&free, &total));
  return (free * 0.95);
}

uint32_t gpu_device_get_num_all()
{
  int32_t numDevices;
  CUDA_ERROR(cudaGetDeviceCount(&numDevices));
  return(numDevices);
}

uint32_t gpu_get_num_supported_devices_(const gpu_dev_arch_t selectedArchitectures)
{
  uint32_t idDevice, numSupportedDevices = 0;
  int32_t numDevices;

  CUDA_ERROR(cudaGetDeviceCount(&numDevices));

  for(idDevice = 0; idDevice < numDevices; ++idDevice){
    const gpu_dev_arch_t deviceArch = gpu_device_get_architecture(idDevice);
    if(deviceArch & selectedArchitectures) numSupportedDevices++;
  }

  return(numSupportedDevices);
}

uint32_t gpu_device_get_stream_configuration(const stream_config_t streamConfig, const uint64_t idThread, const uint32_t idBuffer)
{
  switch (streamConfig){
    case GPU_STREAM_BUFFER_MAPPED:
      return(idBuffer);
    case GPU_STREAM_THREAD_MAPPED:
      return(idThread);
    case GPU_STREAM_APPLICATION_MAPPED:
      return(GPU_DEVICE_STREAM_DEFAULT);
    default: // Original GEM configuration
      return(idBuffer);
  }
}


/************************************************************
Primitives to set devices
************************************************************/

gpu_error_t gpu_device_characterize_all(gpu_device_info_t **dev, const uint32_t numSupportedDevices)
{
  uint32_t idSupportedDevice, allSystemPerformance = 0, allSystemBandwidth = 0;

  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    allSystemPerformance += dev[idSupportedDevice]->absolutePerformance;
    allSystemBandwidth += dev[idSupportedDevice]->absoluteBandwidth;
  }

  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    dev[idSupportedDevice]->numSupportedDevices   = numSupportedDevices;
    dev[idSupportedDevice]->relativeBandwidth     = dev[idSupportedDevice]->absoluteBandwidth / allSystemBandwidth;
    dev[idSupportedDevice]->allSystemBandwidth    = allSystemBandwidth;
    dev[idSupportedDevice]->relativePerformance   = dev[idSupportedDevice]->absolutePerformance / allSystemPerformance;
    dev[idSupportedDevice]->allSystemPerformance  = allSystemPerformance;
  }
  return(SUCCESS);
}

gpu_error_t gpu_device_reset_all(gpu_device_info_t **devices)
{
  const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
  uint32_t idSupportedDevice;

  /* reset all the device environments */
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
    CUDA_ERROR(cudaDeviceReset());
  }
  return(SUCCESS);
}

gpu_error_t gpu_device_synchronize_all(gpu_device_info_t **devices)
{
  const uint32_t  numSupportedDevices = devices[0]->numSupportedDevices;
  uint32_t        idSupportedDevice;

  /* Synchronize all the Devices to the Host */
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
    CUDA_ERROR(cudaDeviceSynchronize());
  }
  return(SUCCESS);
}

void gpu_device_kernel_thread_configuration(const gpu_device_info_t *device, const uint32_t numThreads,
                                            dim3* const blocksPerGrid, dim3* const threadsPerBlock)
{
  if(device->architecture & (GPU_ARCH_FERMI)){
    //We use 2-Dimensional Grid (because Fermi is limited to 65535 Blocks per dim)
    const uint32_t threadsPerRow = gpu_device_get_threads_per_block(device->architecture);
    const uint32_t maxBlocksPerRow = 65535;
    const uint32_t numBlocks = GPU_DIV_CEIL(numThreads, threadsPerRow);
    const uint32_t rowsPerGrid = GPU_DIV_CEIL(numBlocks, maxBlocksPerRow);
    const uint32_t blocksPerRow = (rowsPerGrid > 1) ? maxBlocksPerRow : numBlocks;

    blocksPerGrid->x   = blocksPerRow;
    blocksPerGrid->y   = rowsPerGrid;
    threadsPerBlock->x = threadsPerRow;
  }else{
    const uint32_t threadsPerRow = gpu_device_get_threads_per_block(device->architecture);
    threadsPerBlock->x = threadsPerRow;
    blocksPerGrid->x   = GPU_DIV_CEIL(numThreads, threadsPerRow);
  }
}


/************************************************************
Primitives to manage device driver options
************************************************************/

gpu_error_t gpu_device_set_local_memory_all(gpu_device_info_t **devices, enum cudaFuncCache cacheConfig)
{
  uint32_t idSupportedDevice, numSupportedDevices = devices[0]->numSupportedDevices;

  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
    CUDA_ERROR(cudaDeviceSetCacheConfig(cacheConfig));
  }

  return (SUCCESS);
}

gpu_error_t gpu_device_fast_driver_awake()
{
  //Dummy call to the NVIDIA API to awake earlier the driver.
  void* prt = NULL;
  cudaMalloc(&prt, 0);

  return(SUCCESS);
}


/************************************************************
Primitives to schedule and manage the devices
************************************************************/

gpu_error_t gpu_device_setup_system(gpu_device_info_t **devices)
{
  /* List here all directives to configure all devices*/
  // Fast awake of the driver
  GPU_ERROR(gpu_device_fast_driver_awake());
  // Set preferred all L1 GPU caches
  //enum cudaFuncCache cacheConfig = cudaFuncCachePreferL1;
  GPU_ERROR(gpu_device_set_local_memory_all(devices, cudaFuncCachePreferL1));

  return (SUCCESS);
}

gpu_error_t gpu_device_init(gpu_device_info_t **device, const uint32_t idDevice, const uint32_t idSupportedDevice,
                            const gpu_dev_arch_t selectedArchitectures)
{
  gpu_device_info_t *dev = NULL;
  struct cudaDeviceProp devProp;

  CUDA_ERROR(cudaGetDeviceProperties(&devProp, idDevice));

  dev = (gpu_device_info_t *) malloc(sizeof(gpu_device_info_t));
  if(dev == NULL) return(E_ALLOCATE_MEM);

  dev->numDevices           = gpu_device_get_num_all();
  dev->numSupportedDevices  = gpu_get_num_supported_devices_(selectedArchitectures);
  dev->idSupportedDevice    = idSupportedDevice;
  dev->idDevice             = idDevice;
  dev->architecture         = gpu_device_get_architecture(idDevice);
  dev->cudaCores            = gpu_device_get_cuda_cores(idDevice);
  dev->coreClockRate        = GPU_CONVERT_KHZ_TO_GHZ((float)devProp.clockRate);
  dev->memoryBusWidth       = GPU_CONVERT_BITS_TO_BYTES((float)devProp.memoryBusWidth);
  dev->memoryClockRate      = GPU_CONVERT_KHZ_TO_GHZ((float)devProp.memoryClockRate);

  dev->absolutePerformance  = dev->cudaCores * dev->coreClockRate;
  dev->absoluteBandwidth    = 2.0 * dev->memoryClockRate * dev->memoryBusWidth; // GBytes/s

  (* device) = dev;
  return(SUCCESS);
}

gpu_error_t gpu_device_screen_status(const uint32_t idDevice, const bool deviceArchSupported,
                                     const size_t recomendedMemorySize, const size_t requiredMemorySize)
{
  struct cudaDeviceProp devProp;
  const size_t memoryFree = gpu_device_get_free_memory(idDevice);
  const bool dataFitsMemoryDevice = (requiredMemorySize < memoryFree);

  CUDA_ERROR(cudaGetDeviceProperties(&devProp, idDevice));

  if(!deviceArchSupported){
    fprintf(stderr, "GPU Device %d: %-25s Discarded Device [NOT RUNNING] \n", idDevice, devProp.name);
    fprintf(stderr, "GPU Device %d: %-25s UNSUPPORTED Architecture [Compute Capability < 2.0: UNSUITABLE] \n", idDevice, devProp.name);
    return(SUCCESS);
  }
  if(!dataFitsMemoryDevice){
    fprintf(stderr, "GPU Device %d: %-25s Discarded Device [NOT RUNNING] \n", idDevice, devProp.name);
    fprintf(stderr, "GPU Device %d: %-25s INSUFFICIENT DEVICE MEMORY (Mem Required: %lu MBytes - Mem Available: %lu MBytes) \n",
            idDevice, devProp.name, GPU_CONVERT__B_TO_MB(requiredMemorySize), GPU_CONVERT__B_TO_MB(memoryFree));
    return(SUCCESS);
  }
  if((deviceArchSupported && dataFitsMemoryDevice)){
    fprintf(stderr, "GPU Device %d: %-25s Selected Device [RUNNING] \n", idDevice, devProp.name);
    if (memoryFree < recomendedMemorySize){
      fprintf(stderr, "GPU Device %d: %-25s WARNING PERF. SLOWDOWNS (Mem Recommended: %lu MBytes - Mem Available: %lu MBytes) \n",
              idDevice, devProp.name, GPU_CONVERT__B_TO_MB(recomendedMemorySize), GPU_CONVERT__B_TO_MB(memoryFree));
    }
    return(SUCCESS);
  }
  return(SUCCESS);
}


/************************************************************
Functions to free all the buffer resources (HOST & DEVICE)
************************************************************/

gpu_error_t gpu_device_free_list(gpu_device_info_t ***deviceList)
{
  gpu_device_info_t **device = (* deviceList);

  if(device != NULL){
    free(device);
    device = NULL;
  }

  (* deviceList) = device;
  return(SUCCESS);
}

gpu_error_t gpu_device_free_info_all(gpu_device_info_t **devices)
{
  uint32_t idSupportedDevice, numSupportedDevices;

  numSupportedDevices = devices[0]->numSupportedDevices;
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(devices[idSupportedDevice] != NULL){
      free(devices[idSupportedDevice]);
      devices[idSupportedDevice] = NULL;
    }
  }

  GPU_ERROR(gpu_device_free_list(&devices));
  return(SUCCESS);
}

#endif /* GPU_DEVICES_C_ */
