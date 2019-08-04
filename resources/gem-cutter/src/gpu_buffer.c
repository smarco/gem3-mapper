/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include "../include/gpu_buffer.h"

#ifndef GPU_BUFFER_C_
#define GPU_BUFFER_C_

/************************************************************
Functions to get information from the system
************************************************************/

uint32_t gpu_buffer_get_id_device_(const void* const gpuBuffer)
{
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) gpuBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;
  return(mBuff->device[idSupDevice]->idDevice);
}

uint32_t gpu_buffer_get_id_supported_device_(const void* const gpuBuffer)
{
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) gpuBuffer;
  return(mBuff->idSupportedDevice);
}

gpu_error_t gpu_buffer_free(gpu_buffer_t *mBuff)
{
  if(mBuff->h_rawData != NULL){
    CUDA_ERROR(cudaFreeHost(mBuff->h_rawData));
    mBuff->h_rawData = NULL;
  }

  if(mBuff->d_rawData != NULL){
    CUDA_ERROR(cudaFree(mBuff->d_rawData));
    mBuff->d_rawData = NULL;
  }

  return(SUCCESS);
}

gpu_error_t gpu_buffer_get_min_memory_size(size_t *bytesPerBuffer)
{
  const uint32_t averarageNumPEQEntries = 1;
  const uint32_t candidatesPerQuery     = 1;
  const uint32_t averageQuerySize       = 100;
  const uint32_t averageRegionsPerQuery = 10;

  const size_t bytesPerBPMBuffer     = GPU_BPM_FILTER_MIN_ELEMENTS * gpu_bpm_filter_size_per_candidate(averarageNumPEQEntries,candidatesPerQuery);
  const size_t bytesPerSSearchBuffer = GPU_FMI_SEARCH_MIN_ELEMENTS * gpu_fmi_asearch_size_per_query(averageQuerySize, averageRegionsPerQuery);
  const size_t bytesPerASearchBuffer = GPU_FMI_SEARCH_MIN_ELEMENTS * gpu_fmi_ssearch_input_size();
  const size_t bytesPerDecodeBuffer  = GPU_FMI_DECODE_MIN_ELEMENTS * gpu_fmi_decode_input_size();

  (* bytesPerBuffer) = GPU_MAX(bytesPerBPMBuffer,GPU_MAX(bytesPerDecodeBuffer,GPU_MAX(bytesPerSSearchBuffer,bytesPerASearchBuffer)));
  return (SUCCESS);
}

gpu_error_t gpu_buffer_scheduling(gpu_buffer_t ***gpuBuffer, const uint32_t numBuffers, gpu_device_info_t** const device,
                                  gpu_reference_buffer_t *reference, gpu_index_buffer_t *index, float maxMbPerBuffer)
{
  uint32_t idSupportedDevice, numBuffersPerDevice, idLocalBuffer;
  const size_t maxBytesPerBuffer = GPU_CONVERT_MB_TO__B(maxMbPerBuffer);
  const uint32_t numSupportedDevices = device[0]->numSupportedDevices;
  int32_t remainderBuffers = numBuffers, idGlobalBuffer = 0;

  gpu_buffer_t **buffer = (gpu_buffer_t **) malloc(numBuffers * sizeof(gpu_buffer_t *));
  if (buffer == NULL) GPU_ERROR(E_ALLOCATE_MEM);

  cudaStream_t* const listStreams = (cudaStream_t *) malloc(numBuffers * sizeof(cudaStream_t));
  if (listStreams == NULL) GPU_ERROR(E_ALLOCATE_MEM);

  /* Assigning buffers for each GPU (to adapt the workload) */
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    size_t bytesPerDevice, bytesPerBuffer, minimumMemorySize, minBytesPerBuffer;
    const uint32_t idDevice = device[idSupportedDevice]->idDevice;
    const size_t freeDeviceMemory = gpu_device_get_free_memory(idDevice);

    numBuffersPerDevice = GPU_ROUND(numBuffers * device[idSupportedDevice]->relativePerformance);
    if(idSupportedDevice == numSupportedDevices-1) numBuffersPerDevice = remainderBuffers;

    /* Resize the buffers to the available memory (Best effort system) */
    if(maxBytesPerBuffer != 0) bytesPerDevice = GPU_MIN(numBuffersPerDevice * maxBytesPerBuffer, freeDeviceMemory);
      else bytesPerDevice = freeDeviceMemory;

    bytesPerBuffer = bytesPerDevice / numBuffersPerDevice;

    GPU_ERROR(gpu_module_get_min_memory(reference, index, numBuffers, GPU_NONE_MODULES, &minimumMemorySize));
    if(freeDeviceMemory < minimumMemorySize) GPU_ERROR(E_INSUFFICIENT_MEM_GPU);
    GPU_ERROR(gpu_buffer_get_min_memory_size(&minBytesPerBuffer));
    if(minBytesPerBuffer > bytesPerBuffer) GPU_ERROR(E_INSUFFICIENT_MEM_PER_BUFFER);

    for(idLocalBuffer = 0; idLocalBuffer < numBuffersPerDevice; ++idLocalBuffer){
      buffer[idGlobalBuffer] = (gpu_buffer_t *) malloc(sizeof(gpu_buffer_t));
      CUDA_ERROR(cudaSetDevice(idDevice));
      GPU_ERROR(gpu_buffer_configuration(buffer[idGlobalBuffer],idGlobalBuffer,idSupportedDevice,bytesPerBuffer,numBuffers,device,listStreams,reference,index));
      idGlobalBuffer++;
    }
    remainderBuffers -= numBuffersPerDevice;
  }

  (* gpuBuffer) = buffer;
  return (SUCCESS);
}

gpu_error_t gpu_buffer_configuration(gpu_buffer_t* const mBuff, const uint32_t idBuffer, const uint32_t idSupportedDevice,
                                     const size_t bytesPerBuffer, const uint32_t numBuffers, gpu_device_info_t** const device,
                                     cudaStream_t* const listStreams, gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index)
{
  /* Buffer information */
  mBuff->idBuffer           = idBuffer;
  mBuff->numBuffers         = numBuffers;
  mBuff->idSupportedDevice  = idSupportedDevice;
  mBuff->device             = device;
  mBuff->sizeBuffer         = bytesPerBuffer;
  mBuff->typeBuffer         = GPU_NONE_MODULES;
  mBuff->listStreams        = listStreams;
  /* Module structures */
  mBuff->index              = index;
  mBuff->reference          = reference;
  /* Chunk of RAW memory for the buffer */
  mBuff->h_rawData          = NULL;
  mBuff->d_rawData          = NULL;
  /* Set in which Device we create and initialize the structures */
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupportedDevice]->idDevice));
  /* Create the CUDA stream per each buffer */
  CUDA_ERROR(cudaStreamCreate(&listStreams[idBuffer]));
  /* Suceed */
  return(SUCCESS);
}

void gpu_init_buffers_(gpu_buffers_dto_t* const buff, gpu_index_dto_t* const rawIndex,
                       gpu_reference_dto_t* const rawRef, gpu_info_dto_t* const sys)
{
  /* Buffer info */
  const float               maxMbPerBuffer        = buff->maxMbPerBuffer;
  const uint32_t            numBuffers            = buff->numBuffers;
  const gpu_module_t        activeModules         = buff->activeModules;
  /* System info */
  const gpu_dev_arch_t      selectedArchitectures = sys->selectedArchitectures;
  const gpu_data_location_t userAllocOption       = sys->userAllocOption;
  const uint32_t            numSupportedDevices   = gpu_get_num_supported_devices_(selectedArchitectures);
  /* Internal buffers info */
  gpu_buffer_t              **buffer              = NULL;
  gpu_reference_buffer_t    *reference            = NULL;
  gpu_index_buffer_t        *index                = NULL;
  gpu_device_info_t         **devices             = NULL;

  /* Reference and index initialization */
  GPU_ERROR(gpu_reference_init(&reference, rawRef, numSupportedDevices, activeModules));
  GPU_ERROR(gpu_index_init(&index, rawIndex, numSupportedDevices, activeModules));

  /* Analyze and search for the best module configuration (modules auto-activation) */
  GPU_ERROR(gpu_module_configure_system(reference, index, &devices, numBuffers, selectedArchitectures, userAllocOption,
                                        &sys->activatedModules, &sys->allocatedStructures));
  GPU_ERROR(gpu_device_setup_system(devices));

  /* Allocate, transform and send reference to corresponding DEVICES */
  GPU_ERROR(gpu_reference_load(reference, rawRef, reference->activeModules));
  GPU_ERROR(gpu_reference_transfer_CPU_to_GPUs(reference, devices, reference->activeModules));

  /* Allocate, transform and send index to corresponding DEVICES */
  GPU_ERROR(gpu_index_load(index, rawIndex, index->activeModules));
  GPU_ERROR(gpu_index_transfer_CPU_to_GPUs(index, devices, index->activeModules));

  /* Characterize all the system, create the buffers and balance the work along all DEVICES */
  GPU_ERROR(gpu_buffer_scheduling(&buffer, numBuffers, devices, reference, index, maxMbPerBuffer));

  /* Freeing unused host structures */
  GPU_ERROR(gpu_reference_free_unused_host(reference, devices, reference->activeModules));
  GPU_ERROR(gpu_index_free_unused_host(index, devices, index->activeModules));

  buff->buffer = (void **) buffer;
}

void gpu_destroy_buffers_(gpu_buffers_dto_t* buff)
{
  gpu_buffer_t** mBuff        = (gpu_buffer_t **) buff->buffer;
  gpu_device_info_t** devices = mBuff[0]->device;
  const uint32_t numBuffers   = mBuff[0]->numBuffers;
  uint32_t idBuffer;

  GPU_ERROR(gpu_device_synchronize_all(devices));

  /* Free all the references */
  GPU_ERROR(gpu_reference_free(&mBuff[0]->reference, devices, mBuff[0]->reference->activeModules));
  GPU_ERROR(gpu_index_free(&mBuff[0]->index, devices, mBuff[0]->index->activeModules));

  for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
    const uint32_t idSupDevice = mBuff[idBuffer]->idSupportedDevice;
    CUDA_ERROR(cudaSetDevice(devices[idSupDevice]->idDevice));
    GPU_ERROR(gpu_buffer_free(mBuff[idBuffer]));
    CUDA_ERROR(cudaStreamDestroy(mBuff[idBuffer]->listStreams[idBuffer]));
  }

  /* Deallocate the global streams */
  if(mBuff[0]->listStreams != NULL){
    free(mBuff[0]->listStreams);
    mBuff[0]->listStreams = NULL;
  }

  /* Deallocate each buffer */
  for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
    free(mBuff[idBuffer]);
  }

  GPU_ERROR(gpu_device_reset_all(devices));
  GPU_ERROR(gpu_device_free_info_all(devices));

  if(mBuff != NULL){
    free(mBuff);
    mBuff = NULL;
    buff->buffer = NULL;
  }
}

/************************************************************
Primitives to schedule and manage the buffers
************************************************************/

void gpu_alloc_buffer_(void* const gpuBuffer, const uint64_t idThread)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) gpuBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;
  const uint32_t idBuffer    = mBuff->idBuffer;
  const uint32_t idStream    = gpu_device_get_stream_configuration(GPU_DEVICE_STREAM_CONFIG, idThread - 1, idBuffer);

  mBuff->h_rawData = NULL;
  mBuff->d_rawData = NULL;
  mBuff->idStream  = idStream;

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  //ALLOCATE HOST AND DEVICE BUFFER
  CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
  CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));
}

void gpu_realloc_buffer_(void* const gpuBuffer, const float maxMbPerBuffer)
{
  gpu_buffer_t* const mBuff = (gpu_buffer_t *) gpuBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;
  mBuff->sizeBuffer = GPU_CONVERT_MB_TO__B(maxMbPerBuffer);

  //FREE HOST AND DEVICE BUFFER
  GPU_ERROR(gpu_buffer_free(mBuff));

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  //ALLOCATE HOST AND DEVICE BUFFER
  CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
  CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));
}

#endif /* GPU_BUFFER_C_ */

