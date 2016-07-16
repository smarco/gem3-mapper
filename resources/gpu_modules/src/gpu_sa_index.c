/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_SA_INDEX_C_
#define GPU_SA_INDEX_C_

#include "../include/gpu_sa_index.h"

/************************************************************
Get information functions
************************************************************/

gpu_error_t gpu_sa_index_get_size(const gpu_sa_buffer_t* const sa, size_t* const bytesPerSA)
{
  (* bytesPerSA) = sa->numEntries * sizeof(gpu_sa_entry_t);
  return (SUCCESS);
}


/************************************************************
 GLOBAL METHODS: INPUT / OUPUT Functions
************************************************************/

gpu_error_t gpu_sa_index_read_specs(int fp, gpu_sa_buffer_t* const sa)
{
  size_t result, bytesRequest;
  // Read the specifications of the SA
  bytesRequest = sizeof(uint64_t);
  result = read(fp, (void *)&sa->numEntries, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  bytesRequest = sizeof(uint64_t);
  result = read(fp, (void *)&sa->sampligRate, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_sa_index_read(int fp, gpu_sa_buffer_t* const sa)
{
  const size_t   bytesRequest = sizeof(gpu_sa_entry_t) * sa->numEntries;
  const uint64_t numRequests  = GPU_DIV_CEIL(bytesRequest, GPU_FILE_SIZE_BLOCK);
  size_t result, numBytesRequested = 0;
  uint64_t idRequest;
  // Read the SA in a block iterative way
  for(idRequest = 0; idRequest < numRequests; ++idRequest){
    const size_t requestSize = GPU_MIN(GPU_FILE_SIZE_BLOCK, bytesRequest - numBytesRequested);
    result = read(fp, (void *)sa->h_sa + numBytesRequested, requestSize);
    if (result != requestSize) return (E_READING_FILE);
    numBytesRequested += requestSize;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_sa_index_write_specs(int fp, const gpu_sa_buffer_t* const sa)
{
  size_t result, bytesRequest;
  // Write the specifications of the SA
  bytesRequest = sizeof(uint64_t);
  result = write(fp, (void *)&sa->numEntries, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  bytesRequest = sizeof(uint64_t);
  result = write(fp, (void *)&sa->sampligRate, bytesRequest);
  if (result != bytesRequest) return (E_WRITING_FILE);
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_sa_index_write(int fp, const gpu_sa_buffer_t* const sa)
{
  const size_t   bytesRequest = sizeof(gpu_sa_entry_t) * sa->numEntries;
  const uint64_t numRequests  = GPU_DIV_CEIL(bytesRequest, GPU_FILE_SIZE_BLOCK);
  size_t result, numBytesRequested = 0;
  uint64_t idRequest;
  // Write the SA in a block iterative way
  for(idRequest = 0; idRequest < numRequests; ++idRequest){
    const size_t requestSize = GPU_MIN(GPU_FILE_SIZE_BLOCK, bytesRequest - numBytesRequested);
    result = write(fp, (void *)sa->h_sa + numBytesRequested, requestSize);
    if (result != requestSize) return (E_WRITING_FILE);
    numBytesRequested += requestSize;
  }
  // Succeed
  return (SUCCESS);
}

/************************************************************
 GLOBAL METHODS: Transfer the index (HOST <-> DEVICES)
************************************************************/

gpu_error_t gpu_sa_index_transfer_CPU_to_GPUs(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices)
{
  uint32_t deviceFreeMemory, idSupportedDevice;
  uint32_t numSupportedDevices = devices[0]->numSupportedDevices;

  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(sa->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
      const size_t cpySize = sa->numEntries * sizeof(gpu_sa_entry_t);
      deviceFreeMemory = gpu_device_get_free_memory(devices[idSupportedDevice]->idDevice);
      if ((GPU_CONVERT__B_TO_MB(cpySize)) > deviceFreeMemory) return(E_INSUFFICIENT_MEM_GPU);
      CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
      //Synchronous allocate & transfer the FM-index to the GPU
      CUDA_ERROR(cudaMalloc((void**) &sa->d_sa[idSupportedDevice], cpySize));
      CUDA_ERROR(cudaMemcpy(sa->d_sa[idSupportedDevice], sa->h_sa, cpySize, cudaMemcpyHostToDevice));
    }else{
      sa->d_sa[idSupportedDevice] = sa->h_sa;
    }
  }

  return (SUCCESS);
}


/************************************************************
 GLOBAL METHODS: Functions to transform the index
************************************************************/

gpu_error_t gpu_sa_index_transform_ASCII(const char* const textBWT, gpu_sa_buffer_t* const sa)
{
  return(E_NOT_IMPLEMENTED);
}

gpu_error_t gpu_sa_index_transform_GEM_FULL(const gpu_gem_sa_dto_t* const gpu_gem_sa_dto, gpu_sa_buffer_t* const sa)
{
  uint64_t idEntry = 0;

  sa->sampligRate     = gpu_gem_sa_dto->sa_sampling;
  sa->numEntries      = GPU_DIV_CEIL(gpu_gem_sa_dto->sa_length, sa->sampligRate);
  sa->memorySpace     = NULL;
  sa->d_sa            = NULL;

  for(idEntry = 0; idEntry < sa->numEntries; ++idEntry)
    sa->h_sa[idEntry] = gpu_gem_sa_dto->sa[idEntry];

  return (SUCCESS);
}

gpu_error_t gpu_sa_index_load_specs_MFASTA_FULL(const char* const indexRaw, gpu_sa_buffer_t* const sa)
{
  return (E_NOT_IMPLEMENTED);
}

gpu_error_t gpu_sa_index_transform_MFASTA_FULL(const char* const indexRaw, gpu_sa_buffer_t* const sa)
{
  return (E_NOT_IMPLEMENTED);
}


/************************************************************
 GLOBAL METHODS: Index initialization functions
************************************************************/

gpu_error_t gpu_sa_index_init_dto(gpu_sa_buffer_t* const sa)
{
  //Initialize the SA index structure
  sa->d_sa           = NULL;
  sa->h_sa           = NULL;
  sa->hostAllocStats = GPU_PAGE_UNLOCKED;
  sa->memorySpace    = NULL;
  sa->numEntries     = 0;
  sa->sampligRate    = 0;

  return (SUCCESS);
}

gpu_error_t gpu_sa_index_init(gpu_sa_buffer_t* const sa, const uint64_t saNumEntries,
                              const uint32_t sampligRate, const uint32_t numSupportedDevices)
{
  uint32_t idSupDevice;
  GPU_ERROR(gpu_sa_index_init_dto(sa));

  sa->sampligRate    = sampligRate;
  sa->numEntries     = saNumEntries;

  sa->d_sa = (gpu_sa_entry_t **) malloc(numSupportedDevices * sizeof(gpu_sa_entry_t *));
  if (sa->d_sa == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  sa->memorySpace = (memory_alloc_t *) malloc(numSupportedDevices * sizeof(memory_alloc_t));
  if (sa->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);

  for(idSupDevice = 0; idSupDevice < numSupportedDevices; ++idSupDevice){
    sa->d_sa[idSupDevice]        = NULL;
    sa->memorySpace[idSupDevice] = GPU_NONE_MAPPED;
  }

  return (SUCCESS);
}

gpu_error_t gpu_sa_index_allocate(gpu_sa_buffer_t* const sa)
{
  if(sa->hostAllocStats & GPU_PAGE_LOCKED){
    CUDA_ERROR(cudaHostAlloc((void**) &sa->h_sa, sa->numEntries * sizeof(gpu_sa_entry_t), cudaHostAllocMapped));
  }else{
    sa->h_sa = malloc (sa->numEntries * sizeof(gpu_sa_entry_t));
    if (sa->h_sa == NULL) return (E_ALLOCATE_MEM);
  }

  return(SUCCESS);
}

/************************************************************
 GLOBAL METHODS: Functions to release DEVICE & HOST indexes
************************************************************/

gpu_error_t gpu_sa_index_free_host(gpu_sa_buffer_t* const sa)
{
    if(sa->h_sa != NULL){
      if(sa->hostAllocStats == GPU_PAGE_LOCKED) CUDA_ERROR(cudaFreeHost(sa->h_sa));
      else free(sa->h_sa);
      sa->h_sa = NULL;
    }

    return(SUCCESS);
}

gpu_error_t gpu_sa_index_free_unused_host(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices)
{
  uint32_t idSupportedDevice, numSupportedDevices;
  bool indexInHostSideUsed = false;

  numSupportedDevices = devices[0]->numSupportedDevices;
  //Free all the unused references in the host side
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(sa->memorySpace[idSupportedDevice] == GPU_HOST_MAPPED) indexInHostSideUsed = true;
  }

  if(!indexInHostSideUsed){
    GPU_ERROR(gpu_sa_index_free_host(sa));
  }

  return(SUCCESS);
}

gpu_error_t gpu_sa_index_free_device(gpu_sa_buffer_t* const sa, gpu_device_info_t** const devices)
{
  const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
  uint32_t idSupportedDevice;

  //Free all the references in the devices
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
    if(sa->d_sa[idSupportedDevice] != NULL){
      if(sa->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED)
        CUDA_ERROR(cudaFree(sa->d_sa[idSupportedDevice]));
      sa->d_sa[idSupportedDevice] = NULL;
    }
  }

  //Free the index list
  if(sa->d_sa != NULL){
    free(sa->d_sa);
    sa->d_sa = NULL;
  }

  return(SUCCESS);
}

gpu_error_t gpu_sa_index_free_metainfo(gpu_sa_buffer_t* const sa)
{
  //Free the index list
  if(sa->memorySpace != NULL){
    free(sa->memorySpace);
    sa->memorySpace = NULL;
  }
  return(SUCCESS);
}

#endif /* GPU_SA_INDEX_C_ */
