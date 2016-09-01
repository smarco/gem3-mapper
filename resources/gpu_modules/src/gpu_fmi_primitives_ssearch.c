/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_PRIMITIVES_C_
#define GPU_FMI_PRIMITIVES_C_

#include "../include/gpu_fmi_primitives.h"
#include "../include/gpu_sa_primitives.h"

/************************************************************
Functions to get the GPU FMI buffers
************************************************************/

gpu_fmi_search_seed_t* gpu_fmi_ssearch_buffer_get_seeds_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.ssearch.seeds.h_seeds);
}

gpu_sa_search_inter_t* gpu_fmi_ssearch_buffer_get_sa_intervals_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.ssearch.saIntervals.h_intervals);
}

/************************************************************
Functions to get the maximum elements of the buffers
************************************************************/

uint32_t gpu_fmi_ssearch_buffer_get_max_seeds_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.ssearch.numMaxSeeds);
}

/************************************************************
Functions to initialize the buffers (E. SEARCH)
************************************************************/

size_t gpu_fmi_ssearch_input_size()
{
  return(sizeof(gpu_fmi_search_seed_t) + sizeof(gpu_sa_search_inter_t));
}

void gpu_fmi_ssearch_reallocate_host_buffer_layout(gpu_buffer_t* mBuff)
{
  const void* rawAlloc = mBuff->h_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.ssearch.seeds.h_seeds = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.ssearch.seeds.h_seeds + mBuff->data.ssearch.numMaxSeeds);
  //Adjust the host buffer layout (output)
  mBuff->data.ssearch.saIntervals.h_intervals = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.ssearch.saIntervals.h_intervals + mBuff->data.ssearch.numMaxIntervals);
}

void gpu_fmi_ssearch_reallocate_device_buffer_layout(gpu_buffer_t* mBuff)
{
  const void* rawAlloc = mBuff->d_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.ssearch.seeds.d_seeds = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.ssearch.seeds.d_seeds + mBuff->data.ssearch.numMaxSeeds);
  //Adjust the host buffer layout (output)
  mBuff->data.ssearch.saIntervals.d_intervals = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.ssearch.saIntervals.d_intervals + mBuff->data.ssearch.numMaxIntervals);
}

void gpu_fmi_ssearch_init_buffer_(void* const fmiBuffer)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
  const size_t        sizeBuff   = mBuff->sizeBuffer * 0.95;
  const uint32_t      numInputs  = sizeBuff / gpu_fmi_ssearch_input_size();

  //set the type of the buffer
  mBuff->typeBuffer = GPU_FMI_EXACT_SEARCH;

  //set real size of the input
  mBuff->data.ssearch.numMaxSeeds     = numInputs;
  mBuff->data.ssearch.numMaxIntervals = numInputs;
  gpu_fmi_ssearch_reallocate_host_buffer_layout(mBuff);
  gpu_fmi_ssearch_reallocate_device_buffer_layout(mBuff);
}

void gpu_fmi_ssearch_init_and_realloc_buffer_(void *fmiBuffer, const uint32_t numSeeds)
{
  gpu_buffer_t* const mBuff                   = (gpu_buffer_t *) fmiBuffer;
  const uint32_t      idSupDevice             = mBuff->idSupportedDevice;
  const float         resizeFactor            = 2.0;
  const size_t        bytesPerSearchBuffer    = numSeeds * gpu_fmi_ssearch_input_size();

  //Recalculate the minimum buffer size
  mBuff->sizeBuffer = bytesPerSearchBuffer * resizeFactor;

  //FREE HOST AND DEVICE BUFFER
  GPU_ERROR(gpu_buffer_free(mBuff));

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  //ALLOCATE HOST AND DEVICE BUFFER
  CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
  CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));

  gpu_fmi_ssearch_init_buffer_(fmiBuffer);
}


/************************************************************
Debug functions for the index
************************************************************/

char gpu_fmi_ssearch_bin_to_char(const uint32_t base)
{
  switch(base){
    case GPU_ENC_DNA_CHAR_A:
      return('A');
    case GPU_ENC_DNA_CHAR_C:
      return('C');
    case GPU_ENC_DNA_CHAR_G:
      return('G');
    case GPU_ENC_DNA_CHAR_T:
      return('T');
    default :
      return('X');
  }
}

uint32_t gpu_fmi_ssearch_print_seed(const gpu_fmi_search_seed_t seed, const uint32_t seedSize)
{
  char plainSeed[GPU_FMI_SEED_MAX_CHARS] = {0};
  uint64_t bitmap = seed.hi;
  uint32_t idBase;

  for(idBase = 0; idBase < seedSize; ++idBase){
    uint32_t base = bitmap & 0x3;
    plainSeed[idBase] = gpu_fmi_ssearch_bin_to_char(base);
    bitmap >>= GPU_FMI_SEED_CHAR_LENGTH;
    if(idBase == GPU_UINT32_LENGTH) bitmap = seed.low;
  }

  for(idBase = 0; idBase < seedSize; ++idBase)
    printf("%c", plainSeed[seedSize - idBase - 1]);

  return(SUCCESS);
}

uint32_t flag_print = 0;
uint32_t gpu_fmi_ssearch_print_buffer(const void* const fmiBuffer)
{
  if(flag_print == 0){
    gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
    const uint32_t numSeeds    = mBuff->data.ssearch.seeds.numSeeds;
          uint32_t missMatches = 0;
          uint32_t idSeed;

    printf("Buffer: %d ------------------------------------\n", mBuff->idBuffer);
    for(idSeed = 0; idSeed < numSeeds; ++idSeed){
      const uint64_t hiSeedSection = mBuff->data.ssearch.seeds.h_seeds[idSeed].low;
      const uint32_t seedSize = hiSeedSection >> (GPU_UINT64_LENGTH - GPU_FMI_SEED_FIELD_SIZE);
      printf("[%d] seed=", idSeed);
      gpu_fmi_ssearch_print_seed(mBuff->data.ssearch.seeds.h_seeds[idSeed], seedSize);
      printf("\t size=%d \t (GPU) lo=%lu \t hi=%lu \n",
             seedSize,
             mBuff->data.ssearch.saIntervals.h_intervals[idSeed].low,
             mBuff->data.ssearch.saIntervals.h_intervals[idSeed].hi);
      missMatches++;
    }
    printf("Buffer: %d ------------------------------------\n", mBuff->idBuffer);
    flag_print = 1;
  }
  return (SUCCESS);
}

/************************************************************
Functions to transfer data HOST <-> DEVICE (E. SEARCH)
************************************************************/

gpu_error_t gpu_fmi_ssearch_transfer_CPU_to_GPU(gpu_buffer_t *mBuff)
{
  const gpu_fmi_search_seeds_buffer_t*  seedBuff = &mBuff->data.ssearch.seeds;
  const cudaStream_t                    idStream =  mBuff->listStreams[mBuff->idStream];
  const size_t                          cpySize  =  seedBuff->numSeeds * sizeof(gpu_fmi_search_seed_t);

  //Transfer seeds from CPU to the GPU
  CUDA_ERROR(cudaMemcpyAsync(seedBuff->d_seeds, seedBuff->h_seeds, cpySize, cudaMemcpyHostToDevice, idStream));

  return (SUCCESS);
}

gpu_error_t gpu_fmi_ssearch_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff)
{
  const gpu_fmi_search_sa_inter_buffer_t* interBuff = &mBuff->data.ssearch.saIntervals;
  const cudaStream_t                      idStream  =  mBuff->listStreams[mBuff->idStream];
  const size_t                            cpySize   =  interBuff->numIntervals * sizeof(gpu_sa_search_inter_t);

  //Transfer SA intervals (occurrence results) from CPU to the GPU
  CUDA_ERROR(cudaMemcpyAsync(interBuff->h_intervals, interBuff->d_intervals, cpySize, cudaMemcpyDeviceToHost, idStream));

  return (SUCCESS);
}

void gpu_fmi_ssearch_send_buffer_(void* const fmiBuffer, const uint32_t numSeeds)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;

  //Set real size of the input
  mBuff->data.ssearch.seeds.numSeeds = numSeeds;
  mBuff->data.ssearch.saIntervals.numIntervals = numSeeds;

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
  GPU_ERROR(gpu_fmi_ssearch_transfer_CPU_to_GPU(mBuff));
  GPU_ERROR(gpu_fmi_ssearch_process_buffer(mBuff));
  GPU_ERROR(gpu_fmi_ssearch_transfer_GPU_to_CPU(mBuff));
}

void gpu_fmi_ssearch_receive_buffer_(const void* const fmiBuffer)
{
  const gpu_buffer_t* const mBuff    = (gpu_buffer_t *) fmiBuffer;
  const cudaStream_t        idStream =  mBuff->listStreams[mBuff->idStream];

  //Synchronize Stream (the thread wait for the commands done in the stream)
  CUDA_ERROR(cudaStreamSynchronize(idStream));
}

#endif /* GPU_FMI_PRIMITIVES_C_ */

