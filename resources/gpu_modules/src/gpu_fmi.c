/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Host scheduler for BPM on GPU
 */

#include "../include/gpu_fmi.h"

/************************************************************
Functions to get the GPU FMI buffers
************************************************************/

gpu_fmi_entry_t* gpu_fmi_buffer_get_index_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->index->h_fmi);
}

gpu_fmi_search_seed_t* gpu_fmi_search_buffer_get_seeds_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.search.seeds.h_seeds);
}

gpu_fmi_search_sa_inter_t* gpu_fmi_search_buffer_get_sa_intervals_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.search.saIntervals.h_intervals);
}

gpu_fmi_decode_init_pos_t* gpu_fmi_decode_buffer_get_init_pos_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.decode.initPositions.h_initBWTPos);
}

gpu_fmi_decode_end_pos_t* gpu_fmi_decode_buffer_get_end_pos_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.decode.endPositions.h_endBWTPos);
}


/************************************************************
Functions to get the maximum elements of the buffers
************************************************************/

uint32_t gpu_fmi_search_buffer_get_max_seeds_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.search.numMaxSeeds);
}

uint32_t gpu_fmi_decode_buffer_get_max_positions_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.decode.numMaxInitPositions);
}

/************************************************************
Functions to init the buffers (E. SEARCH)
************************************************************/

size_t gpu_fmi_search_input_size()
{
  return(sizeof(gpu_fmi_search_seed_t) + sizeof(gpu_fmi_search_sa_inter_t));
}

void gpu_fmi_search_reallocate_host_buffer_layout(gpu_buffer_t* mBuff)
{
  const void* rawAlloc = mBuff->h_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.search.seeds.h_seeds = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.search.seeds.h_seeds + mBuff->data.search.numMaxSeeds);
  //Adjust the host buffer layout (output)
  mBuff->data.search.saIntervals.h_intervals = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.search.saIntervals.h_intervals + mBuff->data.search.numMaxIntervals);
}

void gpu_fmi_search_reallocate_device_buffer_layout(gpu_buffer_t* mBuff)
{
  const void* rawAlloc = mBuff->d_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.search.seeds.d_seeds = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.search.seeds.d_seeds + mBuff->data.search.numMaxSeeds);
  //Adjust the host buffer layout (output)
  mBuff->data.search.saIntervals.d_intervals = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.search.saIntervals.d_intervals + mBuff->data.search.numMaxIntervals);
}

void gpu_fmi_search_init_buffer_(void* const fmiBuffer)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
  const size_t        sizeBuff   = mBuff->sizeBuffer * 0.95;
  const uint32_t      numInputs  = sizeBuff / gpu_fmi_search_input_size();

  //set the type of the buffer
  mBuff->typeBuffer = GPU_FMI_EXACT_SEARCH;

  //set real size of the input
  mBuff->data.search.numMaxSeeds     = numInputs;
  mBuff->data.search.numMaxIntervals = numInputs;
  gpu_fmi_search_reallocate_host_buffer_layout(mBuff);
  gpu_fmi_search_reallocate_device_buffer_layout(mBuff);
}

void gpu_fmi_search_init_and_realloc_buffer_(void *fmiBuffer, const uint32_t numSeeds)
{
  gpu_buffer_t* const mBuff                   = (gpu_buffer_t *) fmiBuffer;
  const uint32_t      idSupDevice             = mBuff->idSupportedDevice;
  const float         resizeFactor            = 2.0;
  const size_t        bytesPerSearchBuffer    = numSeeds * gpu_fmi_search_input_size();

  //Recalculate the minimum buffer size
  mBuff->sizeBuffer = bytesPerSearchBuffer * resizeFactor;

  //FREE HOST AND DEVICE BUFFER
  GPU_ERROR(gpu_free_buffer(mBuff));

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  //ALLOCATE HOST AND DEVICE BUFFER
  CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
  CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));

  gpu_fmi_search_init_buffer_(fmiBuffer);
}


/************************************************************
Functions to init the buffers (DECODE)
************************************************************/

size_t gpu_fmi_decode_input_size()
{
  return(sizeof(gpu_fmi_decode_init_pos_t) + sizeof(gpu_fmi_decode_end_pos_t));
}

void gpu_fmi_decode_reallocate_host_buffer_layout(gpu_buffer_t* const mBuff)
{
  const void* rawAlloc = mBuff->h_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.decode.initPositions.h_initBWTPos = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.initPositions.h_initBWTPos + mBuff->data.decode.numMaxInitPositions);
  //Adjust the host buffer layout (output)
  mBuff->data.decode.endPositions.h_endBWTPos   = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.endPositions.h_endBWTPos + mBuff->data.decode.numMaxEndPositions);
}

void gpu_fmi_decode_reallocate_device_buffer_layout(gpu_buffer_t* const mBuff)
{
  const void* rawAlloc = mBuff->d_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.decode.initPositions.d_initBWTPos = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.initPositions.d_initBWTPos + mBuff->data.decode.numMaxInitPositions);
  //Adjust the host buffer layout (output)
  mBuff->data.decode.endPositions.d_endBWTPos   = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.endPositions.d_endBWTPos + mBuff->data.decode.numMaxEndPositions);
}

void gpu_fmi_decode_init_buffer_(void* const fmiBuffer)
{
  gpu_buffer_t* const mBuff            = (gpu_buffer_t *) fmiBuffer;
  const size_t        sizeBuff         = mBuff->sizeBuffer * 0.95;
  const uint32_t      numMaxPositions  = sizeBuff / gpu_fmi_decode_input_size();

  //set the type of the buffer
  mBuff->typeBuffer = GPU_FMI_DECODE_POS;

  //Set real size of the input
  mBuff->data.decode.numMaxInitPositions = numMaxPositions;
  mBuff->data.decode.numMaxEndPositions = numMaxPositions;
  gpu_fmi_decode_reallocate_host_buffer_layout(mBuff);
  gpu_fmi_decode_reallocate_device_buffer_layout(mBuff);
}

void gpu_fmi_decode_init_and_realloc_buffer_(void* const fmiBuffer, const uint32_t numDecodes)
{
  gpu_buffer_t* const mBuff                   = (gpu_buffer_t *) fmiBuffer;
  const uint32_t      idSupDevice             = mBuff->idSupportedDevice;
  const float         resizeFactor            = 2.0;
  const size_t        bytesPerDecodeBuffer    = numDecodes * gpu_fmi_decode_input_size();

  //Recalculate the minimum buffer size
  mBuff->sizeBuffer = bytesPerDecodeBuffer * resizeFactor;

  //FREE HOST AND DEVICE BUFFER
  GPU_ERROR(gpu_free_buffer(mBuff));

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  //ALLOCATE HOST AND DEVICE BUFFER
  CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
  CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));

  gpu_fmi_decode_init_buffer_(fmiBuffer);
}

/************************************************************
Debug functions for the index
************************************************************/

char gpu_fmi_search_bin_to_char(const uint32_t base)
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

uint32_t gpu_fmi_search_print_seed(const gpu_fmi_search_seed_t seed, const uint32_t seedSize)
{
  char plainSeed[GPU_FMI_SEED_MAX_CHARS] = {0};
  uint64_t bitmap = seed.hi;
  uint32_t idBase;

  for(idBase = 0; idBase < seedSize; ++idBase){
    uint32_t base = bitmap & 0x3;
    plainSeed[idBase] = gpu_fmi_search_bin_to_char(base);
    bitmap >>= GPU_FMI_SEED_CHAR_LENGTH;
    if(idBase == GPU_UINT32_LENGTH) bitmap = seed.low;
  }

  for(idBase = 0; idBase < seedSize; ++idBase)
    printf("%c", plainSeed[seedSize - idBase - 1]);

  return(SUCCESS);
}

uint32_t flag_print = 0;
uint32_t gpu_fmi_search_print_buffer(const void* const fmiBuffer)
{
  if(flag_print == 0){
    gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
    const uint32_t maxSeeds    = 100; // Just check and print the first results
    const uint32_t numSeeds    = mBuff->data.search.seeds.numSeeds;
          uint32_t missMatches = 0;
          uint32_t idSeed;

    printf("Buffer: %d ------------------------------------\n", mBuff->idBuffer);
    for(idSeed = 0; idSeed < numSeeds; ++idSeed){
      const uint64_t hiSeedSection = mBuff->data.search.seeds.h_seeds[idSeed].low;
      const uint32_t seedSize = hiSeedSection >> (GPU_UINT64_LENGTH - GPU_FMI_SEED_FIELD_SIZE);
      printf("[%d] seed=", idSeed);
      gpu_fmi_search_print_seed(mBuff->data.search.seeds.h_seeds[idSeed], seedSize);
      printf("\t size=%d \t (GPU) lo=%llu \t hi=%llu \n",
             seedSize,
             mBuff->data.search.saIntervals.h_intervals[idSeed].low,
             mBuff->data.search.saIntervals.h_intervals[idSeed].hi);
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

gpu_error_t gpu_fmi_search_transfer_CPU_to_GPU(gpu_buffer_t *mBuff)
{
  const gpu_fmi_search_seeds_buffer_t*  seedBuff = &mBuff->data.search.seeds;
  const cudaStream_t                    idStream =  mBuff->idStream;
  const size_t                          cpySize  =  seedBuff->numSeeds * sizeof(gpu_fmi_search_seed_t);

  //Transfer seeds from CPU to the GPU
  CUDA_ERROR(cudaMemcpyAsync(seedBuff->d_seeds, seedBuff->h_seeds, cpySize, cudaMemcpyHostToDevice, idStream));

  return (SUCCESS);
}

gpu_error_t gpu_fmi_search_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff)
{
  const gpu_fmi_search_sa_inter_buffer_t* interBuff = &mBuff->data.search.saIntervals;
  const cudaStream_t                      idStream  =  mBuff->idStream;
  const size_t                            cpySize   =  interBuff->numIntervals * sizeof(gpu_fmi_search_sa_inter_t);

  //Transfer SA intervals (occurrence results) from CPU to the GPU
  CUDA_ERROR(cudaMemcpyAsync(interBuff->h_intervals, interBuff->d_intervals, cpySize, cudaMemcpyDeviceToHost, idStream));

  return (SUCCESS);
}

void gpu_fmi_search_send_buffer_(void* const fmiBuffer, const uint32_t numSeeds)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;

  //Set real size of the input
  mBuff->data.search.seeds.numSeeds = numSeeds;
  mBuff->data.search.saIntervals.numIntervals = numSeeds;

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
  GPU_ERROR(gpu_fmi_search_transfer_CPU_to_GPU(mBuff));
  GPU_ERROR(gpu_fmi_search_process_buffer(mBuff));
  GPU_ERROR(gpu_fmi_search_transfer_GPU_to_CPU(mBuff));
}

void gpu_fmi_search_receive_buffer_(const void* const fmiBuffer)
{
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;

  //Synchronize Stream (the thread wait for the commands done in the stream)
  CUDA_ERROR(cudaStreamSynchronize(mBuff->idStream));
}

/************************************************************
Functions to transfer data HOST <-> DEVICE (Decode)
************************************************************/

gpu_error_t gpu_fmi_decode_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff)
{
  const gpu_fmi_decode_init_pos_buffer_t* initPosBuff = &mBuff->data.decode.initPositions;
  const cudaStream_t                      idStream    =  mBuff->idStream;
  const size_t                            cpySize     =  initPosBuff->numDecodings * sizeof(gpu_fmi_decode_init_pos_t);

  //Transfer seeds from CPU to the GPU
  CUDA_ERROR(cudaMemcpyAsync(initPosBuff->d_initBWTPos, initPosBuff->h_initBWTPos, cpySize, cudaMemcpyHostToDevice, idStream));

  return (SUCCESS);
}

gpu_error_t gpu_fmi_decode_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff)
{
  const gpu_fmi_decode_end_pos_buffer_t*  endPosBuff = &mBuff->data.decode.endPositions;
  const cudaStream_t                      idStream   =  mBuff->idStream;
  const size_t                            cpySize    =  endPosBuff->numDecodings * sizeof(gpu_fmi_decode_end_pos_t);

  //Transfer SA intervals (occurrence results) from CPU to the GPU
  CUDA_ERROR(cudaMemcpyAsync(endPosBuff->h_endBWTPos, endPosBuff->d_endBWTPos, cpySize, cudaMemcpyDeviceToHost, idStream));

  return (SUCCESS);
}

void gpu_fmi_decode_send_buffer_(void* const fmiBuffer, const uint32_t numDecodings, const uint32_t samplingRate)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;

  //Set real size of the input
  mBuff->data.decode.initPositions.numDecodings = numDecodings;
  mBuff->data.decode.endPositions.numDecodings  = numDecodings;
  mBuff->data.decode.samplingRate               = samplingRate;

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  GPU_ERROR(gpu_fmi_decode_transfer_CPU_to_GPU(mBuff));
  GPU_ERROR(gpu_fmi_decode_process_buffer(mBuff));
  GPU_ERROR(gpu_fmi_decode_transfer_GPU_to_CPU(mBuff));
}

void gpu_fmi_decode_receive_buffer_(const void* const fmiBuffer)
{
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;

  //Synchronize Stream (the thread wait for the commands done in the stream)
  CUDA_ERROR(cudaStreamSynchronize(mBuff->idStream));
}
