/*
 * PROJECT: Thread-cooperative FM-index on GPU
 * FILE: genIndex.h
 * DATE: 1/9/2015
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: FM-Index DNA backward-search customized for GEM Mapper
 */

#ifndef GPU_SA_DECODE_C_
#define GPU_SA_DECODE_C_

#include "../include/gpu_sa_core.h"

void __global__ gpu_sa_decoding_kernel(const uint64_t* const d_SA, const uint32_t numDecodings,
                                       const ulonglong2* const d_endBWTPos, uint64_t* const d_textPos)
{
  const uint32_t idDecoding = gpu_get_thread_idx();

  if(idDecoding < numDecodings){
    const ulonglong2 saPosition = d_endBWTPos[idDecoding];
    if((saPosition.x < GPU_UINT64_MAX_VALUE) && (saPosition.x < GPU_UINT64_MAX_VALUE))
      d_textPos[idDecoding] = d_SA[saPosition.x] + saPosition.y;
  }
}

extern "C"
gpu_error_t gpu_sa_decoding_launch_kernel(const uint64_t* const d_SA, const uint32_t numDecodings,
                                          const ulonglong2* const d_endBWTPos, uint64_t* const d_textPos)
{
  const uint32_t threads = 128;
  const uint32_t blocks  = GPU_DIV_CEIL(numDecodings, threads);
  const uint32_t nreps   = 10;

  float elapsed_time_ms = 0.0f;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

    for(uint32_t iteration = 0; iteration < nreps; ++iteration)
      gpu_sa_decoding_kernel<<<blocks,threads>>>(d_SA, numDecodings, d_endBWTPos, d_textPos);

  cudaEventRecord(stop, 0);
  cudaThreadSynchronize();
  cudaEventElapsedTime(&elapsed_time_ms, start, stop);
  elapsed_time_ms /= nreps;

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  printf("\t Time Kernel GPU:  %8.2f ms\n", elapsed_time_ms);

  return(SUCCESS);
}

extern "C"
gpu_error_t gpu_sa_decode_process_buffer(gpu_buffer_t* const mBuff)
{
  gpu_index_buffer_t                  *index        =  mBuff->index;
  gpu_fmi_decode_end_pos_buffer_t     *endPos       = &mBuff->data.decode.endPositions;
  gpu_fmi_decode_text_pos_buffer_t    *textPos      = &mBuff->data.decode.textPositions;
  uint32_t                            numDecodings  =  mBuff->data.decode.endPositions.numDecodings;
  cudaStream_t                        idStream      =  mBuff->idStream;
  uint32_t                            idSupDev      =  mBuff->idSupportedDevice;
  gpu_device_info_t                   *device       =  mBuff->device[idSupDev];

  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = numDecodings;
  gpu_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);

  gpu_sa_decoding_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>(index->sa.d_sa[idSupDev], numDecodings, (ulonglong2*) endPos->d_endBWTPos,
                                                                          (uint64_t*) textPos->d_textPos);

  return(SUCCESS);
}

#endif /* GPU_SA_DECODE_C_ */



