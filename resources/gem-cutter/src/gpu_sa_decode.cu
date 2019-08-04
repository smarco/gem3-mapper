/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_SA_DECODE_C_
#define GPU_SA_DECODE_C_

#include "../include/gpu_sa_core.h"

void __global__ gpu_sa_decoding_kernel(const uint64_t* const d_SA, const uint32_t samplingRate,
                                       const uint32_t numDecodings, const ulonglong2* const d_endBWTPos,
                                       uint64_t* const d_textPos)
{
  const uint32_t idDecoding = gpu_get_thread_idx();

  if(idDecoding < numDecodings){
    const ulonglong2 saPosition   = d_endBWTPos[idDecoding];
          uint64_t   textPosition = GPU_UINT64_MAX_VALUE;

    if((saPosition.x < GPU_UINT64_MAX_VALUE) && (saPosition.y < GPU_UINT64_MAX_VALUE))
      textPosition = d_SA[saPosition.x / samplingRate] + saPosition.y;
    
    d_textPos[idDecoding] = textPosition;
  }
}

extern "C"
gpu_error_t gpu_sa_decoding_launch_kernel(const uint64_t* const d_SA, const uint32_t samplingRate,
		                                      const uint32_t numDecodings, const ulonglong2* const d_endBWTPos,
		                                      uint64_t* const d_textPos)
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
      gpu_sa_decoding_kernel<<<blocks,threads>>>(d_SA, samplingRate, numDecodings, d_endBWTPos, d_textPos);

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
  const gpu_index_buffer_t* const               index               =  mBuff->index;
  const gpu_fmi_decode_end_pos_buffer_t* const  endPos              = &mBuff->data.decode.endPositions;
  const gpu_fmi_decode_text_pos_buffer_t* const textPos             = &mBuff->data.decode.textPositions;
  const uint32_t                                numDecodings        =  mBuff->data.decode.textPositions.numDecodings;
  const uint32_t                                numMaxEndPositions  =  mBuff->data.decode.numMaxEndPositions;
  const uint32_t                                numMaxTextPositions =  mBuff->data.decode.numMaxTextPositions;
  const uint32_t                                samplingRate        =  mBuff->data.decode.samplingRate;
  const cudaStream_t                            idStream            =  mBuff->listStreams[mBuff->idStream];
  const uint32_t                                idSupDev            =  mBuff->idSupportedDevice;
  const gpu_device_info_t* const                device              =  mBuff->device[idSupDev];

  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = numDecodings;
  gpu_device_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);
  // Sanity-check (checks buffer overflowing)
  if((numDecodings > numMaxEndPositions) || (numDecodings > numMaxTextPositions)){
    printf("DECODING SA textPositions.numDecodings=%u, endPositions.numDecodings=%u, maxDecodings=%u \n", numDecodings, mBuff->data.decode.endPositions.numDecodings, mBuff->data.decode.numMaxTextPositions);
    return(E_OVERFLOWING_BUFFER);
  }

  gpu_sa_decoding_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>(index->sa.d_sa[idSupDev], samplingRate, numDecodings,
		                                                                      (ulonglong2*) endPos->d_endBWTPos, (uint64_t*) textPos->d_textPos);

  return(SUCCESS);
}

#endif /* GPU_SA_DECODE_C_ */



