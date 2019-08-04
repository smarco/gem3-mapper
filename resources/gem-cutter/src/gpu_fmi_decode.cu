/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_DECODE_CU_
#define GPU_FMI_DECODE_CU_

#include "../include/gpu_fmi_core.h"

void __global__ gpu_fmi_decoding_kernel(const gpu_fmi_device_entry_t* const fmi, const uint64_t bwtSize, const uint32_t numDecodings,
                                        const uint64_t* const d_initBWTPos, ulonglong2* const d_endBWTPos, const uint32_t samplingRate)
{
  const uint32_t globalThreadIdx      = gpu_get_thread_idx();
  const uint32_t localWarpThreadIdx   = globalThreadIdx % GPU_WARP_SIZE;
  const uint32_t idDecoding           = globalThreadIdx / GPU_FMI_DECODE_THREADS_PER_ENTRY;

  if(idDecoding < numDecodings){

    const uint32_t decodeEntryIdx        = localWarpThreadIdx  / GPU_FMI_DECODE_THREADS_PER_ENTRY;
    const uint32_t decodeEntryThreadIdx  = localWarpThreadIdx  % GPU_FMI_DECODE_THREADS_PER_ENTRY;
    const uint32_t fmiEntryThreadIdx     = localWarpThreadIdx  % GPU_FMI_THREADS_PER_ENTRY;
    const uint32_t loadFMIEntryThreadIdx = generate_load_threadIdx(decodeEntryThreadIdx);

          uint4    loadEntry;
          uint64_t interval = d_initBWTPos[idDecoding], idStep = 0;
          uint32_t foundBaseN = 0;

    __shared__ gpu_fmi_exch_bmp_mem_t   exchBMP[GPU_FMI_ENTRIES_PER_BLOCK];
               gpu_fmi_exch_bmp_mem_t * const decExchBMP = &exchBMP[threadIdx.x / GPU_FMI_THREADS_PER_ENTRY];

    while((interval % samplingRate) && (foundBaseN == 0)){
      const uint64_t entryIdx       =  interval / GPU_FMI_ENTRY_SIZE;
      const uint32_t bitmapPosition =  interval % GPU_FMI_ENTRY_SIZE;
            uint32_t bit1, bit0, bit2;

      // Loading FM-index entry in thread cooperative way
      if(( 0 < decodeEntryThreadIdx) && (decodeEntryThreadIdx < GPU_FMI_DECODE_THREADS_PER_LOAD))
        loadEntry = fmi[entryIdx].v[loadFMIEntryThreadIdx];

      // Gathering the base and sharing it with the rest of the threads
      gather_base_from_BWT(loadEntry, decExchBMP, bitmapPosition, fmiEntryThreadIdx, decodeEntryIdx, &bit1, &bit0, &bit2);

      // Gathering the counters
      const uint32_t missedEntry = (entryIdx % GPU_FMI_ALTERNATE_COUNTERS != bit1) ? 1 : 0;
      loadEntry = gather_counters_FMI(loadEntry, missedEntry, decodeEntryIdx, decodeEntryThreadIdx);

      // Compute LF-Mapping th0 of each group contain the result)
      interval = LF_Mapping(loadEntry, decExchBMP, missedEntry, fmiEntryThreadIdx, bitmapPosition, bit1, bit0);

      // Share interval along the thread group
      const uint32_t lane = decodeEntryIdx * GPU_FMI_DECODE_THREADS_PER_ENTRY;
      interval        = shfl_64(interval, lane);

      // Exit condition
      if(bit2 == 0) foundBaseN = 1;

      // Increment for the next Backward-Search
      idStep++;
    }

    ulonglong2 BWTPos = {interval, idStep};
    if(foundBaseN) BWTPos = make_ulonglong2(GPU_UINT64_MAX_VALUE, GPU_UINT64_MAX_VALUE);

    // Save intervals
    if(decodeEntryThreadIdx  == 0) d_endBWTPos[idDecoding] = BWTPos;
  }
}

extern "C"
gpu_error_t gpu_fmi_decoding_launch_kernel(const gpu_fmi_device_entry_t* const d_fmi, const uint64_t bwtSize,
                                           const uint32_t numDecodings, const uint64_t* const d_initBWTPos,
                                           ulonglong2* const d_endBWTPos)
{
  const uint32_t samplingRate = 4;
  const uint32_t threads = 128;
  const uint32_t blocks  = GPU_DIV_CEIL(numDecodings * GPU_FMI_DECODE_THREADS_PER_ENTRY, threads);
  const uint32_t nreps   = 10;

  printf("GPU numDecodings=%u \n", numDecodings);

  float elapsed_time_ms = 0.0f;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

    for(uint32_t iteration = 0; iteration < nreps; ++iteration)
      gpu_fmi_decoding_kernel<<<blocks,threads>>>(d_fmi, bwtSize, numDecodings, d_initBWTPos, d_endBWTPos, samplingRate);

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
gpu_error_t gpu_fmi_decode_process_buffer(gpu_buffer_t* const mBuff)
{
  const gpu_index_buffer_t* const               index               =  mBuff->index;
  const gpu_fmi_decode_buffer_t* const          decBuff             = &mBuff->data.decode;
  const gpu_fmi_decode_init_pos_buffer_t* const initPos             = &mBuff->data.decode.initPositions;
  const gpu_fmi_decode_end_pos_buffer_t* const  endPos              = &mBuff->data.decode.endPositions;
  const uint32_t                                numDecodings        =  mBuff->data.decode.initPositions.numDecodings;
  const uint32_t                                numMaxInitPositions =  mBuff->data.decode.numMaxInitPositions;
  const uint32_t                                numMaxEndPositions  =  mBuff->data.decode.numMaxEndPositions;
  const cudaStream_t                            idStream            =  mBuff->listStreams[mBuff->idStream];
  const uint32_t                                idSupDev            =  mBuff->idSupportedDevice;
  const gpu_device_info_t* const                device              =  mBuff->device[idSupDev];

  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = numDecodings * GPU_FMI_DECODE_THREADS_PER_ENTRY;
  gpu_device_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);
  // Sanity-check (checks buffer overflowing)
  if((numDecodings > numMaxInitPositions) || (numDecodings > numMaxEndPositions)){
    return(E_OVERFLOWING_BUFFER);
  }

  gpu_fmi_decoding_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((gpu_fmi_device_entry_t*) index->fmi.d_fmi[idSupDev], index->fmi.bwtSize,
                                                                           initPos->numDecodings, (uint64_t*) initPos->d_initBWTPos,
                                                                           (ulonglong2*) endPos->d_endBWTPos, decBuff->samplingRate);
  return(SUCCESS);
}

#endif /* GPU_FMI_DECODE_CU_ */


