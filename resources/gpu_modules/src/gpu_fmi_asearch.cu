/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_SEARCH_C_
#define GPU_FMI_SEARCH_C_

#include "../include/gpu_fmi_core.h"

GPU_INLINE __device__ void get_base(const char currentBase, uint32_t* const bit0, uint32_t* const bit1, bool* const foundN)
{
  // Gathering the base of the seed
  (* bit0)   =  currentBase & 0x1L;
  (* bit1)   = (currentBase & 0x2L) >> 1;
  (* foundN) = (currentBase & 0x4L) >> 2;
}

GPU_INLINE __device__ void advance_step_LF_mapping(const gpu_fmi_device_entry_t* const fmi, const uint32_t bit0, const uint32_t bit1,
                                                   uint64_t* const L, uint64_t* const R, gpu_fmi_exch_bmp_mem_t * const seedExchBMP)
{
  const uint32_t globalThreadIdx     = gpu_get_thread_idx();
  const uint32_t localWarpThreadIdx  = globalThreadIdx     % GPU_WARP_SIZE;
  const uint32_t localEntryIdx       = localWarpThreadIdx  / GPU_FMI_THREADS_PER_ENTRY;
  const uint32_t localEntryThreadIdx = localWarpThreadIdx  % GPU_FMI_THREADS_PER_ENTRY;

  // Communicate along the FMI entry group threads the L o R interval
  uint64_t       interval       = (localEntryIdx % GPU_FMI_ENTRIES_PER_QUERY) ? (* R) : (* L);
  const uint64_t entryIdx       =  interval / GPU_FMI_ENTRY_SIZE;
  const uint32_t bitmapPosition =  interval % GPU_FMI_ENTRY_SIZE;

  // Loading FM-index entry in thread cooperative way
  const uint32_t missedEntry   = (entryIdx % GPU_FMI_ALTERNATE_COUNTERS != bit1) ? 1 : 0;
  const uint64_t entryIdxFixed = (localEntryThreadIdx == 0) ? entryIdx + missedEntry : entryIdx;
  uint4 loadEntry              = fmi[entryIdxFixed].v[localEntryThreadIdx];

  // Compute LF-Mapping (th0 of each group contain the result)
  interval = LF_Mapping(loadEntry, seedExchBMP, missedEntry, localEntryThreadIdx, bitmapPosition, bit1, bit0);

  // Update interval & communicate the th0 (L, R) to the rest of group threads
  const uint32_t lane = GPU_SELECT_OFFSET(localWarpThreadIdx, GPU_FMI_THREADS_PER_QUERY);
  (* L) = shfl_64(interval, lane);
  (* R) = shfl_64(interval, lane + GPU_FMI_THREADS_PER_ENTRY);

}

#define GPU_FMI_STEPS             4
#define GPU_FMI_OCC_THRESHOLD     20  //Min number of regions per seed

void __global__ gpu_fmi_asearch_kernel(const gpu_fmi_device_entry_t* const fmi, const uint64_t bwtSize,
                                       const char* const queries, const uint2* const queryInfo, const uint32_t numQueries,
                                       uint2* const regions, ulonglong2* const regIntervals, uint2* const regOffset,
                                       const uint32_t maxExtraSteps, const uint32_t occThreshold, const uint32_t occShrinkFactor,
                                       const uint32_t maxRegionsFactor)
{
  // Thread group-scheduling initializations
  const uint32_t globalThreadIdx     = gpu_get_thread_idx();
  const uint32_t localWarpThreadIdx  = globalThreadIdx     % GPU_WARP_SIZE;
  const uint32_t idQuery             = globalThreadIdx     / GPU_FMI_THREADS_PER_QUERY;

  if (idQuery < numQueries){
    // Setting thread buffer affinity
    const uint32_t    querySize      = queryInfo[idQuery].y;
    const char* const query          = queries + queryInfo[idQuery].x + querySize - 1; //end of query
    ulonglong2* const regionInterval = regIntervals + regions[idQuery].x;
    uint2* const      regionOffset   = regOffset    + regions[idQuery].x;

    //Shared memory space dedicated for the internal FMI entry thread-communications
    __shared__ gpu_fmi_exch_bmp_mem_t   exchBMP[GPU_FMI_ENTRIES_PER_BLOCK];
               gpu_fmi_exch_bmp_mem_t * const seedExchBMP = &exchBMP[threadIdx.x / GPU_FMI_THREADS_PER_ENTRY];

    const uint32_t maxRegions = GPU_MAX(GPU_DIV_CEIL(querySize, maxRegionsFactor), GPU_FMI_MIN_REGIONS);
          uint32_t idBase = 0, idRegion = 0;
    //Extracts and locates each seed
    while ((idBase < querySize) && (idRegion < maxRegions)){
      //Region initializations
      uint64_t L = 0, R = bwtSize;
      uint64_t occ = R - L, initBase = idBase, endBase = idBase;
      uint32_t bit0, bit1;
      bool     foundN = false;
      //Searching for the next seed
      while((occ > occThreshold) && (idBase && querySize) && !foundN){
        // Gathering the base of the seed
        get_base(LDG(query - idBase), &bit0, &bit1, &foundN);
        if(!foundN) advance_step_LF_mapping(fmi, bit0, bit1, &L, &R, seedExchBMP);
        idBase++;
        occ = R - L;
      }
      //Evaluate current seed (discard or continue exploration)
      if(occ < occThreshold){
        uint64_t endL = L, endR = R;
        endBase = idBase;
        if(!foundN){
          // Extension initialization
          uint32_t shrinkOccThreshold = occ >> occShrinkFactor, idStep = 0;
          //Last steps extension (exploration for consecutive 4 bases)
          while((idBase < querySize) && (occ != 0) && (idStep < maxExtraSteps) && !foundN){
            // Gathering the base of the seed
            get_base(LDG(query - idBase), &bit0, &bit1, &foundN);
            if(!foundN) advance_step_LF_mapping(fmi, bit0, bit1, &L, &R, seedExchBMP);
            // Update seed information
            occ = R - L;
            if((occ <= shrinkOccThreshold) && (occ != 0)){
              endL = L; endR = R;
              endBase = idBase;
              shrinkOccThreshold = occ;
            }
            // Move to next base
            shrinkOccThreshold >>= occShrinkFactor;
            idBase++; idStep++;
          }
        }
        // Save extracted region (SA intervals + Query position)
        if((localWarpThreadIdx % GPU_FMI_THREADS_PER_QUERY) == 0){
            regionInterval[idRegion] = make_ulonglong2(endL, endR);
            regionOffset[idRegion]   = make_uint2(initBase, endBase);
        }
        idRegion++;
      }
    }
    // Save region profile info (number of extracted regions)
    if((localWarpThreadIdx % GPU_FMI_THREADS_PER_QUERY) == 0)
      regions[idQuery].y = idRegion;
  }
}

extern "C"
gpu_error_t gpu_fmi_asearch_launch_kernel(const gpu_fmi_device_entry_t* const d_fmi, const uint64_t bwtSize,
                                          const char* const d_queries, const uint2* const d_queryInfo, const uint32_t numQueries,
                                          uint2* const d_regions, ulonglong2* const d_regIntervals, uint2* const d_regOffset)
{
  //GPU thread configuration
  const uint32_t threads = 128;
  const uint32_t blocks  = GPU_DIV_CEIL(numQueries * GPU_FMI_THREADS_PER_QUERY, threads);
  const uint32_t nreps   = 10;
  //Search configuration
  const uint32_t occShrinkFactor  =  GPU_FMI_D;
  const uint32_t occThreshold     =  GPU_FMI_OCC_THRESHOLD;
  const uint32_t maxExtraSteps    =  GPU_FMI_STEPS;
  const uint32_t maxRegionsFactor =  GPU_FMI_RATIO_REGIONS;

  float elapsed_time_ms = 0.0f;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

    for(uint32_t iteration = 0; iteration < nreps; ++iteration)
      gpu_fmi_asearch_kernel<<<blocks,threads>>>(d_fmi, bwtSize, d_queries, d_queryInfo, numQueries,
                                                 d_regions, d_regIntervals, d_regOffset,
                                                 maxExtraSteps, occThreshold, occShrinkFactor, maxRegionsFactor);

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
gpu_error_t gpu_fmi_asearch_process_buffer(gpu_buffer_t* const mBuff)
{
  // Getting internal buffers
  const gpu_index_buffer_t* const           index            =  mBuff->index;
  const gpu_fmi_asearch_queries_buffer_t*   queries          = &mBuff->data.asearch.queries;
  const gpu_fmi_asearch_regions_buffer_t*   regions          = &mBuff->data.asearch.regions;
  // Getting buffer sizes
  const uint32_t                            numQueries       =  mBuff->data.asearch.queries.numQueries;
  const uint32_t                            numMaxQueries    =  mBuff->data.asearch.numMaxQueries;
  const uint32_t                            numBases         =  mBuff->data.asearch.queries.numBases;
  const uint32_t                            numMaxBases      =  mBuff->data.asearch.numMaxBases;
  // Getting device information
  const cudaStream_t                        idStream         =  mBuff->idStream;
  const uint32_t                            idSupDev         =  mBuff->idSupportedDevice;
  const gpu_device_info_t* const            device           =  mBuff->device[idSupDev];
  //Search configuration
  const uint32_t                            occShrinkFactor  =  BASE2LOG(mBuff->data.asearch.alphabetSize);
  const uint32_t                            occThreshold     =  mBuff->data.asearch.occMinThreshold;
  const uint32_t                            maxExtraSteps    =  mBuff->data.asearch.extraSteps;
  const uint32_t                            maxRegionsFactor =  mBuff->data.asearch.maxRegionsFactor;

  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = numQueries * GPU_FMI_THREADS_PER_QUERY;
  gpu_device_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);
  // Sanity-check (checks buffer overflowing)
  if((numQueries > numMaxQueries) || (numBases > numMaxBases))
    return(E_OVERFLOWING_BUFFER);

  gpu_fmi_asearch_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((gpu_fmi_device_entry_t*) index->fmi.d_fmi[idSupDev], index->fmi.bwtSize,
                                                                          (char*) queries->d_queries, (uint2*) queries->d_queryInfo, numQueries,
                                                                          (uint2*) queries->d_regions, (ulonglong2*) regions->d_intervals, (uint2*) regions->d_regionsOffsets,
                                                                          maxExtraSteps, occThreshold, occShrinkFactor, maxRegionsFactor);
  return(SUCCESS);
}

#endif /* GPU_FMI_SEARCH_CU_ */
