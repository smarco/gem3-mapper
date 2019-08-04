/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

extern "C" {
#include "gpu_commons.h"
#include "gpu_index.h"
#include "gpu_buffer.h"
}
#include "gpu_resources.h"


#ifndef GPU_FMI_CORE_H_
#define GPU_FMI_CORE_H_

GPU_INLINE __device__ uint32_t generate_load_threadIdx(const uint64_t decodeEntryThreadIdx)
{
  uint32_t localThreadIdx = (decodeEntryThreadIdx  > (GPU_FMI_THREADS_PER_ENTRY + 1)) ? 0                         : decodeEntryThreadIdx;
           localThreadIdx = (decodeEntryThreadIdx == (GPU_FMI_THREADS_PER_ENTRY + 1)) ? GPU_FMI_THREADS_PER_ENTRY : localThreadIdx;
           localThreadIdx = (decodeEntryThreadIdx ==  GPU_FMI_THREADS_PER_ENTRY)      ? 0                         : localThreadIdx;
  return(localThreadIdx);
}

GPU_INLINE __device__ uint32_t count_bitmap(const uint32_t bitmap, const int32_t shift, const uint32_t idxCounterGroup)
{
  uint32_t mask = GPU_UINT32_ONES << (GPU_UINT32_LENGTH - shift);
           mask = (shift > GPU_UINT32_LENGTH) ? GPU_UINT32_ONES : mask;
           mask = (shift > 0) ? mask : GPU_UINT32_ZEROS;
           mask = (idxCounterGroup) ? ~mask : mask;
  return (__popc(bitmap & mask));
}

GPU_INLINE __device__ void gather_bitmaps(const uint4 loadData, volatile gpu_fmi_exch_bmp_mem_t* const exchBMP, const uint32_t localEntryThreadIdx)
{
  const uint32_t idBMP = (localEntryThreadIdx == 0) ?  GPU_FMI_THREADS_PER_ENTRY - 1 : localEntryThreadIdx - 1;
  exchBMP->v[idBMP] = loadData.w;
}

GPU_INLINE __device__ uint32_t compute_bitmaps(uint4 bitmap, const uint32_t bitmapPosition, const uint32_t bit0, const uint32_t bit1,
                                               const uint32_t localEntryThreadIdx, const uint32_t idxCounterGroup)
{
  const int32_t  relativePosition = bitmapPosition - (localEntryThreadIdx * GPU_UINT32_LENGTH);
  uint32_t resultBitmaps, bmpCollapsed, numCaracters;

  bitmap.x = bit0 ? bitmap.x : ~bitmap.x;
  bitmap.y = bit1 ? bitmap.y : ~bitmap.y;

  bmpCollapsed  = bitmap.x & bitmap.y & bitmap.z;
  resultBitmaps = count_bitmap(bmpCollapsed, relativePosition, idxCounterGroup);
  numCaracters  = reduce_add_warp(resultBitmaps, GPU_FMI_THREADS_PER_ENTRY);

  return (numCaracters);
}

GPU_INLINE __device__ uint64_t select_counter(const uint4 vectorCounters, const uint32_t indexBase)
{
  const uint2 vectorCountersA = {vectorCounters.x, vectorCounters.y};
  const uint2 vectorCountersB = {vectorCounters.z, vectorCounters.w};
  uint64_t scalarCountersA, scalarCountersB;

  V2S_B64(vectorCountersA, scalarCountersA);
  V2S_B64(vectorCountersB, scalarCountersB);

  return((indexBase == 0) ? scalarCountersA : scalarCountersB);
}

GPU_INLINE __device__ void gather_base_from_BWT(uint4 vbitmap, gpu_fmi_exch_bmp_mem_t* const exchBMP, const uint32_t bitmapPosition,
                                                const uint32_t fmiEntryThreadIdx, const uint32_t decodeEntryIdx,
                                                uint32_t* bit1, uint32_t* bit0, uint32_t* bit2)
{
  uint32_t localBit0, localBit1, localBit2;

  gather_bitmaps(vbitmap, exchBMP, fmiEntryThreadIdx);
  if(fmiEntryThreadIdx == 0) vbitmap = exchBMP->s;

  const int32_t  relativePosition = (GPU_UINT32_LENGTH - (bitmapPosition % GPU_UINT32_LENGTH)) - 1;
  localBit0 = (vbitmap.x >> relativePosition) & 1;
  localBit1 = (vbitmap.y >> relativePosition) & 1;
  localBit2 = (vbitmap.z >> relativePosition) & 1;

  const uint32_t lane = (decodeEntryIdx * GPU_FMI_DECODE_THREADS_PER_ENTRY) + (bitmapPosition / GPU_UINT32_LENGTH);
  (* bit0) = shfl_32(localBit0, lane);
  (* bit1) = shfl_32(localBit1, lane);
  (* bit2) = shfl_32(localBit2, lane);
}

GPU_INLINE __device__ uint4 gather_counters_FMI(uint4 loadEntry, const uint32_t missedEntry, const uint32_t decodeEntryIdx, const uint32_t decodeEntryThreadIdx)
{
  uint4 auxData;
  const uint32_t idThreadCounter = (decodeEntryIdx * GPU_FMI_DECODE_THREADS_PER_ENTRY) + GPU_FMI_THREADS_PER_ENTRY;
  const uint32_t lane = (missedEntry) ? idThreadCounter + 1 : idThreadCounter;

  auxData.x = shfl_32(loadEntry.x, lane);
  auxData.y = shfl_32(loadEntry.y, lane);
  auxData.z = shfl_32(loadEntry.z, lane);
  auxData.w = shfl_32(loadEntry.w, lane);

  if(decodeEntryThreadIdx == 0) loadEntry = auxData;

  return(loadEntry);
}

GPU_INLINE __device__ uint64_t LF_Mapping(uint4 loadEntry, gpu_fmi_exch_bmp_mem_t* const exchBMP, const uint32_t missedEntry,
                                          const uint32_t localEntryThreadIdx, const uint32_t bitmapPosition,
                                          const uint32_t bit1, const uint32_t bit0)
{
  // Select the counter candidate
  const uint64_t resultCounters = select_counter(loadEntry, bit0);

  // Reorganize entry layout between threads (send bitmaps to thread 0 using shared mem)
  gather_bitmaps(loadEntry, exchBMP, localEntryThreadIdx);
  if(localEntryThreadIdx == 0) loadEntry = exchBMP->s;

  // Count the number of occ in the bitmap
  const uint32_t resultBitmaps  = compute_bitmaps(loadEntry, bitmapPosition, bit0, bit1, localEntryThreadIdx, missedEntry);

  // Compute interval alternate counters
  const uint64_t interval = (missedEntry) ? resultCounters - resultBitmaps : resultCounters + resultBitmaps;

  return(interval);
}

#endif /* GPU_FMI_CORE_H_ */

