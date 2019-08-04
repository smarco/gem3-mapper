/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_BPM_CORE_H_
#define GPU_BPM_CORE_H_

extern "C" {
#include "gpu_commons.h"
#include "gpu_reference.h"
#include "gpu_buffer.h"
}
#include "gpu_resources.h"

__device__ const gpu_char4_t GPU_CIGAR_INIT = {(char)0, (char)0, GPU_CIGAR_NULL, (char)0};
                                                                                 //(x,y) Matrix DP movement; CIGAR op; Match efficiency;
__device__ const gpu_bpm_align_device_cigar_entry_t GPU_CIGAR_DEVICE_NULL      = {(char) 0, (char) 0, GPU_CIGAR_NULL,      (char) 0};
__device__ const gpu_bpm_align_device_cigar_entry_t GPU_CIGAR_DEVICE_MATCH     = {(char)-1, (char)-1, GPU_CIGAR_MATCH,     (char) 0};
__device__ const gpu_bpm_align_device_cigar_entry_t GPU_CIGAR_DEVICE_MISSMATCH = {(char)-1, (char)-1, GPU_CIGAR_MISSMATCH, (char) 0};
__device__ const gpu_bpm_align_device_cigar_entry_t GPU_CIGAR_DEVICE_INSERTION = {(char)-1, (char) 0, GPU_CIGAR_INSERTION, (char) 1};
__device__ const gpu_bpm_align_device_cigar_entry_t GPU_CIGAR_DEVICE_DELETION  = {(char) 0, (char)-1, GPU_CIGAR_DELETION,  (char)-1};

                                                                                   //List of cigar events
__device__ const gpu_bpm_align_device_cigar_entry_t gpu_bmp_align_cigar_events[] = {{(char) 0, (char) 0, GPU_CIGAR_NULL,      (char) 0},
                                                                                    {(char)-1, (char)-1, GPU_CIGAR_MATCH,     (char) 0},
                                                                                    {(char)-1, (char)-1, GPU_CIGAR_MISSMATCH, (char) 0},
                                                                                    {(char)-1, (char) 0, GPU_CIGAR_INSERTION, (char) 1},
                                                                                    {(char) 0, (char)-1, GPU_CIGAR_DELETION,  (char)-1}};

__device__ const gpu_bpm_align_device_cigar_entry_t gpu_bmp_align_cigar_lut[16] = {//Left Cigar Alignment
                                                                                   {(char)-1, (char)-1, GPU_CIGAR_MISSMATCH, (char) 0},
                                                                                   {(char)-1, (char)-1, GPU_CIGAR_MATCH,     (char) 0},
                                                                                   {(char)-1, (char) 0, GPU_CIGAR_INSERTION, (char) 1},
                                                                                   {(char)-1, (char)-1, GPU_CIGAR_MATCH,     (char) 0},
                                                                                   {(char) 0, (char)-1, GPU_CIGAR_DELETION,  (char)-1},
                                                                                   {(char)-1, (char)-1, GPU_CIGAR_MATCH,     (char) 0},
                                                                                   {(char) 0, (char)-1, GPU_CIGAR_DELETION,  (char)-1},
                                                                                   {(char)-1, (char)-1, GPU_CIGAR_MATCH,     (char) 0},
                                                                                   //Right Cigar Alignment
                                                                                   {(char)-1, (char)-1, GPU_CIGAR_MISSMATCH, (char) 0},
                                                                                   {(char)-1, (char)-1, GPU_CIGAR_MATCH,     (char) 0},
                                                                                   {(char)-1, (char) 0, GPU_CIGAR_INSERTION, (char) 1},
                                                                                   {(char)-1, (char) 0, GPU_CIGAR_INSERTION, (char) 1},
                                                                                   {(char) 0, (char)-1, GPU_CIGAR_DELETION,  (char)-1},
                                                                                   {(char) 0, (char)-1, GPU_CIGAR_DELETION,  (char)-1},
                                                                                   {(char) 0, (char)-1, GPU_CIGAR_DELETION,  (char)-1},
                                                                                   {(char) 0, (char)-1, GPU_CIGAR_DELETION,  (char)-1}};

GPU_INLINE __device__ void cooperative_shift(uint32_t *value, const uint32_t shiftedBits,
                                             const uint32_t localThreadIdx, const uint32_t BMPS_PER_THREAD)
{
  const uint32_t laneIdx = threadIdx.x % GPU_WARP_SIZE;
  uint32_t carry;

  carry = shfl_32(value[BMPS_PER_THREAD-1], laneIdx - 1);
  carry = (localThreadIdx) ? carry : 0;
  #pragma unroll
  for(uint32_t idBMP = BMPS_PER_THREAD-1; idBMP > 0; --idBMP){
    value[idBMP] = funnelshift_lc_32(value[idBMP-1], value[idBMP], shiftedBits);
  }
  value[0] = funnelshift_lc_32(carry, value[0], shiftedBits);
}

GPU_INLINE __device__ void cooperative_sum(const uint32_t *A, const uint32_t *B, uint32_t *C,
                                           const uint32_t localThreadIdx, const uint32_t BMPS_PER_THREAD)

{
  const uint32_t laneIdx = threadIdx.x % GPU_WARP_SIZE;
  uint32_t carry;

  UADD__CARRY_OUT(C[0], A[0], B[0]);
  #pragma unroll
  for(uint32_t idBMP = 1; idBMP < BMPS_PER_THREAD; ++idBMP){
    UADD__IN_CARRY_OUT(C[idBMP], A[idBMP], B[idBMP]);
  }
  UADD__IN_CARRY(carry, 0, 0);

  while(any_32(carry)){
    carry = shfl_32(carry, laneIdx - 1);
    carry = (localThreadIdx) ? carry : 0;
    UADD__CARRY_OUT(C[0], C[0], carry);
    #pragma unroll
    for(uint32_t idBMP = 1; idBMP < BMPS_PER_THREAD; ++idBMP){
      UADD__IN_CARRY_OUT(C[idBMP], C[idBMP], 0);
    }
    UADD__IN_CARRY(carry, 0, 0);
  }
}

#endif /* GPU_BPM_CORE_H_ */

