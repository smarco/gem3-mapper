/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_RESOURCES_H_
#define GPU_RESOURCES_H_

extern "C" {
#include "gpu_commons.h"
#include "gpu_devices.h"
}

#define GPU_CC_FERMI_1G   200
#define GPU_CC_FERMI_2G   210
#define GPU_CC_KEPLER_1G  300
#define GPU_CC_KEPLER_2G  350
#define GPU_CC_KEPLER_3G  370
#define GPU_CC_MAXWELL_1G 500
#define GPU_CC_MAXWELL_2G 520

/* Defines related to GPU Architecture */
#if   (__CUDA_ARCH__ <= GPU_CC_FERMI_2G)
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_FERMI
#elif (__CUDA_ARCH__ <= GPU_CC_KEPLER_3G)
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_KEPLER
#elif (__CUDA_ARCH__ <= GPU_CC_MAXWELL_2G)
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_MAXWELL
#else
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_NEWGEN
#endif

/*************************************
GPU Side defines (ASM instructions)
**************************************/

// output temporal carry in internal register
#define UADD__CARRY_OUT(c, a, b) \
  asm volatile("add.cc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b))

// add & output with temporal carry of internal register
#define UADD__IN_CARRY_OUT(c, a, b) \
  asm volatile("addc.cc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b))

// add with temporal carry of internal register
#define UADD__IN_CARRY(c, a, b) \
  asm volatile("addc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b))

// packing and unpacking: from uint64_t to uint2
#define V2S_B64(v,s) \
  asm("mov.b64 %0, {%1,%2};" : "=l"(s) : "r"(v.x), "r"(v.y))

// packing and unpacking: from uint2 to uint64_t
#define S2V_B64(s,v) \
  asm("mov.b64 {%0,%1}, %2;" : "=r"(v.x), "=r"(v.y) : "l"(s))


/*************************************
DEVICE side basic block primitives
**************************************/

#if (__CUDA_ARCH__ < GPU_CC_KEPLER_2G)
  #define LDG(ptr)  (*(ptr))
#else
  #define LDG(ptr)  __ldg(ptr)
#endif

#if (__CUDA_ARCH__ < GPU_CC_KEPLER_1G)
__shared__ uint32_t interBuff[GPU_THREADS_PER_BLOCK];
GPU_INLINE __device__ uint32_t __emulated_shfl(const int scalarValue, const uint32_t source_lane)
{
  const uint32_t warpIdx = threadIdx.x / GPU_WARP_SIZE;
  const uint32_t laneIdx = threadIdx.x % GPU_WARP_SIZE;
  volatile uint32_t *interShuffle = interBuff + (warpIdx * GPU_WARP_SIZE);
  interShuffle[laneIdx] = scalarValue;
  return(interShuffle[source_lane % GPU_WARP_SIZE]);
}
#endif

GPU_INLINE __device__ uint32_t shfl_32(uint32_t scalarValue, const int lane)
{
  #if (__CUDA_ARCH__ < GPU_CC_KEPLER_1G)
    return (__emulated_shfl(scalarValue, (uint32_t)lane));
  #else
    return (__shfl((int)scalarValue, (int)lane, GPU_WARP_SIZE));
  #endif
}

GPU_INLINE __device__ uint64_t shfl_64(uint64_t scalarValue, const int lane)
{
  uint2 vectorValue;
  S2V_B64(scalarValue, vectorValue);
  vectorValue.x = shfl_32(vectorValue.x, lane);
  vectorValue.y = shfl_32(vectorValue.y, lane);
  V2S_B64(vectorValue, scalarValue);
  return (scalarValue);
}

GPU_INLINE __device__ uint32_t reduce_add_warp(uint32_t resultBitmaps, const uint32_t numThreads)
{
  const uint32_t lane_id = threadIdx.x % GPU_WARP_SIZE;
  for (int32_t i = 1; i < numThreads; i *= 2){
    int32_t n = shfl_32(resultBitmaps, lane_id + i);
    resultBitmaps += n;
  }
  return(resultBitmaps);
}

GPU_INLINE __device__ uint64_t funnelshift_left_64(const uint64_t scalarValueA,
                                                   const uint64_t scalarValueB,
                                                   const uint32_t shiftedBits)
{   /* Concat A with B and left-shifts nbits:  (A | B) << nbits */
  const uint32_t complementShiftedBits = GPU_UINT64_LENGTH - shiftedBits;
  return ( (scalarValueA << complementShiftedBits)
         | (scalarValueB >> shiftedBits));
}

GPU_INLINE __device__ uint32_t __emulated_funnelshift_right_32(const uint32_t scalarValueA,
                                                               const uint32_t scalarValueB,
                                                               const uint32_t shiftedBits)
{ /* (scalarValueB | scalarValueA) >> 1 : (stores higher 32 bits) */
  const uint32_t complementShiftedBits = GPU_UINT32_LENGTH - shiftedBits;
  return ((scalarValueB << shiftedBits)
        | (scalarValueA >> complementShiftedBits));
}

GPU_INLINE __device__ uint32_t funnelshift_lc_32(const uint32_t scalarValueA,
                                                 const uint32_t scalarValueB,
                                                 const uint32_t shiftedBits)
{
  #if (__CUDA_ARCH__ < GPU_CC_KEPLER_2G)
    return(__emulated_funnelshift_right_32(scalarValueA, scalarValueB, shiftedBits));
  #else
    return(__funnelshift_lc(scalarValueA, scalarValueB, shiftedBits));
  #endif
}

GPU_INLINE __device__ uint32_t gpu_get_thread_idx()
{
  #if (__CUDA_ARCH__ < GPU_CC_KEPLER_1G)
    return(blockIdx.y * blockDim.y + (blockIdx.x * blockDim.x + threadIdx.x));
  #else
    return(blockIdx.x * blockDim.x + threadIdx.x);
  #endif
}

GPU_INLINE __device__ uint32_t gpu_get_sm_idx(){
  uint32_t smid;
  asm volatile("mov.u32 %0, %%smid;" : "=r"(smid));
  return(smid);
}

GPU_INLINE __device__ uint32_t gpu_get_warp_idx(){
  uint32_t warpid;
  asm volatile("mov.u32 %0, %%warpid;" : "=r"(warpid));
  return(warpid);
}

GPU_INLINE __device__ uint32_t gpu_get_lane_idx(){
  uint32_t laneid;
  asm volatile("mov.u32 %0, %%laneid;" : "=r"(laneid));
  return(laneid);
}

#endif /* GPU_RESOURCES_H_ */
