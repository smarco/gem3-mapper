/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
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

#define GPU_NVCC_COLLECTIVES_REQ_SYNC 9000

#define GPU_CC_FERMI_1G   200
#define GPU_CC_FERMI_2G   210
#define GPU_CC_KEPLER_1G  300
#define GPU_CC_KEPLER_2G  350
#define GPU_CC_KEPLER_3G  370
#define GPU_CC_MAXWELL_1G 500
#define GPU_CC_MAXWELL_2G 520
#define GPU_CC_PASCAL_1G  610
#define GPU_CC_PASCAL_2G  600
#define GPU_CC_VOLTA_1G   720
#define GPU_CC_VOLTA_2G   700

/* Defines related to GPU Architecture */
#if   (__CUDA_ARCH__ <= GPU_CC_FERMI_2G)
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_FERMI
#elif (__CUDA_ARCH__ <= GPU_CC_KEPLER_3G)
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_KEPLER
#elif (__CUDA_ARCH__ <= GPU_CC_MAXWELL_2G)
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_MAXWELL
#elif (__CUDA_ARCH__ <= GPU_CC_PASCAL_2G)
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_PASCAL
#elif (__CUDA_ARCH__ <= GPU_CC_VOLTA_2G)
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_VOLTA
#else
  #define GPU_THREADS_PER_BLOCK   GPU_THREADS_PER_BLOCK_NEWGEN
#endif

// Setting threads in the warp that participate in a collective operation
#define GPU_COOPERATIVE_THREADS_ALL   0xFFFFFFFFu

/*************************************
GPU Side conversion types
**************************************/
typedef union{
  char4    v4;
  uint32_t s;
} gpu_char4_t;

/*************************************
GPU Side definitions (ASM instructions)
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

GPU_INLINE __device__ uint32_t ballot_32(const uint32_t threadCondition)
{
  #if (CUDART_VERSION < GPU_NVCC_COLLECTIVES_REQ_SYNC)
    //#warning ("Cooperative operations: Using native __ballot")
    return (__ballot(threadCondition));
  #else
    //#warning ("Cooperative operations: Using native sync __ballot")
    return (__ballot_sync(GPU_COOPERATIVE_THREADS_ALL, threadCondition));
  #endif
}

GPU_INLINE __device__ uint32_t any_32(const uint32_t threadCondition)
{
  #if (CUDART_VERSION < GPU_NVCC_COLLECTIVES_REQ_SYNC)
    //#warning ("Cooperative operations: Using native __any")
    return (__any(threadCondition));
  #else
    //#warning ("Cooperative operations: Using native sync __any")
    return (__any_sync(GPU_COOPERATIVE_THREADS_ALL, threadCondition));
  #endif
}

GPU_INLINE __device__ uint32_t shfl_32(uint32_t scalarValue, const int lane)
{
  #if (__CUDA_ARCH__ < GPU_CC_KEPLER_1G)
    //#warning ("Cooperative operations: Using emulated __shuffle")
    return (__emulated_shfl(scalarValue, (uint32_t)lane));
  #elif (CUDART_VERSION < GPU_NVCC_COLLECTIVES_REQ_SYNC)
    //#warning ("Cooperative operations: Using native __shuffle")
    return (__shfl((int)scalarValue, (int)lane, GPU_WARP_SIZE));
  #else
    //#warning ("Cooperative operations: Using native sync __shuffle")
    return (__shfl_sync(GPU_COOPERATIVE_THREADS_ALL, (int)scalarValue, (int)lane, GPU_WARP_SIZE));
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

GPU_INLINE __device__ uint64_t funnelshift_right_64(const uint64_t scalarValueA,
                                                    const uint64_t scalarValueB,
                                                    const uint32_t shiftedBits)
{   /* Concat A with B and right-shifts nbits:  (A | B) >> nbits */
  const uint32_t complementShiftedBits = GPU_UINT64_LENGTH - shiftedBits;
  return ( (scalarValueA >> shiftedBits)
         | (scalarValueB << complementShiftedBits));
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

// Compose
GPU_INLINE __device__ void gpu_decompose_uintv4(uint32_t* const arrayV4, const uint4 scalarV4)
{
	arrayV4[0] = scalarV4.x;
	arrayV4[1] = scalarV4.y;
	arrayV4[2] = scalarV4.z;
	arrayV4[3] = scalarV4.w;
}

// Decompose
GPU_INLINE __device__ uint4 gpu_compose_uintv4(const uint32_t* const arrayV4)
{
  uint4 scalarV4 = {arrayV4[0], arrayV4[1], arrayV4[2], arrayV4[3]};
  return(scalarV4);
}

// Extract
GPU_INLINE __device__ uint32_t gpu_extract_uintv4(const int32_t idElement, const uint32_t* const arrayV4)
{
  const uint32_t VECTOR_LENGHT = 4;
  uint32_t       value         = arrayV4[0];
  #pragma unroll
  for(uint32_t i = 1; i < VECTOR_LENGHT; ++i){
    const uint32_t tmp = arrayV4[i];
    value = (idElement == i) ? tmp : value;
  }
  return value;
}

// Select
GPU_INLINE __device__ uint32_t gpu_select_uintv4(const int32_t idElement, const uint4 scalarV4)
{
  uint32_t value = scalarV4.x;
  value = (idElement == 1) ? scalarV4.y : value;
  value = (idElement == 2) ? scalarV4.z : value;
  value = (idElement == 3) ? scalarV4.w : value;
  return value;
}

#endif /* GPU_RESOURCES_H_ */
