
#ifndef GPU_RESOURCES_H_
#define GPU_RESOURCES_H_

#include "gpu_commons.h"
#include "gpu_devices.h"

/* Defines related to GPU Architecture */
#if   (__CUDA_ARCH__ < GPU_CC_KEPLER_1G)
	#define	GPU_THREADS_PER_BLOCK		256
#elif (__CUDA_ARCH__ < GPU_CC_MAXWELL_1G)
	#define	GPU_THREADS_PER_BLOCK		128
#else
	#define	GPU_THREADS_PER_BLOCK		64
#endif

#define GPU_CC_FERMI_1G		200
#define GPU_CC_FERMI_2G		210
#define GPU_CC_KEPLER_1G	300
#define GPU_CC_KEPLER_2G	350
#define GPU_CC_MAXWELL_1G	500
#define GPU_CC_MAXWELL_2G	520

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
		return (__shfl(scalarValue, lane));
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
	return ( (scalarValueA <<  complementShiftedBits)
			| (scalarValueB >> shiftedBits));
}

GPU_INLINE __device__ uint32_t __emulated_funnelshift_right_32(const uint32_t scalarValueA,
				    				   			   	   	   	   const uint32_t scalarValueB,
				    				   			   	   	   	   const uint32_t shiftedBits)
{	/* (scalarValueB | scalarValueA) >> 1 : (stores higher 32 bits) */
	const uint32_t complementShiftedBits = GPU_UINT32_LENGTH - shiftedBits;
	return ((scalarValueB <<  shiftedBits)
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

GPU_INLINE void gpu_kernel_thread_configuration(const uint32_t numThreads, dim3 *blocksPerGrid, dim3 *threadsPerBlock)
{
	#if (__CUDA_ARCH__ < GPU_CC_KEPLER_1G)
		//We use 2-Dimensional Grid (because Fermi is limited to 65535 Blocks per dim)
		const uint32_t maxBlocksPerRow = 65535;
		const uint32_t numBlocks = GPU_DIV_CEIL(numThreads, GPU_THREADS_PER_BLOCK);
		const uint32_t rowsPerGrid = GPU_DIV_CEIL(numBlocks, maxBlocksPerRow);
		const uint32_t blocksPerRow = (rowsPerGrid > 1) ? maxBlocksPerRow : numBlocks;

		blocksPerGrid->x   = blocksPerRow;
		blocksPerGrid->y   = rowsPerGrid;
		threadsPerBlock->x = GPU_THREADS_PER_BLOCK;
	#else
		threadsPerBlock->x = GPU_THREADS_PER_BLOCK;
		blocksPerGrid->x   = GPU_DIV_CEIL(numThreads, GPU_THREADS_PER_BLOCK);
	#endif
}

#endif /* GPU_RESOURCES_H_ */
