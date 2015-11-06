
#ifndef GPU_RESOURCES_H_
#define GPU_RESOURCES_H_

#include "gpu_commons.h"

/*************************************
GPU Side defines (ASM instructions)
**************************************/

// output temporal carry in internal register
#define UADD__CARRY_OUT(c, a, b) \
    asm volatile("add.cc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b));

// add & output with temporal carry of internal register
#define UADD__IN_CARRY_OUT(c, a, b) \
    asm volatile("addc.cc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b));

// add with temporal carry of internal register
#define UADD__IN_CARRY(c, a, b) \
    asm volatile("addc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b));

// packing and unpacking: from uint64_t to uint2
#define V2S_B64(v,s) \
	asm("mov.b64 %0, {%1,%2};" : "=l"(s) : "r"(v.x), "r"(v.y))

// packing and unpacking: from uint2 to uint64_t
#define S2V_B64(s,v) \
	asm("mov.b64 {%0,%1}, %2;" : "=r"(v.x), "=r"(v.y) : "l"(s))


/*************************************
DEVICE side basic block primitives
**************************************/

inline __device__ uint64_t shfl_64(uint64_t scalarValue, const int lane)
{
	uint2 vectorValue;
	S2V_B64(scalarValue, vectorValue);
	vectorValue.x = __shfl(vectorValue.x, lane);
	vectorValue.y = __shfl(vectorValue.y, lane);
	V2S_B64(vectorValue, scalarValue);
	return (scalarValue);
}

inline __device__ uint32_t shfl_32(uint32_t scalarValue, const int lane)
{
	return (__shfl(scalarValue, lane));
}

inline __device__ uint32_t reduceAddWarp(uint32_t resultBitmaps, const uint32_t numThreads)
{
	for (int32_t i = 1; i < numThreads; i *= 2){
		int32_t n = __shfl_down((int) resultBitmaps, i, GPU_WARP_SIZE);
			resultBitmaps += n;
	}
	return(resultBitmaps);
}

#endif /* GPU_RESOURCES_H_ */
