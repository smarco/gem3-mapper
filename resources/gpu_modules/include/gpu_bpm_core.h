
extern "C" {
#include "gpu_commons.h"
#include "gpu_reference.h"
#include "gpu_buffers.h"
}
#include "gpu_resources.h"


#ifndef GPU_BPM_CORE_H_
#define GPU_BPM_CORE_H_

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

	while(__any(carry)){
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

GPU_INLINE __device__ void set_BMP(uint32_t *BMP, const uint4 BMPv4)
{
	BMP[0] = BMPv4.x;
	BMP[1] = BMPv4.y;
	BMP[2] = BMPv4.z;
	BMP[3] = BMPv4.w;
}

inline __device__ uint32_t select(const int32_t indexWord, uint32_t *BMP, uint32_t value, const uint32_t BMPS_PER_THREAD)
{
	#pragma unroll
	for(uint32_t i = 0; i < BMPS_PER_THREAD; ++i){
		const uint32_t tmp = BMP[i];
		value = (indexWord == i) ? tmp : value;
	}
	return value;
}

#endif /* GPU_BPM_CORE_H_ */

