
#ifndef GPU_BPM_CORE_H_
#define GPU_BPM_CORE_H_

#include "gpu_resources.h"
#include "gpu_commons.h"
#include "gpu_scheduler.h"

GPU_INLINE __device__ void shuffle_collaborative_shift(uint32_t *value, const uint32_t shiftedBits,
													   const uint32_t localThreadIdx, const uint32_t BMPS_PER_THREAD)
{
	uint32_t carry;

	carry = __shfl_up((int) value[BMPS_PER_THREAD-1], 1);
	carry = (localThreadIdx) ? carry : 0;
	#pragma unroll
	for(uint32_t idBMP = BMPS_PER_THREAD-1; idBMP > 0; --idBMP){
		value[idBMP] = __funnelshift_lc(value[idBMP-1], value[idBMP], shiftedBits);
	}
	value[0] = __funnelshift_lc(carry, value[0], shiftedBits);
}

GPU_INLINE __device__ void shuffle_collaborative_sum(const uint32_t *A, const uint32_t *B, uint32_t *C,
												 const uint32_t localThreadIdx, const uint32_t BMPS_PER_THREAD)

{
	uint32_t carry;

	UADD__CARRY_OUT(C[0], A[0], B[0])
	#pragma unroll
	for(uint32_t idBMP = 1; idBMP < BMPS_PER_THREAD; ++idBMP){
		UADD__IN_CARRY_OUT(C[idBMP], A[idBMP], B[idBMP])
	}
	UADD__IN_CARRY(carry, 0, 0)

	while(__any(carry)){
		carry = __shfl_up((int) (carry), 1);
		carry = (localThreadIdx) ? carry : 0;
		UADD__CARRY_OUT(C[0], C[0], carry)
		#pragma unroll
		for(uint32_t idBMP = 1; idBMP < BMPS_PER_THREAD; ++idBMP){
			UADD__IN_CARRY_OUT(C[idBMP], C[idBMP], 0)
		}
		UADD__IN_CARRY(carry, 0, 0)
	}
}

GPU_INLINE __device__ void setBMP(uint32_t *BMP, const uint4 BMPv4)
{
	BMP[0] = BMPv4.x;
	BMP[1] = BMPv4.y;
	BMP[2] = BMPv4.z;
	BMP[3] = BMPv4.w;
}

GPU_INLINE __device__ uint64_t funnelShiftL(const uint64_t currentCandidateEntry,
				    				   		 const uint64_t lastCandidateEntry,
				    				   		 const uint32_t shiftedBits)
{
	const uint32_t complementShiftedBits = GPU_UINT64_LENGTH - shiftedBits;

	return ((lastCandidateEntry >> shiftedBits) |
			(currentCandidateEntry <<  complementShiftedBits));
}

#endif /* GPU_BPM_CORE_H_ */

