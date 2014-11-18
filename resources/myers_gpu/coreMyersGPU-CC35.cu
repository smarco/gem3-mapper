/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: BPM implementation for CUDA GPUs with compute capability 3.5 
 */

#include <stdio.h>
#include "myers-common.h"

inline __device__ void shuffle_collaborative_shift_CC35(uint32_t value_A, uint32_t value_B, uint32_t value_C, uint32_t value_D,  
					       						   		const uint32_t localThreadIdx, 
					       						   		uint32_t* res_A, uint32_t* res_B, uint32_t* res_C, uint32_t* res_D)
{
	uint32_t carry;

	carry = __shfl_up((int) value_D, 1);
	carry = (localThreadIdx) ? carry : 0;
	value_D = __funnelshift_lc(value_C, value_D, 1);
	value_C = __funnelshift_lc(value_B, value_C, 1);
	value_B = __funnelshift_lc(value_A, value_B, 1);
	value_A = __funnelshift_lc(carry,   value_A, 1);

	(* res_A) = value_A;
	(* res_B) = value_B;
	(* res_C) = value_C;
	(* res_D) = value_D;
}

inline __device__ void shuffle_collaborative_sum_CC35(const uint32_t a_A, const uint32_t a_B, const uint32_t a_C, const uint32_t a_D, 
					 							 	  const uint32_t b_A, const uint32_t b_B, const uint32_t b_C, const uint32_t b_D, 
					 							 	  const uint32_t localThreadIdx,
					 							 	  uint32_t* sum_A, uint32_t* sum_B, uint32_t* sum_C, uint32_t* sum_D)
{

	uint32_t carry, c_A, c_B, c_C, c_D;

	UADD__CARRY_OUT   (c_A, a_A, b_A)
	UADD__IN_CARRY_OUT(c_B, a_B, b_B)
	UADD__IN_CARRY_OUT(c_C, a_C, b_C)
	UADD__IN_CARRY_OUT(c_D, a_D, b_D)
	UADD__IN_CARRY    (carry, 0,   0)

	while(__any(carry)){
		carry = __shfl_up((int) (carry), 1);
		carry = (localThreadIdx) ? carry : 0;
		UADD__CARRY_OUT   (c_A, c_A, carry)
		UADD__IN_CARRY_OUT(c_B, c_B, 0)
		UADD__IN_CARRY_OUT(c_C, c_C, 0)
		UADD__IN_CARRY_OUT(c_D, c_D, 0)
		UADD__IN_CARRY    (carry, 0, 0) 
	}

	(* sum_A) = c_A;
	(* sum_B) = c_B;
	(* sum_C) = c_C;
	(* sum_D) = c_D;
}

inline __device__ uint32_t select_CC35(const uint32_t indexWord, 
				    				   const uint32_t A, const uint32_t B, 
				    				   const uint32_t C, const uint32_t D)
{
	uint32_t value = A;

	value = (indexWord == 1) ? B : value;
	value = (indexWord == 2) ? C : value;
	value = (indexWord == 3) ? D : value;

	return value;
}

inline __device__ uint64_t funnelShiftL_CC35(const uint64_t currentCandidateEntry, 
				    				   		const uint64_t lastCandidateEntry,
				    				   		const uint32_t shiftedBits)
{
	const uint32_t complementShiftedBits = BMP_GPU_UINT64_LENGTH - shiftedBits;
	
	return ((lastCandidateEntry >> shiftedBits) |
			(currentCandidateEntry <<  complementShiftedBits));
}

__device__ void myerslocalKeplerKernel_CC35( const d_qryEntry_t * __restrict d_queries, const uint64_t * __restrict d_reference, const bpm_gpu_cand_info_t *d_candidates,
											 const uint32_t *d_reorderBuffer, bpm_gpu_res_entry_t *d_reorderResults, const bpm_gpu_qry_info_t *d_qinfo,
								 			 const uint32_t idCandidate, const uint64_t sizeRef, const uint32_t numReorderedResults, 
											 const uint32_t intraQueryThreadIdx, const uint32_t threadsPerQuery)
{
	if (idCandidate < numReorderedResults){

		const uint64_t * __restrict localCandidate;

		uint4 	Ph, Mh, Pv, Mv, Xv, Xh, tEq, Eq, sum;
		uint32_t PH, MH, indexWord;

		const uint32_t originalCandidate = d_reorderBuffer[idCandidate];
		const uint64_t positionRef = d_candidates[originalCandidate].position;
		const uint32_t sizeQuery = d_qinfo[d_candidates[originalCandidate].query].size;
		const uint32_t entry = d_qinfo[d_candidates[originalCandidate].query].posEntry + intraQueryThreadIdx;
		const uint32_t sizeCandidate = d_candidates[originalCandidate].size;
		const uint32_t candidateAlignment = (positionRef % REFERENCE_CHARS_PER_ENTRY) * REFERENCE_CHAR_LENGTH;

		uint64_t candidate, lastCandidateEntry, currentCandidateEntry;

		const uint32_t mask = ((sizeQuery % BMP_GPU_UINT32_LENGTH) == 0) ? UINT32_ONE_LAST_MASK : 1 << ((sizeQuery % BMP_GPU_UINT32_LENGTH) - 1);
		int32_t  score = sizeQuery, minScore = sizeQuery;
		uint32_t idColumn = 0, minColumn = 0, idEntry = 0;
		
		indexWord = ((sizeQuery - 1) & (PEQ_LENGTH_PER_CUDA_THREAD - 1)) / BMP_GPU_UINT32_LENGTH;

		if((positionRef < sizeRef) && ((sizeRef - positionRef) > sizeCandidate)){

			localCandidate = d_reference + (positionRef / REFERENCE_CHARS_PER_ENTRY);

			Pv.x = UINT32_ONES;
			Mv.x = 0;

			Pv.y = UINT32_ONES;
			Mv.y = 0;

			Pv.z = UINT32_ONES;
			Mv.z = 0;

			Pv.w = UINT32_ONES;
			Mv.w = 0;

			lastCandidateEntry = localCandidate[idEntry];

			for(idColumn = 0; idColumn < sizeCandidate; idColumn++){

				if((idColumn % REFERENCE_CHARS_PER_ENTRY) == 0){
						idEntry++;
						currentCandidateEntry = localCandidate[idEntry];
						//candidate = __funnelshift_l(currentCandidateEntry, lastCandidateEntry, candidateAlignment); 
						candidate = funnelShiftL_CC35(currentCandidateEntry, lastCandidateEntry, candidateAlignment);
						lastCandidateEntry = currentCandidateEntry;
				}

				Eq = __ldg(&d_queries[entry].bitmap[candidate & 0x07]);

				Xv.x = Eq.x | Mv.x;
				Xv.y = Eq.y | Mv.y;
				Xv.z = Eq.z | Mv.z;
				Xv.w = Eq.w | Mv.w;

				tEq.x = Eq.x & Pv.x;
				tEq.y = Eq.y & Pv.y;
				tEq.z = Eq.z & Pv.z;
				tEq.w = Eq.w & Pv.w;

				//TODO: Review nvcc code generation using inline functions + param. by reference
				shuffle_collaborative_sum_CC35(tEq.x, tEq.y, tEq.z, tEq.w, Pv.x, Pv.y, Pv.z, Pv.w, 
										  	   intraQueryThreadIdx, 
										  	   &sum.x, &sum.y, &sum.z, &sum.w);

				Xh.x = (sum.x ^ Pv.x) | Eq.x;
				Xh.y = (sum.y ^ Pv.y) | Eq.y;
				Xh.z = (sum.z ^ Pv.z) | Eq.z;
				Xh.w = (sum.w ^ Pv.w) | Eq.w;

				Ph.x = Mv.x | ~(Xh.x | Pv.x);
				Ph.y = Mv.y | ~(Xh.y | Pv.y);
				Ph.z = Mv.z | ~(Xh.z | Pv.z);
				Ph.w = Mv.w | ~(Xh.w | Pv.w);

				Mh.x = Pv.x & Xh.x;
				Mh.y = Pv.y & Xh.y;
				Mh.z = Pv.z & Xh.z;
				Mh.w = Pv.w & Xh.w;

				PH = select_CC35(indexWord, Ph.x, Ph.y, Ph.z, Ph.w);
				MH = select_CC35(indexWord, Mh.x, Mh.y, Mh.z, Mh.w);
				score += (((PH & mask) != 0) - ((MH & mask) != 0));

				shuffle_collaborative_shift_CC35(Ph.x, Ph.y, Ph.z, Ph.w, 
												 intraQueryThreadIdx,
										   		 &Ph.x, &Ph.y, &Ph.z, &Ph.w);
				shuffle_collaborative_shift_CC35(Mh.x, Mh.y, Mh.z, Mh.w,
												 intraQueryThreadIdx,
										   		 &Mh.x, &Mh.y, &Mh.z, &Mh.w);

				Pv.x = Mh.x | ~(Xv.x | Ph.x);
				Pv.y = Mh.y | ~(Xv.y | Ph.y);
				Pv.z = Mh.z | ~(Xv.z | Ph.z);
				Pv.w = Mh.w | ~(Xv.w | Ph.w);

				Mv.x = Ph.x & Xv.x;
				Mv.y = Ph.y & Xv.y;
				Mv.z = Ph.z & Xv.z;
				Mv.w = Ph.w & Xv.w;

				candidate >>= REFERENCE_CHAR_LENGTH;
				minColumn = (score < minScore) ? idColumn : minColumn;
				minScore  = (score < minScore) ? score    : minScore;
			}

			if(intraQueryThreadIdx  == (threadsPerQuery - 1)){
				d_reorderResults[idCandidate].column = minColumn;
				d_reorderResults[idCandidate].score = minScore;
			}
		}
	}
}

__global__ void myersKeplerKernel_CC35(const d_qryEntry_t *d_queries, const uint64_t * d_reference, const bpm_gpu_cand_info_t *d_candidates, const uint32_t *d_reorderBuffer,
						    		   bpm_gpu_res_entry_t *d_reorderResults, const bpm_gpu_qry_info_t *d_qinfo, const uint64_t sizeRef,  const uint32_t numReorderedResults,
						    		   uint32_t *d_initPosPerBucket, uint32_t *d_initWarpPerBucket, uint32_t numWarps)
{
	uint32_t bucketIdx = 0;
	uint32_t globalThreadIdx = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t globalWarpIdx = globalThreadIdx / WARP_SIZE;
	uint32_t localThreadInTheBucket, idCandidate, intraQueryThreadIdx, threadsPerQuery, queriesPerWarp, localIdCandidateInTheBucket;

	while((bucketIdx != (WARP_SIZE + 1)) && (d_initWarpPerBucket[bucketIdx] <= globalWarpIdx)){
		bucketIdx++;
	}
	bucketIdx--;

	localThreadInTheBucket = globalThreadIdx - (d_initWarpPerBucket[bucketIdx] * WARP_SIZE);
	threadsPerQuery = bucketIdx + 1;
	queriesPerWarp = WARP_SIZE / threadsPerQuery;
	localIdCandidateInTheBucket = ((localThreadInTheBucket / WARP_SIZE) * queriesPerWarp) + ((threadIdx.x % WARP_SIZE) / threadsPerQuery);
	idCandidate = d_initPosPerBucket[bucketIdx] + localIdCandidateInTheBucket;
	intraQueryThreadIdx = (threadIdx.x % WARP_SIZE) % threadsPerQuery;

	myerslocalKeplerKernel_CC35(d_queries, d_reference, d_candidates, d_reorderBuffer, d_reorderResults, d_qinfo,
	 				 			idCandidate, sizeRef, numReorderedResults, intraQueryThreadIdx, threadsPerQuery);
}

extern "C"
myersError_t processMyersBufferOnKepler2ndGen(buffer_t *mBuff)
{
	reference_buffer_t 	*ref 		= mBuff->reference;
	queries_buffer_t 	*qry 		= mBuff->queries;
	candidates_buffer_t	*cand 		= mBuff->candidates;
	reorder_buffer_t 	*rebuff 	= mBuff->reorderBuffer;
	results_buffer_t 	*res 		= mBuff->results;
	cudaStream_t 		idStream	= mBuff->idStream;
	uint32_t			idSupDev	= mBuff->device->idSupportedDevice;


	uint32_t threadsPerBlock = CUDA_NUM_THREADS;
	uint32_t numThreads = rebuff->numWarps * WARP_SIZE;
	uint32_t blocksPerGrid = DIV_CEIL(numThreads, threadsPerBlock);

	//printf("KEPLER 2ndGen: LAUNCH KERNEL -- Bloques: %d - Th_block %d\n", blocksPerGrid, threadsPerBlock);
	myersKeplerKernel_CC35<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((d_qryEntry_t *)qry->d_queries, ref->d_reference[idSupDev], cand->d_candidates, rebuff->d_reorderBuffer,
																		res->d_reorderResults, qry->d_qinfo, ref->size, res->numReorderedResults,
																		rebuff->d_initPosPerBucket, rebuff->d_initWarpPerBucket, rebuff->numWarps);
	return(SUCCESS);
}
