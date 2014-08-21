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

inline __device__ uint32_t selectEq_CC35(const uint32_t indexBase, 
				    					 const uint32_t Eq0, const uint32_t Eq1, 
				    					 const uint32_t Eq2, const uint32_t Eq3, 
				    					 const uint32_t Eq4)
{
	uint32_t Eq = Eq0;

	Eq = (indexBase == 1) ? Eq1 : Eq;
	Eq = (indexBase == 2) ? Eq2 : Eq;
	Eq = (indexBase == 3) ? Eq3 : Eq;
	Eq = (indexBase == 4) ? Eq4 : Eq;

	return Eq;
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

__device__ void myerslocalKeplerKernel_CC35( const d_qryEntry_t *d_queries, const uint32_t * __restrict d_reference, const candInfo_t *d_candidates, 
											 const uint32_t *d_reorderBuffer, resEntry_t *d_reorderResults, const qryInfo_t *d_qinfo, 
								 			 const uint32_t idCandidate, const uint32_t sizeRef, const uint32_t numReorderedResults, 
											 const uint32_t intraQueryThreadIdx, const uint32_t threadsPerQuery)
{
	if (idCandidate < numReorderedResults){

		const uint32_t * __restrict localCandidate;

		uint32_t Ph_A, Mh_A, Pv_A, Mv_A, Xv_A, Xh_A, Eq_A, tEq_A;
		uint32_t Ph_B, Mh_B, Pv_B, Mv_B, Xv_B, Xh_B, Eq_B, tEq_B;
		uint32_t Ph_C, Mh_C, Pv_C, Mv_C, Xv_C, Xh_C, Eq_C, tEq_C;
		uint32_t Ph_D, Mh_D, Pv_D, Mv_D, Xv_D, Xh_D, Eq_D, tEq_D;
		uint4 	 Eq0, Eq1, Eq2, Eq3, Eq4;
		uint32_t PH, MH, indexWord;
		uint32_t sum_A, sum_B, sum_C, sum_D;

		const uint32_t originalCandidate = d_reorderBuffer[idCandidate];
		const uint64_t positionRef = d_candidates[originalCandidate].position;
		const uint32_t sizeQuery = d_qinfo[d_candidates[originalCandidate].query].size;
		const uint32_t entry = d_qinfo[d_candidates[originalCandidate].query].posEntry + intraQueryThreadIdx;
		const uint32_t sizeCandidate = d_candidates[originalCandidate].size; /* sizeQuery * (1 + 2 * distance)*/
		const uint32_t numEntriesPerCandidate = (sizeCandidate / REFERENCE_CHARS_PER_ENTRY) + ((sizeCandidate % REFERENCE_CHARS_PER_ENTRY) ? 2 : 1);
		uint32_t candidate;

		const uint32_t mask = ((sizeQuery % UINT32_LENGTH) == 0) ? UINT32_ONE_LAST_MASK : 1 << ((sizeQuery % UINT32_LENGTH) - 1);
		int32_t  score = sizeQuery, minScore = sizeQuery;
		uint32_t idColumn = 0, minColumn = 0, indexBase;
		uint32_t intraBase, idEntry;
		
		indexWord = ((sizeQuery - 1) & (PEQ_LENGTH_PER_CUDA_THREAD - 1)) / UINT32_LENGTH;

		if((positionRef < sizeRef) && ((sizeRef - positionRef) > sizeCandidate)){

			localCandidate = d_reference + (positionRef / REFERENCE_CHARS_PER_ENTRY);

			Pv_A = UINT32_ONES;
			Mv_A = 0;

			Pv_B = UINT32_ONES;
			Mv_B = 0;

			Pv_C = UINT32_ONES;
			Mv_C = 0;

			Pv_D = UINT32_ONES;
			Mv_D = 0;

			Eq0 = d_queries[entry].bitmap[0];
			Eq1 = d_queries[entry].bitmap[1];
			Eq2 = d_queries[entry].bitmap[2];
			Eq3 = d_queries[entry].bitmap[3];
			Eq4 = d_queries[entry].bitmap[4];

			for(idEntry = 0; idEntry < numEntriesPerCandidate; idEntry++){

				candidate = localCandidate[idEntry];

				for(intraBase = 0; intraBase < REFERENCE_CHARS_PER_ENTRY; intraBase++){
					
					indexBase = candidate & 0x07;
					Eq_A = selectEq_CC35(indexBase, Eq0.x, Eq1.x, Eq2.x, Eq3.x, Eq4.x);
					Eq_B = selectEq_CC35(indexBase, Eq0.y, Eq1.y, Eq2.y, Eq3.y, Eq4.y);
					Eq_C = selectEq_CC35(indexBase, Eq0.z, Eq1.z, Eq2.z, Eq3.z, Eq4.z);
					Eq_D = selectEq_CC35(indexBase, Eq0.w, Eq1.w, Eq2.w, Eq3.w, Eq4.w);

					Xv_A = Eq_A | Mv_A;
					Xv_B = Eq_B | Mv_B;
					Xv_C = Eq_C | Mv_C;
					Xv_D = Eq_D | Mv_D;

					tEq_A = Eq_A & Pv_A;
					tEq_B = Eq_B & Pv_B;
					tEq_C = Eq_C & Pv_C;
					tEq_D = Eq_D & Pv_D;

					shuffle_collaborative_sum_CC35(tEq_A, tEq_B, tEq_C, tEq_D, Pv_A, Pv_B, Pv_C, Pv_D, 
											  	   intraQueryThreadIdx, 
											  	   &sum_A, &sum_B, &sum_C, &sum_D);

					Xh_A = (sum_A ^ Pv_A) | Eq_A;
					Xh_B = (sum_B ^ Pv_B) | Eq_B;
					Xh_C = (sum_C ^ Pv_C) | Eq_C;
					Xh_D = (sum_D ^ Pv_D) | Eq_D;

					Ph_A = Mv_A | ~(Xh_A | Pv_A);
					Ph_B = Mv_B | ~(Xh_B | Pv_B);
					Ph_C = Mv_C | ~(Xh_C | Pv_C);
					Ph_D = Mv_D | ~(Xh_D | Pv_D);

					Mh_A = Pv_A & Xh_A;
					Mh_B = Pv_B & Xh_B;
					Mh_C = Pv_C & Xh_C;
					Mh_D = Pv_D & Xh_D;

					PH = select_CC35(indexWord, Ph_A, Ph_B, Ph_C, Ph_D);
					MH = select_CC35(indexWord, Mh_A, Mh_B, Mh_C, Mh_D);
					score += (((PH & mask) != 0) - ((MH & mask) != 0));

					shuffle_collaborative_shift_CC35(Ph_A, Ph_B, Ph_C, Ph_D, 
													 intraQueryThreadIdx,
													 &Ph_A, &Ph_B, &Ph_C, &Ph_D);
					shuffle_collaborative_shift_CC35(Mh_A, Mh_B, Mh_C, Mh_D, 
													 intraQueryThreadIdx,
													 &Mh_A, &Mh_B, &Mh_C, &Mh_D);

					Pv_A = Mh_A | ~(Xv_A | Ph_A);
					Pv_B = Mh_B | ~(Xv_B | Ph_B);
					Pv_C = Mh_C | ~(Xv_C | Ph_C);
					Pv_D = Mh_D | ~(Xv_D | Ph_D);

					Mv_A = Ph_A & Xv_A;
					Mv_B = Ph_B & Xv_B;
					Mv_C = Ph_C & Xv_C;
					Mv_D = Ph_D & Xv_D;

					candidate >>= REFERENCE_CHAR_LENGTH;
					minColumn = (score < minScore) ? idColumn : minColumn;
					minScore  = (score < minScore) ? score    : minScore;
					if(intraQueryThreadIdx  == (threadsPerQuery - 1))
					idColumn++;
				}
			}

			if(intraQueryThreadIdx  == (threadsPerQuery - 1)){
	    		d_reorderResults[idCandidate].column = minColumn - (positionRef % REFERENCE_CHARS_PER_ENTRY);
	    		d_reorderResults[idCandidate].score = minScore;
			}
		}
	}
}

__global__ void myersKeplerKernel_CC35(const d_qryEntry_t *d_queries, const uint32_t * d_reference, const candInfo_t *d_candidates, const uint32_t *d_reorderBuffer, 
						    		   resEntry_t *d_reorderResults, const qryInfo_t *d_qinfo, const uint32_t sizeRef,  const uint32_t numReorderedResults,
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
