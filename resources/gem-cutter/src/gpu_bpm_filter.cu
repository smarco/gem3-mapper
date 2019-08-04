/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_BPM_FILTER_CU_
#define GPU_BPM_FILTER_CU_

#include "../include/gpu_bpm_core.h"
#include "../include/gpu_text_core.h"


GPU_INLINE __device__ void gpu_bpm_filter_local_kernel_4bases(const gpu_bpm_filter_device_qry_entry_t * __restrict d_queries, const uint64_t * __restrict d_reference,
		                                                          const uint64_t* __restrict d_reference_masked, const gpu_bpm_filter_cand_info_t *d_candidates, const uint32_t *d_reorderBuffer,
                                                              gpu_bpm_filter_alg_entry_t* const d_results, const gpu_bpm_filter_qry_info_t *d_qinfo,
                                                              const uint32_t idCandidate, const uint64_t sizeRef, const uint32_t numResults,
                                                              const uint32_t intraQueryThreadIdx, const uint32_t threadsPerQuery, const bool binning)
{
  if (idCandidate < numResults){
    const uint64_t *localCandidate;
    const uint32_t PEQS_PER_THREAD = 1;
    const uint32_t BMPS_SIZE       = GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD * PEQS_PER_THREAD;
    const uint32_t BMPS_PER_THREAD = BMPS_SIZE / GPU_UINT32_LENGTH;

    const uint32_t originalCandidate  = (binning) ? d_reorderBuffer[idCandidate] : idCandidate;
    const uint64_t positionRef        = d_candidates[originalCandidate].position;
    const uint32_t sizeQuery          = d_qinfo[d_candidates[originalCandidate].query].tileSize;
    const uint32_t entry              = d_qinfo[d_candidates[originalCandidate].query].posEntry + (intraQueryThreadIdx * PEQS_PER_THREAD);
    const uint32_t sizeCandidate      = d_candidates[originalCandidate].size;
    const uint32_t candidateAlignment = (positionRef % GPU_REFERENCE_CHARS_PER_ENTRY) * GPU_REFERENCE_CHAR_LENGTH;

    uint64_t candidate, lastCandidateEntry, currentCandidateEntry;

    const int32_t indexWord = ((sizeQuery - 1) % BMPS_SIZE) / GPU_UINT32_LENGTH;
    const uint32_t mask = ((sizeQuery % GPU_UINT32_LENGTH) == 0) ? GPU_UINT32_MASK_ONE_HIGH : 1 << ((sizeQuery % GPU_UINT32_LENGTH) - 1);
    int32_t  score = sizeQuery, minScore = sizeQuery;
    uint32_t idColumn = 0, minColumn = 0, idEntry = 0;

    if((positionRef < sizeRef) && ((sizeRef - positionRef) > sizeCandidate)){

      uint32_t  Ph[BMPS_PER_THREAD], Mh[BMPS_PER_THREAD],  Pv[BMPS_PER_THREAD], Mv[BMPS_PER_THREAD];
      uint32_t  Xv[BMPS_PER_THREAD], Xh[BMPS_PER_THREAD], tEq[BMPS_PER_THREAD], Eq[BMPS_PER_THREAD];
      uint32_t sum[BMPS_PER_THREAD];

      localCandidate = d_reference + (positionRef / GPU_REFERENCE_CHARS_PER_ENTRY);

      #pragma unroll
      for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP){
        Pv[idBMP] = GPU_UINT32_ONES;
        Mv[idBMP] = 0;
      }

      lastCandidateEntry = localCandidate[idEntry];

      for(idColumn = 0; idColumn < sizeCandidate; idColumn++){
        uint32_t PH, MH;

        if((idColumn % GPU_REFERENCE_CHARS_PER_ENTRY) == 0){
          idEntry++;
          currentCandidateEntry = localCandidate[idEntry];
          candidate = funnelshift_left_64(currentCandidateEntry, lastCandidateEntry, candidateAlignment);
          lastCandidateEntry = currentCandidateEntry;
        }

        #pragma unroll
        for(uint32_t idPEQ = 0, idBMP = 0; idPEQ < PEQS_PER_THREAD; ++idPEQ, idBMP+=4){
          uint4 Eqv4 = LDG(&d_queries[entry + idPEQ].bitmap[candidate & GPU_REFERENCE_MASK_BASE]);
          gpu_decompose_uintv4(Eq + idBMP, Eqv4);
        }

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Xv[idBMP] = Eq[idBMP] | Mv[idBMP];

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          tEq[idBMP] = Eq[idBMP] & Pv[idBMP];

        cooperative_sum(tEq, Pv, sum, intraQueryThreadIdx, BMPS_PER_THREAD);

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Xh[idBMP] = (sum[idBMP] ^ Pv[idBMP]) | Eq[idBMP];

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Ph[idBMP] = Mv[idBMP] | ~(Xh[idBMP] | Pv[idBMP]);

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Mh[idBMP] = Pv[idBMP] & Xh[idBMP];

        PH = gpu_extract_uintv4(indexWord, Ph);
        MH = gpu_extract_uintv4(indexWord, Mh);
        score += (((PH & mask) != 0) - ((MH & mask) != 0));

        cooperative_shift(Ph, 1, intraQueryThreadIdx, BMPS_PER_THREAD);
        cooperative_shift(Mh, 1, intraQueryThreadIdx, BMPS_PER_THREAD);

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Pv[idBMP] = Mh[idBMP] | ~(Xv[idBMP] | Ph[idBMP]);

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Mv[idBMP] = Ph[idBMP] & Xv[idBMP];

        candidate >>= GPU_REFERENCE_CHAR_LENGTH;
        minColumn = (score < minScore) ? idColumn : minColumn;
        minScore  = (score < minScore) ? score    : minScore;
      }

      if(intraQueryThreadIdx  == (threadsPerQuery - 1)){
        d_results[idCandidate].column = minColumn;
        d_results[idCandidate].score  = minScore;
      }
    }
  }
}

GPU_INLINE __device__ void gpu_bpm_filter_local_kernel(const gpu_bpm_filter_device_qry_entry_t * __restrict d_queries, const uint64_t * __restrict referencePlain,
		                                                   const uint64_t* __restrict referenceMasked, const gpu_bpm_filter_cand_info_t *d_candidates, const uint32_t *d_reorderBuffer,
                                                       gpu_bpm_filter_alg_entry_t* const d_results, const gpu_bpm_filter_qry_info_t *d_qinfo,
                                                       const uint32_t idCandidate, const uint64_t sizeRef, const uint32_t numResults,
                                                       const uint32_t intraQueryThreadIdx, const uint32_t threadsPerQuery, const bool binning)
{
  if (idCandidate < numResults){
    const uint32_t PEQS_PER_THREAD = 1;
    const uint32_t BMPS_SIZE       = GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD * PEQS_PER_THREAD;
    const uint32_t BMPS_PER_THREAD = BMPS_SIZE / GPU_UINT32_LENGTH;

    const uint32_t originalCandidate  = (binning) ? d_reorderBuffer[idCandidate] : idCandidate;
    const uint64_t posCandidate       = d_candidates[originalCandidate].position;
    const uint32_t sizeQuery          = d_qinfo[d_candidates[originalCandidate].query].tileSize;
    const uint32_t entry              = d_qinfo[d_candidates[originalCandidate].query].posEntry + (intraQueryThreadIdx * PEQS_PER_THREAD);
    const uint32_t sizeCandidate      = d_candidates[originalCandidate].size;

    ulong2   infoCandidatePlain  = GPU_TEXT_INIT, infoCandidateMasked = GPU_TEXT_INIT;
    const int32_t indexWord = ((sizeQuery - 1) % BMPS_SIZE) / GPU_UINT32_LENGTH;
    const uint32_t mask = ((sizeQuery % GPU_UINT32_LENGTH) == 0) ? GPU_UINT32_MASK_ONE_HIGH : 1 << ((sizeQuery % GPU_UINT32_LENGTH) - 1);
    int32_t  score = sizeQuery, minScore = sizeQuery;
    uint32_t idColumn = 0, minColumn = 0;

    if((posCandidate < sizeRef) && ((sizeRef - posCandidate) > sizeCandidate)){

      uint32_t  Ph[BMPS_PER_THREAD], Mh[BMPS_PER_THREAD],  Pv[BMPS_PER_THREAD], Mv[BMPS_PER_THREAD];
      uint32_t  Xv[BMPS_PER_THREAD], Xh[BMPS_PER_THREAD], tEq[BMPS_PER_THREAD], Eq[BMPS_PER_THREAD];
      uint32_t sum[BMPS_PER_THREAD];

      #pragma unroll
      for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP){
        Pv[idBMP] = GPU_UINT32_ONES;
        Mv[idBMP] = 0;
      }

      for(idColumn = 0; idColumn < sizeCandidate; idColumn++){
        uint32_t PH, MH;
        const uint8_t encBasePlain  = gpu_text_lookup(referencePlain, posCandidate + idColumn, &infoCandidatePlain, GPU_REFERENCE_PLAIN__CHAR_LENGTH);
        const uint8_t encBaseMasked = gpu_text_lookup(referenceMasked, posCandidate + idColumn, &infoCandidateMasked, GPU_REFERENCE_MASKED__CHAR_LENGTH);
        const uint8_t encBase       = (encBaseMasked << GPU_REFERENCE_PLAIN__CHAR_LENGTH) | encBasePlain;

        #pragma unroll
        for(uint32_t idPEQ = 0, idBMP = 0; idPEQ < PEQS_PER_THREAD; ++idPEQ, idBMP+=4){
          uint4 Eqv4 = LDG(&d_queries[entry + idPEQ].bitmap[encBase]);
          gpu_decompose_uintv4(Eq + idBMP, Eqv4);
        }

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Xv[idBMP] = Eq[idBMP] | Mv[idBMP];

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          tEq[idBMP] = Eq[idBMP] & Pv[idBMP];

        cooperative_sum(tEq, Pv, sum, intraQueryThreadIdx, BMPS_PER_THREAD);

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Xh[idBMP] = (sum[idBMP] ^ Pv[idBMP]) | Eq[idBMP];

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Ph[idBMP] = Mv[idBMP] | ~(Xh[idBMP] | Pv[idBMP]);

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Mh[idBMP] = Pv[idBMP] & Xh[idBMP];

        PH = gpu_extract_uintv4(indexWord, Ph);
        MH = gpu_extract_uintv4(indexWord, Mh);
        score += (((PH & mask) != 0) - ((MH & mask) != 0));

        cooperative_shift(Ph, 1, intraQueryThreadIdx, BMPS_PER_THREAD);
        cooperative_shift(Mh, 1, intraQueryThreadIdx, BMPS_PER_THREAD);

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Pv[idBMP] = Mh[idBMP] | ~(Xv[idBMP] | Ph[idBMP]);

        #pragma unroll
        for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP)
          Mv[idBMP] = Ph[idBMP] & Xv[idBMP];

        minColumn = (score < minScore) ? idColumn : minColumn;
        minScore  = (score < minScore) ? score    : minScore;
      }

      if(intraQueryThreadIdx  == (threadsPerQuery - 1)){
        d_results[idCandidate].column = minColumn;
        d_results[idCandidate].score  = minScore;
      }
    }
  }
}

__global__ void gpu_bpm_filter_kernel(const gpu_bpm_filter_device_qry_entry_t* const d_queries, const uint64_t* const d_reference, const uint64_t* const d_reference_masked, const gpu_bpm_filter_cand_info_t* const d_candidates,
                                      const uint32_t* const d_reorderBuffer, gpu_bpm_filter_alg_entry_t* const d_reorderResults, const gpu_bpm_filter_qry_info_t* const d_qinfo,
                                      const uint64_t sizeRef,  const uint32_t numResults, const uint32_t* const d_initPosPerBucket, const uint32_t* const d_initWarpPerBucket,
                                      const uint32_t numWarps, const bool binning)
{
  const uint32_t globalThreadIdx = gpu_get_thread_idx();
  const uint32_t globalWarpIdx   = globalThreadIdx / GPU_WARP_SIZE;

  uint32_t bucketIdx = 0;
  uint32_t localThreadInTheBucket, idCandidate, intraQueryThreadIdx, threadsPerQuery, queriesPerWarp, localIdCandidateInTheBucket;

  //Scan in which bucket is matched this warp
  while((bucketIdx != (GPU_WARP_SIZE + 1)) && (d_initWarpPerBucket[bucketIdx] <= globalWarpIdx)){
    bucketIdx++;
  }
  bucketIdx--;

  threadsPerQuery             = bucketIdx + 1;
  queriesPerWarp              = GPU_WARP_SIZE / threadsPerQuery;
  localThreadInTheBucket      = globalThreadIdx - (d_initWarpPerBucket[bucketIdx] * GPU_WARP_SIZE);
  localIdCandidateInTheBucket = ((localThreadInTheBucket / GPU_WARP_SIZE) * queriesPerWarp) + ((threadIdx.x % GPU_WARP_SIZE) / threadsPerQuery);
  idCandidate                 = d_initPosPerBucket[bucketIdx] + localIdCandidateInTheBucket;
  intraQueryThreadIdx         = (threadIdx.x % GPU_WARP_SIZE) % threadsPerQuery;

  gpu_bpm_filter_local_kernel(d_queries, d_reference, d_reference_masked, d_candidates, d_reorderBuffer, d_reorderResults, d_qinfo,
                              idCandidate, sizeRef, numResults, intraQueryThreadIdx, threadsPerQuery, binning);
}

extern "C"
gpu_error_t gpu_bpm_filter_process_buffer(gpu_buffer_t *mBuff)
{
  const gpu_reference_buffer_t* const             ref           =  mBuff->reference;
  const gpu_bpm_filter_queries_buffer_t* const    qry           = &mBuff->data.fbpm.queries;
  const gpu_bpm_filter_candidates_buffer_t* const cand          = &mBuff->data.fbpm.candidates;
  const gpu_scheduler_buffer_t* const             rebuff        = &mBuff->data.fbpm.reorderBuffer;
  const gpu_bpm_filter_alignments_buffer_t* const res           = &mBuff->data.fbpm.alignments;
  const cudaStream_t                              idStream      =  mBuff->listStreams[mBuff->idStream];
  const uint32_t                                  idSupDev      =  mBuff->idSupportedDevice;
  const gpu_device_info_t* const                  device        =  mBuff->device[idSupDev];
  const uint32_t                                  numAlignments =  res->numAlignments;
  const uint32_t                                  maxCandidates =  mBuff->data.fbpm.maxCandidates;
  const uint32_t                                  maxAlignments =  mBuff->data.fbpm.maxAlignments;
  const uint32_t                                  numResults    = (mBuff->data.fbpm.queryBinning) ? res->numReorderedAlignments : res->numAlignments;
  gpu_bpm_filter_alg_entry_t* const               d_results     = (mBuff->data.fbpm.queryBinning) ? res->d_reorderAlignments : res->d_alignments;
  // Work scheduling for the GPU kernel
  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = rebuff->numWarps * GPU_WARP_SIZE;
  gpu_device_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);
  // Sanity-check (checks buffer overflowing)
  if((numAlignments > maxCandidates) || (numAlignments > maxAlignments)) 
    return(E_OVERFLOWING_BUFFER);
  // Kernel Launcher
  gpu_bpm_filter_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((gpu_bpm_filter_device_qry_entry_t *)qry->d_queries, ref->d_reference_plain[idSupDev], ref->d_reference_masked[idSupDev],
                                                                          cand->d_candidates, rebuff->threadMapScheduler.d_reorderBuffer, d_results,
                                                                          qry->d_qinfo, ref->size, numResults,
                                                                          rebuff->d_initPosPerBucket, rebuff->d_initWarpPerBucket,
                                                                          rebuff->numWarps, mBuff->data.fbpm.queryBinning);
  return(SUCCESS);
}

#endif /* GPU_BPM_FILTER_CU_ */
