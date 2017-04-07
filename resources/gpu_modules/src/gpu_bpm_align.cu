/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_BPM_ALIGN_CU_
#define GPU_BPM_ALIGN_CU_

#include "../include/gpu_bpm_core.h"
#include "../include/gpu_text_core.h"
#include "../include/gpu_scheduler_core.h"

#define GPU_BMP_ALIGN_CIGAR_LUT_OFFSET      8
#define GPU_BMP_ALIGN_CIGAR_LUT_SIZE        16
#define GPU_BMP_ALIGN_BASE_CANDIDATE_LENGTH 8
#define GPU_BMP_ALIGN_BASE_QUERY_LENGTH     8


GPU_INLINE __device__ void gpu_bpm_align_backtrace(const uint32_t* const dpPV, const uint32_t* const dpMV,
												   const uint64_t* const query, const uint64_t* const candidate, gpu_bpm_align_cigar_entry_t* const dpCIGAR,
												   const bool leftGapAlign, const uint32_t minColumn, const uint32_t sizeQuery,
												   const uint32_t intraQueryThreadIdx, const uint32_t threadsPerQuery,
												   gpu_bpm_align_coord_t* const initCoodRes, uint32_t* const matchEffLenghtRes, uint32_t* const cigarLenghtRes)
{
  // Initializing back-trace threading variables
  const uint32_t threadColumnEntries  = GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH / GPU_UINT32_LENGTH;
  const uint32_t offsetQueryThreadIdx = gpu_get_lane_idx() - intraQueryThreadIdx;
  // Initializing back-trace data variables
  ulong2   infoCandidate = GPU_TEXT_INIT, infoQuery = GPU_TEXT_INIT;
  int32_t x = minColumn, y = sizeQuery - 1, matchEffLenght = 0;
  // Decomposing the CIGAR LUT
  const char4* const globalCigarTableLeft  = (char4*)  gpu_bmp_align_cigar_lut;
  const char4* const globalCigarTableRight = (char4*) (gpu_bmp_align_cigar_lut + GPU_BMP_ALIGN_CIGAR_LUT_OFFSET);
  // Initialization for CIGAR back-trace iteration variables
  gpu_char4_t cigarOP  = GPU_CIGAR_INIT;
  gpu_bpm_align_cigar_event_t accEvent = GPU_CIGAR_NULL, event = GPU_CIGAR_NULL;
  uint32_t cigarLenght = 0, accNum = 0;
  // Performing the back-trace to extract the cigar string
  while ((y >= 0) && (x >= 0)){
    // Thread managing
    const uint32_t dpActiveThread     =  y / GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH;
    const uint32_t dpLocalThreadEntry = (y % GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH) / GPU_UINT32_LENGTH;
    if(intraQueryThreadIdx == dpActiveThread){
      // Read query and candidate bases
      const uint8_t  encBaseCandidate = gpu_text_lookup(candidate, x, &infoCandidate, GPU_BMP_ALIGN_BASE_CANDIDATE_LENGTH);
      const uint8_t  encBaseQuery     = gpu_text_lookup(query, y, &infoQuery, GPU_BMP_ALIGN_BASE_QUERY_LENGTH);
      // Indexation for the dpMatrix element
      const uint32_t idBMP   = ((x + 1) * threadColumnEntries) + dpLocalThreadEntry;
      const uint32_t maskBMP = GPU_UINT32_ONE_MASK << (y % GPU_UINT32_LENGTH);
      // Select CIGAR operation on LUT
      const uint32_t deletion  = ((dpPV[idBMP] & maskBMP) != 0) << 2;
      const uint32_t insertion = ((dpMV[(idBMP - threadColumnEntries)] & maskBMP) != 0) << 1;
      const uint32_t match = (encBaseCandidate == encBaseQuery) && (encBaseQuery != GPU_ENC_DNA_CHAR_N);
      const uint32_t cigarEventEntry = deletion | insertion | match;
      // Recover CIGAR from DP matrix
      const char4* const localCigarTable = (leftGapAlign) ? globalCigarTableLeft : globalCigarTableRight;
      cigarOP.v4 = LDG(&localCigarTable[cigarEventEntry]);
    }
    // Communicate OP variable to the rest of the group threads
    cigarOP.s = shfl_32(cigarOP.s, offsetQueryThreadIdx + dpActiveThread);
    x += cigarOP.v4.x; y += cigarOP.v4.y; event = cigarOP.v4.z; matchEffLenght += cigarOP.v4.w;
    accNum++;
    // Save CIGAR string from end to start position & Resetting the CIGAR stats
    if((intraQueryThreadIdx == 0) && (event != accEvent) && (event != GPU_CIGAR_NULL)){
      const gpu_bpm_align_cigar_entry_t cigarEntry = {accEvent, accNum};
      dpCIGAR[sizeQuery - cigarLenght] = cigarEntry;
      accEvent = event; event = GPU_CIGAR_NULL; accNum = 0; cigarLenght++;
    }
  }
  // Master thread saves the last part of the cigar
  if (intraQueryThreadIdx  == (threadsPerQuery - 1)){
  const gpu_bpm_align_coord_t initCood =  {x, y};
	// Saving the last CIGAR status event
	if(event != GPU_CIGAR_NULL){
    const gpu_bpm_align_cigar_entry_t cigarEntry = {event, accNum};
	  dpCIGAR[sizeQuery - cigarLenght] = cigarEntry;
	  cigarLenght++;
	}
    // Saving the remainder semi-global deletion events
    if(y >= 0){
      const uint32_t numEvents = y + 1;
      const gpu_bpm_align_cigar_entry_t cigarEntry = {GPU_CIGAR_DELETION, numEvents};
	    dpCIGAR[sizeQuery - cigarLenght] = cigarEntry;
	    matchEffLenght -= numEvents;
    }
    // Returning back-trace results
    (* initCoodRes)       = initCood;
    (* matchEffLenghtRes) = matchEffLenght;
    (* cigarLenghtRes)    = cigarLenght;
  }
}

GPU_INLINE __device__ void gpu_bpm_align_dp_matrix(uint4* const dpPV, uint4* const dpMV,
												                           const gpu_bpm_align_device_qry_entry_t* const PEQs, const uint64_t* const candidate,
												                           const uint32_t sizeQuery, const uint32_t sizeCandidate,
												                           const uint32_t intraQueryThreadIdx, const uint32_t threadsPerQuery,
												                           uint32_t* const endMinColumn, uint32_t* const endMinScore)
{
  const uint32_t BMPS_SIZE       = GPU_BPM_ALIGN_PEQ_LENGTH_PER_CUDA_THREAD;
  const uint32_t BMPS_PER_THREAD = BMPS_SIZE / GPU_UINT32_LENGTH;

  const uint32_t laneIdx   = gpu_get_lane_idx();
  const int32_t  indexWord = ((sizeQuery - 1) % BMPS_SIZE) / GPU_UINT32_LENGTH;
  const uint32_t mask      = ((sizeQuery % GPU_UINT32_LENGTH) == 0) ? GPU_UINT32_MASK_ONE_HIGH : 1 << ((sizeQuery % GPU_UINT32_LENGTH) - 1);

  ulong2   infoCandidate = GPU_TEXT_INIT;
  int32_t  score = sizeQuery, minScore = sizeQuery;
  uint32_t idColumn = 0, minColumn = 0;

  if(GPU_BPM_ALIGN_MAX_SIZE_CANDIDATE > sizeCandidate){

    uint32_t Ph[BMPS_PER_THREAD], Mh[BMPS_PER_THREAD],  Pv[BMPS_PER_THREAD], Mv[BMPS_PER_THREAD];
    uint32_t Xv[BMPS_PER_THREAD], Xh[BMPS_PER_THREAD], tEq[BMPS_PER_THREAD], Eq[BMPS_PER_THREAD];
    uint32_t sum[BMPS_PER_THREAD];

    #pragma unroll
    for(uint32_t idBMP = 0; idBMP < BMPS_PER_THREAD; ++idBMP){
      Pv[idBMP] = GPU_UINT32_ONES;
      Mv[idBMP] = 0;
    }

    for(idColumn = 0; idColumn < sizeCandidate; idColumn++){
      uint32_t PH, MH;
      const uint8_t encBase = gpu_text_lookup(candidate, idColumn, &infoCandidate, GPU_BMP_ALIGN_BASE_CANDIDATE_LENGTH);
      const uint4 Eqv4 = LDG(&PEQs->bitmap[encBase]);
      gpu_decompose_uintv4(Eq, Eqv4);

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

      dpPV[idColumn] = gpu_compose_uintv4(Pv);
      dpMV[idColumn] = gpu_compose_uintv4(Mv);

      minColumn = (score < minScore) ? idColumn : minColumn;
      minScore  = (score < minScore) ? score    : minScore;
    }
    //Communicate the score and minimum values along thread group
    (* endMinColumn) = shfl_32(minColumn, laneIdx + (threadsPerQuery - intraQueryThreadIdx - 1));
    (* endMinScore)  = shfl_32(minScore,  laneIdx + (threadsPerQuery - intraQueryThreadIdx - 1));
  }
}

GPU_INLINE __device__ void gpu_bpm_align_local_kernel(const gpu_bpm_align_qry_entry_t* const d_queries,  const gpu_bpm_align_device_qry_entry_t * const d_PEQs, const gpu_bpm_align_qry_info_t* const d_queryInfo,
                                                      const gpu_bpm_align_cand_entry_t * const d_candidates, const gpu_bpm_align_cand_info_t* const d_candidateInfo,
                                                      gpu_bpm_align_cigar_entry_t * const d_cigars, gpu_bpm_align_cigar_info_t* const d_cigarInfo,
                                                      const uint32_t idCandidate, const uint32_t intraQueryThreadIdx, const uint32_t threadsPerQuery)
{
  // Data characterization
  const uint32_t idQuery                             = d_candidateInfo[idCandidate].idQuery;
  const uint32_t idCigar                             = idCandidate;
  const uint32_t sizeQuery                           = d_queryInfo[idQuery].size;
  const uint32_t sizeCandidate                       = d_candidateInfo[idCandidate].size;
  const bool leftGapAlign                            = d_candidateInfo[idCandidate].leftGapAlign;
  // Data Buffers
  const uint64_t* const candidate   			           = (uint64_t*) (d_candidates + d_candidateInfo[idCandidate].posEntryBase);
  const uint64_t* const query                        = (uint64_t*) (d_queries + d_queryInfo[idQuery].posEntryBase);
  const gpu_bpm_align_device_qry_entry_t* const PEQs = d_PEQs + d_queryInfo[idQuery].posEntryPEQ + intraQueryThreadIdx;
  gpu_bpm_align_cigar_entry_t* const cigar           = d_cigars + d_cigarInfo[idCigar].offsetCigarStart + d_cigarInfo[idCigar].offsetCigarLength - 1;
  gpu_bpm_align_cigar_info_t* const cigarInfo        = d_cigarInfo + idCigar;

  // Local Memory (DP Matrix allocated in the CUDA stack)
  uint4 dpPV[GPU_BPM_ALIGN_MAX_SIZE_CANDIDATE];
  uint4 dpMV[GPU_BPM_ALIGN_MAX_SIZE_CANDIDATE];
  // DP matrix conversion for back-trace data-layout specialization
  const uint32_t* const dpPV4 = (uint32_t*) dpPV;
  const uint32_t* const dpMV4 = (uint32_t*) dpMV;

  //Return values for align DP matrix
  uint32_t minColumn = 0, minScore = sizeQuery;
  //Return values for align back-trace
  gpu_bpm_align_coord_t initCood = {0,0};
  uint32_t matchEffLenght = 0, cigarLenght = 0;

  gpu_bpm_align_dp_matrix(dpPV,  dpMV, PEQs, candidate, sizeQuery, sizeCandidate, intraQueryThreadIdx, threadsPerQuery,
                          &minColumn, &minScore);
  gpu_bpm_align_backtrace(dpPV4, dpMV4, query, candidate, cigar, leftGapAlign, minColumn, sizeQuery, intraQueryThreadIdx, threadsPerQuery,
                          &initCood, &matchEffLenght, &cigarLenght);

  // Return the cigar results
  cigarInfo->initCood       = initCood;
  cigarInfo->endCood.x      = minColumn;
  cigarInfo->endCood.y      = sizeQuery - 1;
  cigarInfo->matchEffLenght = matchEffLenght;
  cigarInfo->cigarStartPos  = sizeQuery - cigarLenght - 1;
  cigarInfo->cigarLenght    = cigarLenght;
}

__global__ void gpu_bpm_align_kernel(const gpu_bpm_align_qry_entry_t* const d_queries,  const gpu_bpm_align_device_qry_entry_t * const d_PEQs, const gpu_bpm_align_qry_info_t* const d_queryInfo,
                                     const gpu_bpm_align_cand_entry_t * const d_candidates, const gpu_bpm_align_cand_info_t* const d_candidateInfo, const uint32_t* const d_reorderBuffer,
                                     gpu_bpm_align_cigar_entry_t * const d_cigars, gpu_bpm_align_cigar_info_t* const d_cigarInfo, const uint32_t numCigars,
                                     const uint32_t* const d_initPosPerBucket, const uint32_t* const d_initWarpPerBucket, const bool updateScheduling)
{
  // Thread Identification
  const uint32_t globalThreadIdx = gpu_get_thread_idx();
  uint32_t idCandidate = 0, intraQueryThreadIdx = 0, threadsPerQuery = 0;
  // Rescheduling thread mapping and thread set distribution
  gpu_scheduler_scatter_work(globalThreadIdx, d_initWarpPerBucket, d_initPosPerBucket,
                             &idCandidate, &intraQueryThreadIdx, &threadsPerQuery);
  // Update the buffer input/output for the thread re-scheduling
  idCandidate  = (updateScheduling) ? d_reorderBuffer[idCandidate] : idCandidate;
  // Call to the device align BPM process for the active threads
  if (idCandidate < numCigars){
    // Update the buffer input/output for the thread re-scheduling
    gpu_bpm_align_local_kernel(d_queries, d_PEQs, d_queryInfo, d_candidates, d_candidateInfo, d_cigars, d_cigarInfo,
                               idCandidate, intraQueryThreadIdx, threadsPerQuery);
  }
}

extern "C"
gpu_error_t gpu_bpm_align_process_buffer(gpu_buffer_t *mBuff)
{
  // Internal buffer handles
  const gpu_bpm_align_queries_buffer_t* const     qry               = &mBuff->data.abpm.queries;
  const gpu_bpm_align_candidates_buffer_t* const  cand              = &mBuff->data.abpm.candidates;
  const gpu_scheduler_buffer_t* const             rebuff            = &mBuff->data.abpm.reorderBuffer;
  const gpu_bpm_align_cigars_buffer_t* const      cigar             = &mBuff->data.abpm.cigars;
  // Device properties
  const cudaStream_t                              idStream          =  mBuff->listStreams[mBuff->idStream];
  const uint32_t                                  idSupDev          =  mBuff->idSupportedDevice;
  const gpu_device_info_t* const                  device            =  mBuff->device[idSupDev];
  // Buffer size parameters information
  const uint32_t                                  numQueries        =  qry->numQueries;
  const uint32_t                                  numQueryBases     =  qry->totalQueriesBases;
  const uint32_t                                  numQueryPEQs      =  qry->totalQueriesPEQs;
  const uint32_t                                  numCandidates     =  cand->numCandidates;
  const uint32_t                                  numCandidateBases =  cand->numCandidatesBases;
  // Buffer size parameters for maximal threshold
  const uint32_t                                  maxQueries        =  mBuff->data.abpm.maxQueries;
  const uint32_t                                  maxQueryBases     =  mBuff->data.abpm.maxQueryBases;
  const uint32_t                                  maxQueryPEQs      =  mBuff->data.abpm.maxPEQEntries;
  const uint32_t                                  maxCandidates     =  mBuff->data.abpm.maxCandidates;
  const uint32_t                                  maxCandidateBases =  mBuff->data.abpm.maxCandidateBases;
  const uint32_t                                  maxCigars         =  mBuff->data.abpm.maxCigars;
  // Cigar results information
  const uint32_t                                  numCigarSched     = (mBuff->data.abpm.queryBinning) ? cigar->numReorderedCigars : cigar->numCigars;
  gpu_bpm_align_cigar_info_t*                     cigarsInfo        = (mBuff->data.abpm.queryBinning) ? cigar->d_reorderCigars : cigar->d_cigarsInfo;

  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = rebuff->numWarps * GPU_WARP_SIZE;
  gpu_device_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);
  // Sanity-check (checks buffer overflowing)
  if((numQueries > maxQueries) || (numCandidates > maxCandidates) || (numQueryPEQs > maxQueryPEQs) ||
     (numQueryBases > maxQueryBases) || (numCandidateBases > maxCandidateBases) || (numCigarSched > maxCigars))
    return(E_OVERFLOWING_BUFFER);
  // Launching the BPM align kernel on device
  gpu_bpm_align_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>(qry->d_queries, (gpu_bpm_align_device_qry_entry_t *) qry->d_peq, qry->d_qinfo,
                                                                        cand->d_candidates, cand->d_candidatesInfo, rebuff->d_reorderBuffer, cigar->d_cigars, cigarsInfo, numCigarSched,
                                                                        rebuff->d_initPosPerBucket, rebuff->d_initWarpPerBucket, mBuff->data.abpm.queryBinning);
  return(SUCCESS);
}

#endif /* GPU_BPM_ALIGN_CU_ */
