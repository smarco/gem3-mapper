/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_KMER_FILTER_CU_
#define GPU_KMER_FILTER_CU_

#include "../include/gpu_kmer_core.h"

// Compile Pattern
GPU_INLINE __device__ void gpu_kmer_filter_compile_pattern(uint16_t* const kmerCountQuery, const uint8_t* const query, const uint32_t queryLength)
{
  // Count kmers in query
  uint32_t pos = 0, kmerIdx = 0;
  // Compose the first k-mer
  for (pos = 0; pos < (GPU_KMER_COUNTING_LENGTH - 1); ++pos) {
    const uint8_t encBase = query[pos];
    GPU_KMER_COUNTING_ADD_INDEX(kmerIdx, encBase); // Update kmer-index
  }
  // Compile all k-mers
  for (pos = (GPU_KMER_COUNTING_LENGTH - 1); pos < queryLength; ++pos) {
    const uint8_t encBase = query[pos];
    GPU_KMER_COUNTING_ADD_INDEX(kmerIdx, encBase); // Update kmer-index
    ++(kmerCountQuery[kmerIdx]);                   // Increment kmer-count
  }
}

//      lastCandidateEntry = localCandidate[idEntry];
// for(idColumn = 0; idColumn < sizeCandidate; idColumn++){

//  if((idColumn % GPU_REFERENCE_CHARS_PER_ENTRY) == 0){
//    idEntry++;
//    currentCandidateEntry = localCandidate[idEntry];
//    candidate = funnelshift_left_64(currentCandidateEntry, lastCandidateEntry, candidateAlignment);
//    lastCandidateEntry = currentCandidateEntry;
//  }


GPU_INLINE __device__ void gpu_kmer_query_lookup(const uint64_t* const query, const uint32_t idBase, char* const globalBase,
                                                 uint64_t* const globalInfoBase, uint32_t* const globalInfoQuery)
{
  //Query meta information to reduce the main memory requests
  uint64_t infoBase  = (* globalInfoBase);  // Packet 8 bases (cached)
  uint32_t infoQuery = (* globalInfoQuery); // Position cached
  //Address the cached block request (hash function)
  const uint32_t idEntryQuery = idBase / GPU_KMER_BASES_PER_QUERY_ENTRY;
  const uint32_t idIntraQuery = idBase % GPU_KMER_BASES_PER_QUERY_ENTRY;
  uint32_t shiftBits, base;
  //Make the query request if data was not cached
  if((infoQuery != idEntryQuery) || (infoQuery == GPU_UINT32_ONES)){
    infoBase  = LDG(query + idEntryQuery);
    infoQuery = idEntryQuery;
  }
  //logical base extraction from cached request
  shiftBits = idIntraQuery * GPU_KMER_BASE_QUERY_LENGTH;
  base      = (infoBase >> shiftBits) & GPU_KMER_BASE_QUERY_MASK; //0xFF;
  //return the requested data
  (* globalBase)      = base;      //requested base
  (* globalInfoBase)  = infoBase;  //contain query cached data
  (* globalInfoQuery) = infoQuery; //contain the query position cached data
}

GPU_INLINE __device__ void gpu_kmer_candidate_lookup(const uint64_t* const reference, const uint32_t refPosition, char* const globalBase,
                                                     uint64_t* const globalInfoBase, uint32_t* const globalInfoCandidate)
{
  //Query meta information to reduce the main memory requests
  uint64_t infoBase      = (* globalInfoBase);      // Packet 8 bases (cached)
  uint32_t infoCandidate = (* globalInfoCandidate); // Position cached
  //Address the cached block request (hash function)
  const uint32_t idEntryCandidate = refPosition / GPU_KMER_BASES_PER_CANDIDATE_ENTRY;
  const uint32_t idIntraCandidate = refPosition % GPU_KMER_BASES_PER_CANDIDATE_ENTRY;
  uint32_t shiftBits, base;
  //Make the query request if data was not cached
  if((infoCandidate != idEntryCandidate) || (infoCandidate == GPU_UINT32_ONES)){
    infoBase      = reference[idEntryReference];
    infoCandidate = idEntryCandidate;
  }
  //logical base extraction from cached request
  shiftBits = idIntraCandidate * GPU_KMER_BASE_CANDIDATE_LENGTH;
  base      = (infoBase >> shiftBits) & GPU_KMER_BASE_CANDIDATE_MASK; // 0x03;
  //return the requested data
  (* globalBase)          = base;      //requested base
  (* globalInfoBase)      = infoBase;  //contain query cached data
  (* globalInfoCandidate) = infoQuery; //contain the query position cached data
}

// Filter text region
GPU_INLINE __device__ void gpu_kmer_filter_candidate(const uint16_t* const kmerCountQuery, const uint8_t * const query, const uint32_t queryLength,
                                                     uint16_t* const kmerCountCandidate, const uint64_t* const candidate, const uint32_t candidateLength,
                                                     const uint32_t maxError, uint32_t* const filterDistance)
{
  // Prepare filter
  const uint32_t kmersRequired = queryLength - (GPU_KMER_COUNTING_LENGTH - 1) - GPU_KMER_COUNTING_LENGTH * maxError;
  const uint32_t totalKmersCandidate = GPU_MAX(candidateLength, queryLength);
  const uint32_t initChunk = GPU_MIN(candidateLength, queryLength);
  uint32_t kmersLeft = totalKmersCandidate, kmersInCandidate = 0;
  uint32_t beginPos = 0, kmerIdxBegin = 0;
  uint32_t endPos = 0, kmerIdxEnd = 0;
  // Initial window fill (Composing the first k-mer)
  while (endPos < (GPU_KMER_COUNTING_LENGTH - 1)){
    const uint8_t encChar = candidate[endPos];
    GPU_KMER_COUNTING_ADD_INDEX(kmerIdxEnd, encChar); // Update kmer-index
    ++endPos;
  }
  // Initial window fill
  while ((endPos < initChunk) && !cutOffResult) {
    const uint8_t encBaseEnd = candidate[endPos];
      GPU_KMER_COUNTING_ADD_INDEX(kmerIdxEnd, encBaseEnd); // Update kmer-index
      const uint16_t countPattern = kmerCountQuery[kmerIdxEnd];
      if (countQuery > 0 && (kmerCountCandidate[kmerIdxEnd])++ < countQuery) ++kmersInCandidate;
      // Check filter condition
      if (kmersInCandidate >= kmersRequired)
        cutOffResult = GPU_ALIGN_DISTANCE_ZERO; // Don't filter
      if ((kmersInCandidate < kmersRequired) && (kmersRequired - kmersInCandidate) > kmersLeft)
        cutOffResult = GPU_ALIGN_DISTANCE_INF;  // Filter out
      ++endPos, --kmersLeft;
  }
  // Sliding window count (Composing the first k-mer)
  while (beginPos < (GPU_KMER_COUNTING_LENGTH - 1) && !cutOffResult) {
    const uint8_t encBaseBegin = candidate[beginPos];
    GPU_KMER_COUNTING_ADD_INDEX(kmerIdxBegin, encBaseBegin);
    ++beginPos;
  }
  // Sliding window count (End processing)
  while ((endPos < candidateLength) && !cutOffResult) {
    // Begin (Decrement kmer-count)
    const uint8_t encBaseBegin = candidate[beginPos];
    GPU_KMER_COUNTING_ADD_INDEX(kmerIdxBegin, encBaseBegin);
    const uint16_t countQueryBegin = kmerCountQuery[kmerIdxBegin];
    if (countQueryBegin > 0 && kmerCountCandidate[kmerIdxBegin]-- <= countQueryBegin) --kmersInCandidate;
    // End (Increment kmer-count)
    const uint8_t encBaseEnd = candidate[endPos];
    GPU_KMER_COUNTING_ADD_INDEX(kmerIdxEnd, encBaseEnd);
    const uint16_t countPatternEnd = kmerCountQuery[kmerIdxEnd];
    if (countQueryEnd > 0 && kmerCountCandidate[kmerIdxEnd]++ < countQueryEnd) ++kmersInCandidate;
    // Check filter condition
    if (kmersInCandidate >= kmersRequired)
      cutOffResult = GPU_ALIGN_DISTANCE_ZERO; // Don't filter
    if ((kmersInCandidate < kmersRequired) && (kmersRequired - kmersInCandidate) > kmersLeft)
      cutOffResult = GPU_ALIGN_DISTANCE_INF;  // Filter out
    ++endPos; ++beginPos; --kmersLeft;
  }
  // Not passing the filter
  if (cutOffResult) (* filterDistance) = cutOffResult;
    else (* filterDistance) = GPU_ALIGN_DISTANCE_INF; // Filter out
}

//k-mer filter device kernel
__global__ void gpu_kmer_filter_kernel(const gpu_kmer_qry_entry_t* const d_queries, const gpu_kmer_qry_info_t* const d_qinfo,
                                      const uint64_t* const d_reference, const gpu_kmer_cand_info_t* const d_candidates, const uint64_t sizeRef,
                                      gpu_kmer_alg_entry_t* const d_filterResult, const uint32_t maxError,  const uint32_t numAlignments)
{
  const uint32_t idAlignment = gpu_get_thread_idx();
  if(idAlignment < numAlignments){
    const uint64_t positionRef     = d_candidates[idAlignment].position;
    const uint32_t sizeCandidate   = d_candidates[idAlignment].size;
    const uint32_t sizeQuery       = d_qinfo[d_candidates[idAlignment].query].query_size;
    const uint8_t* const  query    = d_queries + d_qinfo[d_candidates[idAlignment].query].init_offset;
    // Filter initialization
    uint32_t filterDistance = 0;
    uint16_t kmerCountCandidate[GPU_KMER_COUNTING_NUM_KMERS] = {0};
    uint16_t kmerCountQuery[GPU_KMER_COUNTING_NUM_KMERS]     = {0};
    // Compile all k-mers from pattern
    gpu_kmer_filter_compile_pattern(kmerCountQuery, query, queryLength);
    // Compile all k-mers from text
    gpu_kmer_filter_candidate(kmerCountQuery, query, sizeQuery, kmerCountCandidate, candidate, sizeCandidate, maxError, &filterDistance);
    if(filterDistance == GPU_ALIGN_DISTANCE_ZERO) filterDistance = 0;
    d_filterResult[idAlignment] = filterDistance;
  }
}


extern "C"
gpu_error_t gpu_kmer_process_buffer(gpu_buffer_t *mBuff)
{
  const gpu_reference_buffer_t* const       ref           =  mBuff->reference;
  const gpu_kmer_queries_buffer_t* const    qry           = &mBuff->data.kmer.queries;
  const gpu_kmer_candidates_buffer_t* const cand          = &mBuff->data.kmer.candidates;
  const gpu_kmer_alignments_buffer_t* const res           = &mBuff->data.kmer.alignments;
  const cudaStream_t                        idStream      =  mBuff->idStream;
  const uint32_t                            idSupDev      =  mBuff->idSupportedDevice;
  const gpu_device_info_t* const            device        =  mBuff->device[idSupDev];
  const uint32_t                            numAlignments =  res->numAlignments;
  const uint32_t                            maxBases      =  mBuff->data.kmer.maxBases;
  const uint32_t                            maxQueries    =  mBuff->data.kmer.maxQueries;
  const uint32_t                            maxCandidates =  mBuff->data.kmer.maxCandidates;
  const uint32_t                            maxAlignments =  mBuff->data.kmer.maxAlignments;
  const uint32_t                            maxError      =  mBuff->data.kmer.maxError;
  const uint32_t                            numResults    =  res->numAlignments;
  gpu_bpm_alg_entry_t*                      d_results     =  res->d_alignments;

  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = numResults;
  gpu_device_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);
  // Sanity-check (checks buffer overflowing)
  if((numBases > maxBases) || (numQueries > maxQueries) || (numAlignments > maxCandidates) || (numAlignments > maxAlignments))
    return(E_OVERFLOWING_BUFFER);

  gpu_kmer_filter_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>(qry->d_queries, qry->d_queryInfo, ref->d_reference[idSupDev],
                                                                          cand->d_candidates, ref->size, res->d_alignments,
                                                                          maxError, res->numAlignments);
  return(SUCCESS);
}

#endif /* GPU_KMER_FILTER_CU_ */
