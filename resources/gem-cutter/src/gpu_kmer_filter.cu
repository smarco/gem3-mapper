/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_KMER_FILTER_CU_
#define GPU_KMER_FILTER_CU_

#include "../include/gpu_kmer_core.h"
#include "../include/gpu_text_core.h"


// Compile Pattern
GPU_INLINE __device__ void gpu_kmer_filter_compile_pattern(uint16_t* const kmerCountQuery, const uint64_t* const query, const uint32_t queryLength)
{
  // Count kmers in query
  const uint32_t baseQueryLength = GPU_KMER_FILTER_BASE_QUERY_LENGTH;
  uint32_t pos = 0, kmerIdx = 0;
  ulong2   infoQuery = GPU_TEXT_INIT;
  // Compose the first k-mer
  for (pos = 0; pos < (GPU_KMER_FILTER_COUNTING_LENGTH - 1); ++pos) {
    const uint8_t encBase = gpu_text_lookup(query, pos, &infoQuery, baseQueryLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdx, encBase); // Update kmer-index
  }
  // Compile all k-mers
  for (pos = (GPU_KMER_FILTER_COUNTING_LENGTH - 1); pos < queryLength; ++pos) {
    const uint8_t encBase = gpu_text_lookup(query, pos, &infoQuery, baseQueryLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdx, encBase); // Update kmer-index
    kmerCountQuery[kmerIdx]++;                            // Increment kmer-count
  }
}

// Filter text region
GPU_INLINE __device__ void gpu_kmer_filter_candidate_cutoff(const uint16_t* const kmerCountQuery, uint16_t* const kmerCountCandidate,
                                                            const uint32_t queryLength, const uint32_t candidateLength,
                                                            const uint64_t* const candidates, const uint64_t offsetCandPos,
                                                            const uint32_t maxErrorRatio, uint32_t* const filterDistance)
{
  // Prepare filter
  const uint32_t baseCandidateLength = GPU_KMER_FILTER_BASE_CANDIDATE_LENGTH;
  const uint32_t maxError = queryLength * (maxErrorRatio / 100);
  const uint32_t kmersRequired = queryLength - (GPU_KMER_FILTER_COUNTING_LENGTH - 1) - GPU_KMER_FILTER_COUNTING_LENGTH * maxError;
  const uint32_t totalKmersCandidate = GPU_MAX(candidateLength, queryLength);
  const uint32_t initChunk = GPU_MIN(candidateLength, queryLength);
  uint32_t kmersLeft = totalKmersCandidate, kmersInCandidate = 0;
  uint32_t beginPos = 0, kmerIdxBegin = 0, endPos = 0, kmerIdxEnd = 0, cutOffResult = GPU_UINT32_ZEROS;
  ulong2   infoCandidate = GPU_TEXT_INIT;
  // Initial window fill (Composing the first k-mer)
  while (endPos < (GPU_KMER_FILTER_COUNTING_LENGTH - 1)){
    const uint8_t encChar = gpu_text_lookup(candidates, offsetCandPos + endPos, &infoCandidate, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxEnd, encChar); // Update kmer-index
    ++endPos;
  }
  // Initial window fill
  while ((endPos < initChunk) && !cutOffResult) {
      const uint8_t encBaseEnd = gpu_text_lookup(candidates, offsetCandPos + endPos, &infoCandidate, baseCandidateLength);
      GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxEnd, encBaseEnd); // Update kmer-index
      const uint16_t countQuery = kmerCountQuery[kmerIdxEnd];
      if (countQuery > 0 && kmerCountCandidate[kmerIdxEnd]++ < countQuery) ++kmersInCandidate;
      // Check filter condition
      if (kmersInCandidate >= kmersRequired)
        cutOffResult = GPU_ALIGN_DISTANCE_ZERO; // Don't filter
      if ((kmersInCandidate < kmersRequired) && (kmersRequired - kmersInCandidate) > kmersLeft)
        cutOffResult = GPU_ALIGN_DISTANCE_INF;  // Filter out
      ++endPos, --kmersLeft;
  }
  // Sliding window count (Composing the first k-mer)
  while (beginPos < (GPU_KMER_FILTER_COUNTING_LENGTH - 1) && !cutOffResult) {
    const uint8_t encBaseBegin = gpu_text_lookup(candidates, offsetCandPos + beginPos, &infoCandidate, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxBegin, encBaseBegin);
    ++beginPos;
  }
  // Sliding window count (End processing)
  while ((endPos < candidateLength) && !cutOffResult) {
    // Begin (Decrement kmer-count)
    const uint8_t encBaseBegin = gpu_text_lookup(candidates, offsetCandPos + beginPos, &infoCandidate, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxBegin, encBaseBegin);
    const uint16_t countQueryBegin = kmerCountQuery[kmerIdxBegin];
    if (countQueryBegin > 0 && kmerCountCandidate[kmerIdxBegin]-- <= countQueryBegin) --kmersInCandidate;
    // End (Increment kmer-count)
    const uint8_t encBaseEnd = gpu_text_lookup(candidates, offsetCandPos + endPos, &infoCandidate, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxEnd, encBaseEnd);
    const uint16_t countQueryEnd = kmerCountQuery[kmerIdxEnd];
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

// Filter text region
GPU_INLINE __device__ void gpu_kmer_filter_candidate_sliding_window(const uint16_t* const kmerCountQuery, uint16_t* const kmerCountCandidate,
                                                                    const uint32_t queryLength, const uint32_t candidateLength,
                                                                    const uint64_t* const candidates, const uint64_t offsetCandPos,
                                                                    uint32_t* const filterDistance)
{
  // Prepare filter
  const uint32_t baseCandidateLength = GPU_KMER_FILTER_BASE_CANDIDATE_LENGTH;
  const uint32_t initChunk = GPU_MIN(candidateLength, queryLength);
  const uint32_t maxKmers  = GPU_MIN(candidateLength, queryLength) - (GPU_KMER_FILTER_COUNTING_LENGTH - 1);
  uint32_t beginPos = 0, kmerIdxBegin = 0, endPos = 0, kmerIdxEnd = 0;
  ulong2   infoCandidateBegin = GPU_TEXT_INIT, infoCandidateEnd = GPU_TEXT_INIT;
  int32_t  kmersInCandidate = 0, kmerDistance = 0;
  // Initial window fill (Composing the first k-mer)
  while (endPos < (GPU_KMER_FILTER_COUNTING_LENGTH - 1)){
    const uint8_t encChar = gpu_text_lookup(candidates, offsetCandPos + endPos, &infoCandidateEnd, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxEnd, encChar); // Update kmer-index
    ++endPos;
  }
  beginPos = endPos; kmerIdxBegin = kmerIdxEnd;
  // Initial window fill
  while (endPos < initChunk) {
      const uint8_t encBaseEnd = gpu_text_lookup(candidates, offsetCandPos + endPos, &infoCandidateEnd, baseCandidateLength);
      GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxEnd, encBaseEnd); // Update kmer-index
      const uint16_t countQuery = kmerCountQuery[kmerIdxEnd];
      if ((countQuery > 0) && (kmerCountCandidate[kmerIdxEnd]++ <= countQuery)) ++kmersInCandidate;
      ++endPos;
  }
  // Sliding window count (End processing)
  kmerDistance = kmersInCandidate;
  while (endPos < candidateLength) {
    // Begin (Decrement kmer-count)
    const uint8_t encBaseBegin = gpu_text_lookup(candidates, offsetCandPos + beginPos, &infoCandidateBegin, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxBegin, encBaseBegin);
    const uint16_t countQueryBegin = kmerCountQuery[kmerIdxBegin];
    if ((countQueryBegin > 0) && ((kmerCountCandidate[kmerIdxBegin]--) <= countQueryBegin)) --kmersInCandidate;
    // End (Increment kmer-count)
    const uint8_t encBaseEnd = gpu_text_lookup(candidates, offsetCandPos + endPos, &infoCandidateEnd, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxEnd, encBaseEnd);
    const uint16_t countQueryEnd = kmerCountQuery[kmerIdxEnd];
    if ((countQueryEnd > 0) && (kmerCountCandidate[kmerIdxEnd]++ <= countQueryEnd)) ++kmersInCandidate;
    kmerDistance = GPU_MAX(kmerDistance, kmersInCandidate);
    // Advance sliding window
    ++endPos;
    ++beginPos;
  }
  // kmer filtering distance
  (* filterDistance) = (maxKmers - kmerDistance) / GPU_KMER_FILTER_COUNTING_LENGTH;
}


// Filter text region
GPU_INLINE __device__ void gpu_kmer_filter_candidate(const uint16_t* const kmerCountQuery, uint16_t* const kmerCountCandidate,
                                                     const uint32_t queryLength, const uint32_t candidateLength,
                                                     const uint64_t* const candidate, const uint64_t offsetCandPos,
                                                     uint32_t* const filterDistance)
{
  // Prepare filter
  const uint32_t baseCandidateLength = GPU_KMER_FILTER_BASE_CANDIDATE_LENGTH;
  const uint32_t maxKmers = GPU_MIN(candidateLength, queryLength) - (GPU_KMER_FILTER_COUNTING_LENGTH - 1);
  uint32_t endPos = 0, kmerIdxEnd = 0;
  ulong2   infoCandidateEnd = GPU_TEXT_INIT;
  int32_t  kmersInCandidate = 0, kmerDistance = 0;
  // Initial window fill (Composing the first k-mer)
  while (endPos < (GPU_KMER_FILTER_COUNTING_LENGTH - 1)){
    const uint8_t encChar = gpu_text_lookup(candidate, offsetCandPos + endPos, &infoCandidateEnd, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxEnd, encChar); // Update kmer-index
    endPos++;
  }
  // Sliding window count (End processing)
  while (endPos < candidateLength) {
    const uint8_t encBaseEnd = gpu_text_lookup(candidate, offsetCandPos + endPos, &infoCandidateEnd, baseCandidateLength);
    GPU_KMER_FILTER_COUNTING_ADD_INDEX(kmerIdxEnd, encBaseEnd);
    const uint16_t countQueryEnd     = kmerCountQuery[kmerIdxEnd];
    const uint16_t countCandidateEnd = kmerCountCandidate[kmerIdxEnd];
    if ((countQueryEnd > 0) && (countCandidateEnd < countQueryEnd)){
      kmersInCandidate++;
      kmerDistance = GPU_MAX(kmerDistance, kmersInCandidate);
    }
    // Advance sliding window and increment tiling
    kmerCountCandidate[kmerIdxEnd]++;
    endPos++;
  }
  // kmer filtering distance
  (* filterDistance) = (maxKmers - kmerDistance) / GPU_KMER_FILTER_COUNTING_LENGTH;
}

//k-mer filter device kernel
__global__ void gpu_kmer_filter_kernel(const gpu_kmer_filter_qry_entry_t* const d_queries, const gpu_kmer_filter_qry_info_t* const d_qinfo,
                                       const uint64_t* const reference, const gpu_kmer_filter_cand_info_t* const d_candidates, const uint64_t sizeRef,
                                       gpu_kmer_filter_alg_entry_t* const d_filterResult, const uint32_t maxError,  const uint32_t numAlignments)
{
  const uint32_t idAlignment = gpu_get_thread_idx();
  if(idAlignment < numAlignments){
    const uint64_t positionRef     = d_candidates[idAlignment].position;
    const uint32_t sizeCandidate   = d_candidates[idAlignment].size;
    const uint32_t sizeQuery       = d_qinfo[d_candidates[idAlignment].query].query_size;
    const uint64_t* const query    = (uint64_t*) (d_queries + d_qinfo[d_candidates[idAlignment].query].init_offset);
    // Filter initialization
    uint32_t filterDistance = 0;
    uint16_t kmerCountCandidate[GPU_KMER_FILTER_COUNTING_NUM_KMERS] = {0};
    uint16_t kmerCountQuery[GPU_KMER_FILTER_COUNTING_NUM_KMERS]     = {0};
    // Compile all k-mers from pattern
    gpu_kmer_filter_compile_pattern(kmerCountQuery, query, sizeQuery);
    // Compile all k-mers from text
    gpu_kmer_filter_candidate(kmerCountQuery, kmerCountCandidate, sizeQuery, sizeCandidate, reference, positionRef, &filterDistance);
    d_filterResult[idAlignment] = filterDistance;
  }
}


extern "C"
gpu_error_t gpu_kmer_filter_process_buffer(gpu_buffer_t *mBuff)
{
  const gpu_reference_buffer_t* const              ref           =  mBuff->reference;
  const gpu_kmer_filter_queries_buffer_t* const    qry           = &mBuff->data.fkmer.queries;
  const gpu_kmer_filter_candidates_buffer_t* const cand          = &mBuff->data.fkmer.candidates;
  const gpu_kmer_filter_alignments_buffer_t* const res           = &mBuff->data.fkmer.alignments;
  const cudaStream_t                               idStream      =  mBuff->listStreams[mBuff->idStream];
  const uint32_t                                   idSupDev      =  mBuff->idSupportedDevice;
  const gpu_device_info_t* const                   device        =  mBuff->device[idSupDev];
  const uint32_t                                   maxBases      =  mBuff->data.fkmer.maxBases;
  const uint32_t                                   maxQueries    =  mBuff->data.fkmer.maxQueries;
  const uint32_t                                   maxCandidates =  mBuff->data.fkmer.maxCandidates;
  const uint32_t                                   maxAlignments =  mBuff->data.fkmer.maxAlignments;
  const uint32_t                                   maxError      =  mBuff->data.fkmer.maxError;

  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = res->numAlignments;
  gpu_device_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);
  // Sanity-check (checks buffer overflowing)
  if((qry->numBases > maxBases) || (qry->numQueries > maxQueries) || (res->numAlignments > maxCandidates) || (res->numAlignments > maxAlignments))
    return(E_OVERFLOWING_BUFFER);

  gpu_kmer_filter_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>(qry->d_queries, qry->d_queryInfo, ref->d_reference_plain[idSupDev],
                                                                          cand->d_candidates, ref->size, res->d_alignments,
                                                                          maxError, res->numAlignments);
  return(SUCCESS);
}

#endif /* GPU_KMER_FILTER_CU_ */
