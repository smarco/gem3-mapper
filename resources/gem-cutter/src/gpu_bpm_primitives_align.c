/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_BPM_PRIMITIVES_ALIGN_C_
#define GPU_BPM_PRIMITIVES_ALIGN_C_

#include "../include/gpu_bpm_primitives.h"

/************************************************************
Functions to get the GPU BPM buffer sizes
************************************************************/

uint32_t gpu_bpm_align_buffer_get_max_peq_entries_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.maxPEQEntries);
}

uint32_t gpu_bpm_align_buffer_get_max_candidates_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.maxCandidates);
}

uint32_t gpu_bpm_align_buffer_get_max_candidate_size_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.maxCandidateSize);
}

uint32_t gpu_bpm_align_buffer_get_max_queries_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.maxQueries);
}

uint32_t gpu_bpm_align_buffer_get_max_query_bases_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.maxQueryBases);
}

uint32_t gpu_buffer_bpm_align_get_max_cigar_entries_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.maxCigarEntries);
}

/************************************************************
Functions to get the GPU BPM buffers
************************************************************/

gpu_bpm_align_qry_entry_t* gpu_bpm_align_buffer_get_queries_(const void* const bpmBuffer){
  const gpu_buffer_t * const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.queries.h_queries);
}

gpu_bpm_align_peq_entry_t* gpu_bpm_align_buffer_get_peq_entries_(const void* const bpmBuffer){
  const gpu_buffer_t * const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.queries.h_peq);
}

gpu_bpm_align_qry_info_t* gpu_bpm_align_buffer_get_queries_info_(const void* const bpmBuffer){
  const gpu_buffer_t * const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.queries.h_qinfo);
}

gpu_bpm_align_cand_info_t* gpu_bpm_align_buffer_get_candidates_info_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.candidates.h_candidatesInfo);
}

gpu_bpm_align_cigar_entry_t* gpu_bpm_align_buffer_get_cigars_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.cigars.h_cigars);
}

gpu_bpm_align_cigar_info_t* gpu_bpm_align_buffer_get_cigars_info_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.abpm.cigars.h_cigarsInfo);
}



/************************************************************
Functions to init all the BPM resources
************************************************************/

float gpu_bpm_align_size_per_candidate(const uint32_t averageQuerySize, const uint32_t candidatesPerQuery)
{
  const size_t averageNumPEQEntries = GPU_DIV_CEIL(averageQuerySize, GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH);
  const size_t averageCigarSize     = averageQuerySize + 1;
  const size_t bytesPerQueryPEQ     = averageNumPEQEntries * sizeof(gpu_bpm_align_peq_entry_t);
  const size_t bytesPerQueryRaw     = GPU_ROUND_TO(averageQuerySize * sizeof(gpu_bpm_align_qry_entry_t), 8);
  const size_t bytesPerQuery        = bytesPerQueryRaw + bytesPerQueryPEQ + sizeof(gpu_bpm_align_qry_info_t);
  const size_t bytesCigar           = (averageCigarSize * sizeof(gpu_bpm_align_cigar_entry_t)) + sizeof(gpu_bpm_align_cigar_info_t);
  const size_t bytesBinningProcess  = sizeof(uint32_t);
  // Calculate the necessary bytes for each BPM Align operation
  const size_t sizePerCandidate     = (bytesPerQuery / (float)candidatesPerQuery) + bytesCigar + bytesBinningProcess;
  return(sizePerCandidate);
}

uint32_t gpu_bpm_align_candidates_for_binning_padding()
{
  uint32_t idBucket;
  uint32_t bucketPaddingCandidates = 0;
  /* Worst number of dummy candidates added for padding the binning */
  for(idBucket = 1; idBucket < GPU_BPM_ALIGN_NUM_BUCKETS_FOR_BINNING-1; ++idBucket)
    bucketPaddingCandidates += (GPU_WARP_SIZE / idBucket);
  /* Increased bytes per buffer taking account the padding*/
  return(bucketPaddingCandidates);
}

void gpu_bpm_align_reallocate_host_buffer_layout(gpu_buffer_t* mBuff)
{
  void* rawAlloc = mBuff->h_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.abpm.queries.h_queries = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.queries.h_queries + mBuff->data.abpm.maxQueryBases);
  mBuff->data.abpm.queries.h_peq = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.queries.h_peq + mBuff->data.abpm.maxPEQEntries);
  mBuff->data.abpm.queries.h_qinfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.queries.h_qinfo + mBuff->data.abpm.maxQueries);
  //Allocate space for the candidate structures
  mBuff->data.abpm.candidates.h_candidatesInfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.candidates.h_candidatesInfo + mBuff->data.abpm.maxCandidates);
  //Allocate space for the reorder buffer structures
  mBuff->data.abpm.reorderBuffer.threadMapScheduler.h_reorderBuffer = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.reorderBuffer.threadMapScheduler.h_reorderBuffer + mBuff->data.abpm.maxReorderBuffer);
  mBuff->data.abpm.reorderBuffer.h_initPosPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.reorderBuffer.h_initPosPerBucket + mBuff->data.abpm.maxBuckets);
  mBuff->data.abpm.reorderBuffer.h_initWarpPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.reorderBuffer.h_initWarpPerBucket + mBuff->data.abpm.maxBuckets);
  mBuff->data.abpm.reorderBuffer.h_endPosPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.reorderBuffer.h_endPosPerBucket + mBuff->data.abpm.maxBuckets);
  //Adjust the host buffer layout (output)
  mBuff->data.abpm.cigars.h_cigarsInfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.cigars.h_cigarsInfo + mBuff->data.abpm.maxCigars);
  mBuff->data.abpm.cigars.h_cigars = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.cigars.h_cigars + mBuff->data.abpm.maxCigarEntries);
}

void gpu_bpm_align_reallocate_device_buffer_layout(gpu_buffer_t* mBuff)
{
  void* rawAlloc = mBuff->d_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.abpm.queries.d_queries = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.queries.d_queries + mBuff->data.abpm.maxQueryBases);
  mBuff->data.abpm.queries.d_peq = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.queries.d_peq + mBuff->data.abpm.maxPEQEntries);
  mBuff->data.abpm.queries.d_qinfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.queries.d_qinfo + mBuff->data.abpm.maxQueries);
  //Allocate space for the candidate structures
  mBuff->data.abpm.candidates.d_candidatesInfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.candidates.d_candidatesInfo + mBuff->data.abpm.maxCandidates);
  //Allocate space for the reorder buffer structures
  mBuff->data.abpm.reorderBuffer.threadMapScheduler.d_reorderBuffer = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.reorderBuffer.threadMapScheduler.d_reorderBuffer + mBuff->data.abpm.maxReorderBuffer);
  mBuff->data.abpm.reorderBuffer.d_initPosPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.reorderBuffer.d_initPosPerBucket + mBuff->data.abpm.maxBuckets);
  mBuff->data.abpm.reorderBuffer.d_initWarpPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.reorderBuffer.d_initWarpPerBucket + mBuff->data.abpm.maxBuckets);
  mBuff->data.abpm.reorderBuffer.d_endPosPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.reorderBuffer.d_endPosPerBucket + mBuff->data.abpm.maxBuckets);
  //Adjust the device buffer layout (output)
  mBuff->data.abpm.cigars.d_cigarsInfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.cigars.d_cigarsInfo + mBuff->data.abpm.maxCigars);
  mBuff->data.abpm.cigars.d_cigars = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.abpm.cigars.d_cigars + mBuff->data.abpm.maxCigarEntries);
}


void gpu_bpm_align_init_buffer_(void* const bpmBuffer, const uint32_t averageQuerySize, const uint32_t candidatesPerQuery)
{
  gpu_buffer_t* const mBuff                   = (gpu_buffer_t *) bpmBuffer;
  const double        sizeBuff                = mBuff->sizeBuffer * 0.95;
  const uint32_t      averageNumPEQEntries    = GPU_DIV_CEIL(averageQuerySize, GPU_BPM_ALIGN_PEQ_ENTRY_LENGTH);
  const uint32_t      averageCigarSize        = averageQuerySize + 1;
  const uint32_t      numInputs               = (uint32_t)(sizeBuff / gpu_bpm_align_size_per_candidate(averageQuerySize, candidatesPerQuery));
  const uint32_t      bucketPaddingCandidates = gpu_bpm_align_candidates_for_binning_padding();
  const uint32_t      maxCandidates           = numInputs;
  //set the type of the buffer
  mBuff->typeBuffer = GPU_BPM_ALIGN;
  //Set real size of the input
  mBuff->data.abpm.maxCandidates     = maxCandidates;
  mBuff->data.abpm.maxCigars         = maxCandidates;
  mBuff->data.abpm.maxQueries        = (maxCandidates / candidatesPerQuery);
  mBuff->data.abpm.maxCandidateSize  = GPU_BPM_ALIGN_MAX_SIZE_CANDIDATE;
  // Set internal data buffers sizes
  mBuff->data.abpm.maxPEQEntries     = mBuff->data.abpm.maxQueries * averageNumPEQEntries;
  mBuff->data.abpm.maxCigarEntries   = mBuff->data.abpm.maxCigars  * averageCigarSize;
  mBuff->data.abpm.maxQueryBases     = mBuff->data.abpm.maxQueries * averageQuerySize;
  // Set the reorder buffer size
  mBuff->data.abpm.maxReorderBuffer  = mBuff->data.abpm.maxCandidates + bucketPaddingCandidates;
  mBuff->data.abpm.maxBuckets        = GPU_BPM_ALIGN_NUM_BUCKETS_FOR_BINNING;
  mBuff->data.abpm.queryBinSize      = 0;
  mBuff->data.abpm.queryBinning      = true;
  // Set the corresponding buffer layout
  gpu_bpm_align_reallocate_host_buffer_layout(mBuff);
  gpu_bpm_align_reallocate_device_buffer_layout(mBuff);
}

void gpu_bpm_align_init_and_realloc_buffer_(void *bpmBuffer, const uint32_t totalPEQEntries, const uint32_t totalQueryBases,
                                            const uint32_t totalQueries, const uint32_t totalCandidates)
{
  // Buffer re-initialization
  gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  const uint32_t averageQuerySize       = GPU_DIV_CEIL(totalQueryBases, totalQueries);
  const uint32_t candidatesPerQuery     = GPU_DIV_CEIL(totalCandidates, totalQueries);
  // Re-map the buffer layout with new information trying to fit better
  gpu_bpm_align_init_buffer_(bpmBuffer, averageQuerySize, candidatesPerQuery);
  // Checking if we need to reallocate a bigger buffer
  if( (totalPEQEntries     > gpu_bpm_align_buffer_get_max_peq_entries_(bpmBuffer))     ||
      (totalQueryBases     > gpu_bpm_align_buffer_get_max_query_bases_(bpmBuffer))     ||
      (totalCandidates     > gpu_bpm_align_buffer_get_max_candidates_(bpmBuffer))      ||
      (totalQueries        > gpu_bpm_align_buffer_get_max_queries_(bpmBuffer))){
    // Resize the GPU buffer to fit the required input
    const uint32_t  idSupDevice             = mBuff->idSupportedDevice;
    const float     resizeFactor            = 2.0;
    const size_t    bytesPerBPMBuffer       = totalCandidates * gpu_bpm_align_size_per_candidate(averageQuerySize, candidatesPerQuery);
    //printf("RESIZE[BPM_ALIGN] %d %d \n",  mBuff->sizeBuffer, bytesPerBPMBuffer * resizeFactor);
    //Recalculate the minimum buffer size
    mBuff->sizeBuffer = bytesPerBPMBuffer * resizeFactor;
    //FREE HOST AND DEVICE BUFFER
    GPU_ERROR(gpu_buffer_free(mBuff));
    //Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
    //ALLOCATE HOST AND DEVICE BUFFER
    CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
    CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));
    // Re-map the buffer layout with the new size
    gpu_bpm_align_init_buffer_(bpmBuffer, averageQuerySize, candidatesPerQuery);
  }
}

/************************************************************
Functions to send & process a BPM buffer to GPU
************************************************************/

gpu_error_t gpu_bpm_align_reorder_process(const gpu_bpm_align_queries_buffer_t* const qry, const gpu_bpm_align_candidates_buffer_t* const cand,
                                    	    gpu_scheduler_buffer_t* const rebuff, gpu_bpm_align_cigars_buffer_t* const res)
{
  const uint32_t numBuckets = GPU_BPM_ALIGN_NUM_BUCKETS_FOR_BINNING;
  uint32_t idBucket, idCandidate;
  uint32_t numThreadsPerQuery, numQueriesPerWarp;
  uint32_t tmpBuckets[numBuckets], numCandidatesPerBucket[numBuckets], numWarpsPerBucket[numBuckets];
  // Init buckets (32 buckets => max 4096 bases)
  rebuff->numWarps = 0;
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
    numCandidatesPerBucket[idBucket]      = 0;
    numWarpsPerBucket[idBucket]           = 0;
  }
  // Fill buckets with elements per bucket
  for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
    idBucket = (qry->h_qinfo[cand->h_candidatesInfo[idCandidate].idQuery].size - 1) / GPU_BPM_ALIGN_PEQ_LENGTH_PER_CUDA_THREAD;
    idBucket = (idBucket < (rebuff->numBuckets - 1)) ? idBucket : (rebuff->numBuckets - 1);
    numCandidatesPerBucket[idBucket]++;
  }
  // Number of warps per bucket
  rebuff->threadMapScheduler.elementsPerBuffer = 0;
  for(idBucket = 0; idBucket < rebuff->numBuckets - 1; idBucket++){
    numThreadsPerQuery = idBucket + 1;
    numQueriesPerWarp = GPU_WARP_SIZE / numThreadsPerQuery;
    numWarpsPerBucket[idBucket] = GPU_DIV_CEIL(numCandidatesPerBucket[idBucket], numQueriesPerWarp);
    rebuff->h_initPosPerBucket[idBucket] = rebuff->threadMapScheduler.elementsPerBuffer;
    rebuff->h_endPosPerBucket[idBucket] = rebuff->threadMapScheduler.elementsPerBuffer + numCandidatesPerBucket[idBucket];
    rebuff->threadMapScheduler.elementsPerBuffer += numWarpsPerBucket[idBucket] * numQueriesPerWarp;
  }
  // Fill the start position warps for each bucket
  for(idBucket = 1; idBucket < rebuff->numBuckets; idBucket++)
    rebuff->h_initWarpPerBucket[idBucket] = rebuff->h_initWarpPerBucket[idBucket-1] + numWarpsPerBucket[idBucket-1];
  // Reorder by size the candidates
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++)
    tmpBuckets[idBucket] = rebuff->h_initPosPerBucket[idBucket];
  for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
    idBucket = (qry->h_qinfo[cand->h_candidatesInfo[idCandidate].idQuery].size - 1) / GPU_BPM_ALIGN_PEQ_LENGTH_PER_CUDA_THREAD;
    if (idBucket < (rebuff->numBuckets - 1)){ // Sanity check (discards binning on outlayer buckets)
      rebuff->threadMapScheduler.h_reorderBuffer[tmpBuckets[idBucket]] = idCandidate;
      tmpBuckets[idBucket]++;
    }
  }
  // Calculate the number of warps necessaries in the GPU
  for(idBucket = 0; idBucket < (rebuff->numBuckets - 1); idBucket++)
    rebuff->numWarps += numWarpsPerBucket[idBucket];
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_align_reordering_buffer(gpu_buffer_t *mBuff)
{
  const uint32_t                           		 binSize  = mBuff->data.abpm.queryBinSize;
  const bool                               		 binning  = mBuff->data.abpm.queryBinning;
  const gpu_bpm_align_queries_buffer_t* const    qry      = &mBuff->data.abpm.queries;
  const gpu_bpm_align_candidates_buffer_t* const cand     = &mBuff->data.abpm.candidates;
  gpu_bpm_align_cigars_buffer_t* const           res      = &mBuff->data.abpm.cigars;
  gpu_scheduler_buffer_t* const                  rebuff   = &mBuff->data.abpm.reorderBuffer;
  uint32_t                           			 idBucket;
  //Re-initialize the reorderBuffer (to reuse the buffer)
  rebuff->numBuckets                             = GPU_BPM_ALIGN_NUM_BUCKETS_FOR_BINNING;
  rebuff->threadMapScheduler.elementsPerBuffer   = 0;
  rebuff->numWarps                               = 0;
  //Initialize buckets (32 buckets => max 4096 bases)
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
    rebuff->h_initPosPerBucket[idBucket]  = 0;
    rebuff->h_initWarpPerBucket[idBucket] = 0;
    rebuff->h_endPosPerBucket[idBucket]   = 0;
  }
  if(binning){
    GPU_ERROR(gpu_bpm_align_reorder_process(qry, cand, rebuff, res));
  }else{
    //Calculate the number of warps necessaries in the GPU
    rebuff->numWarps = GPU_DIV_CEIL(binSize * cand->numCandidates, GPU_WARP_SIZE);
    //Fill the current bucket
    rebuff->h_initPosPerBucket[binSize - 1]  = 0;
    rebuff->h_initWarpPerBucket[binSize - 1] = 0;
    rebuff->h_endPosPerBucket[binSize - 1]   = cand->numCandidates;
    //Fill the start warp_id & candidate for each bucket
    for(idBucket = binSize; idBucket < rebuff->numBuckets; idBucket++){
      rebuff->h_initPosPerBucket[idBucket]  = cand->numCandidates;
      rebuff->h_initWarpPerBucket[idBucket] = rebuff->numWarps;
      rebuff->h_endPosPerBucket[idBucket]   = cand->numCandidates;
    }
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_align_transfer_CPU_to_GPU(gpu_buffer_t *mBuff)
{
  gpu_bpm_align_queries_buffer_t    *qry      = &mBuff->data.abpm.queries;
  gpu_bpm_align_candidates_buffer_t *cand     = &mBuff->data.abpm.candidates;
  gpu_scheduler_buffer_t            *rebuff   = &mBuff->data.abpm.reorderBuffer;
  gpu_bpm_align_cigars_buffer_t     *res      = &mBuff->data.abpm.cigars;
  cudaStream_t                      idStream  = mBuff->listStreams[mBuff->idStream];
  size_t                            cpySize   = 0;
  float                             bufferUtilization;
  // Defining buffer offsets
  cpySize += qry->totalQueriesPEQs * sizeof(gpu_bpm_align_peq_entry_t);
  cpySize += qry->totalQueriesBases * sizeof(gpu_bpm_align_qry_entry_t);
  cpySize += qry->numQueries * sizeof(gpu_bpm_align_qry_info_t);
  cpySize += cand->numCandidates * sizeof(gpu_bpm_align_cand_info_t);
  cpySize += rebuff->threadMapScheduler.elementsPerBuffer * sizeof(uint32_t);
  cpySize += rebuff->numBuckets * sizeof(uint32_t);
  cpySize += rebuff->numBuckets * sizeof(uint32_t);
  cpySize += rebuff->numBuckets * sizeof(uint32_t);
  cpySize += res->numCigars * sizeof(gpu_bpm_align_cigar_info_t);
  cpySize += res->numCigarEntries * sizeof(gpu_bpm_align_cigar_entry_t);
  bufferUtilization = (double)cpySize / (double)mBuff->sizeBuffer;
  // Compacting transferences with high buffer occupation
  if(bufferUtilization > 0.15){
    cpySize  = ((void *) (res->d_cigarsInfo + res->numCigars)) - ((void *) qry->d_queries);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_queries, qry->h_queries, cpySize, cudaMemcpyHostToDevice, idStream));
  }else{
    // Transfer Binary Queries to GPU
    cpySize = qry->totalQueriesBases * sizeof(gpu_bpm_align_qry_entry_t);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_queries, qry->h_queries, cpySize, cudaMemcpyHostToDevice, idStream));
    // Transfer to GPU the information associated with raw Queries
    cpySize = qry->totalQueriesPEQs * sizeof(gpu_bpm_align_peq_entry_t);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_peq, qry->h_peq, cpySize, cudaMemcpyHostToDevice, idStream));
    // Transfer to GPU the information associated with query information
    cpySize = qry->numQueries * sizeof(gpu_bpm_align_qry_info_t);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_qinfo, qry->h_qinfo, cpySize, cudaMemcpyHostToDevice, idStream));
    // Transfer Candidates to GPU info
    cpySize = cand->numCandidates * sizeof(gpu_bpm_align_cand_info_t);
    CUDA_ERROR(cudaMemcpyAsync(cand->d_candidatesInfo, cand->h_candidatesInfo, cpySize, cudaMemcpyHostToDevice, idStream));
    // Transfer reordered buffer to GPU
    cpySize = rebuff->threadMapScheduler.elementsPerBuffer * sizeof(uint32_t);
    CUDA_ERROR(cudaMemcpyAsync(rebuff->threadMapScheduler.d_reorderBuffer, rebuff->threadMapScheduler.h_reorderBuffer, cpySize, cudaMemcpyHostToDevice, idStream));
    // Transfer bucket information to GPU
    cpySize = rebuff->numBuckets * sizeof(uint32_t);
    CUDA_ERROR(cudaMemcpyAsync(rebuff->d_initPosPerBucket, rebuff->h_initPosPerBucket, cpySize, cudaMemcpyHostToDevice, idStream));
    CUDA_ERROR(cudaMemcpyAsync(rebuff->d_initWarpPerBucket, rebuff->h_initWarpPerBucket, cpySize, cudaMemcpyHostToDevice, idStream));
    CUDA_ERROR(cudaMemcpyAsync(rebuff->d_endPosPerBucket, rebuff->h_endPosPerBucket, cpySize, cudaMemcpyHostToDevice, idStream));
    // Transfer intermediate cigar info
    cpySize = res->numCigars * sizeof(gpu_bpm_align_cigar_info_t);
    CUDA_ERROR(cudaMemcpyAsync(res->d_cigarsInfo, res->h_cigarsInfo, cpySize, cudaMemcpyHostToDevice, idStream));
  }
  // Succeed
  return (SUCCESS);
}

//Adjust the host buffer layout (output)
gpu_error_t gpu_bpm_align_transfer_GPU_to_CPU(gpu_buffer_t *mBuff)
{
  cudaStream_t                    idStream  =  mBuff->listStreams[mBuff->idStream];
  gpu_bpm_align_cigars_buffer_t   *res      = &mBuff->data.abpm.cigars;
  size_t                          cpySize;
  // Avoiding transferences of the intermediate results (binning input work regularization)
  cpySize = ((void *) (res->d_cigars + res->numCigarEntries)) - ((void *) res->d_cigarsInfo);
  CUDA_ERROR(cudaMemcpyAsync(res->h_cigarsInfo, res->d_cigarsInfo, cpySize, cudaMemcpyDeviceToHost, idStream));
  // Succeed
  return (SUCCESS);
}

void gpu_bpm_align_send_buffer_(void* const bpmBuffer, const uint32_t numPEQEntries, const uint32_t numQueryBases,
                                const uint32_t numQueries, const uint32_t numCandidates, const uint32_t queryBinSize)
{
  gpu_buffer_t* const mBuff                      = (gpu_buffer_t *) bpmBuffer;
  const uint32_t    idSupDevice                  = mBuff->idSupportedDevice;
  uint32_t          idCandidate, numCigarEntries = 0;
  const uint32_t    bucketPaddingCandidates      = gpu_bpm_align_candidates_for_binning_padding();
  // Set real size of the things
  mBuff->data.abpm.queryBinning                  = !GPU_ISPOW2(queryBinSize);
  mBuff->data.abpm.queryBinSize                  = mBuff->data.abpm.queryBinning ? 0 : queryBinSize;
  mBuff->data.abpm.queries.totalQueriesPEQs      = numPEQEntries;
  mBuff->data.abpm.queries.totalQueriesBases     = numQueryBases;
  mBuff->data.abpm.queries.numQueries            = numQueries;
  mBuff->data.abpm.candidates.numCandidates      = numCandidates;
  mBuff->data.abpm.cigars.numCigars              = numCandidates;
  // Mapping the maximum number of task re-scheduled by the number of threads.
  mBuff->data.abpm.reorderBuffer.threadMapScheduler.elementsPerBuffer = numCandidates + bucketPaddingCandidates;
  // Setting real output number of elements
  for(idCandidate = 0; idCandidate < numCandidates; idCandidate++){
	const uint32_t idQuery = mBuff->data.abpm.candidates.h_candidatesInfo[idCandidate].idQuery;
	mBuff->data.abpm.cigars.h_cigarsInfo[idCandidate].offsetCigarStart = numCigarEntries;
    numCigarEntries += mBuff->data.abpm.queries.h_qinfo[idQuery].size + 1;
  }
  mBuff->data.abpm.cigars.numCigarEntries = numCigarEntries;
  // Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
  // CPU->GPU Transfers & Process Kernel in Asynchronous way
  GPU_ERROR(gpu_bpm_align_reordering_buffer(mBuff));
  GPU_ERROR(gpu_bpm_align_transfer_CPU_to_GPU(mBuff));
  /* INCLUDED SUPPORT for future GPUs with PTX ASM code (JIT compiling) */
  GPU_ERROR(gpu_bpm_align_process_buffer(mBuff));
  // GPU->CPU Transfers
  GPU_ERROR(gpu_bpm_align_transfer_GPU_to_CPU(mBuff));
}

/************************************************************
Functions to receive & process a BPM buffer from GPU
************************************************************/

void gpu_bpm_align_receive_buffer_(void* const bpmBuffer)
{
  gpu_buffer_t* const mBuff       = (gpu_buffer_t *) bpmBuffer;
  const uint32_t      idSupDevice = mBuff->idSupportedDevice;
  const cudaStream_t  idStream    = mBuff->listStreams[mBuff->idStream];
  // Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
  // Synchronize Stream (the thread wait for the commands done in the stream)
  CUDA_ERROR(cudaStreamSynchronize(idStream));
}

#endif /* GPU_BPM_PRIMITIVES_ALIGN_C_ */
