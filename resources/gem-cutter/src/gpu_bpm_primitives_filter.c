/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_BPM_PRIMITIVES_FILTER_C_
#define GPU_BPM_PRIMITIVES_FILTER_C_

#include "../include/gpu_bpm_primitives.h"

/************************************************************
Functions to get the GPU BPM buffers
************************************************************/

uint32_t gpu_bpm_filter_buffer_get_max_peq_entries_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.fbpm.maxPEQEntries);
}

uint32_t gpu_bpm_filter_buffer_get_max_candidates_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.fbpm.maxCandidates);
}

uint32_t gpu_bpm_filter_buffer_get_max_queries_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.fbpm.maxQueries);
}

gpu_bpm_filter_qry_entry_t* gpu_bpm_filter_buffer_get_peq_entries_(const void* const bpmBuffer){
  const gpu_buffer_t * const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.fbpm.queries.h_queries);
}

gpu_bpm_filter_cand_info_t* gpu_bpm_filter_buffer_get_candidates_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.fbpm.candidates.h_candidates);
}

gpu_bpm_filter_qry_info_t* gpu_bpm_filter_buffer_get_peq_info_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.fbpm.queries.h_qinfo);
}

gpu_bpm_filter_alg_entry_t* gpu_bpm_filter_buffer_get_alignments_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.fbpm.alignments.h_alignments);
}


/************************************************************
Functions to init all the BPM resources
************************************************************/

float gpu_bpm_filter_size_per_candidate(const uint32_t averageQuerySize, const uint32_t candidatesPerQuery)
{
  const size_t averageNumPEQEntries = GPU_DIV_CEIL(averageQuerySize, GPU_BPM_FILTER_PEQ_ENTRY_LENGTH);
  const size_t bytesPerQuery        = averageNumPEQEntries * sizeof(gpu_bpm_filter_qry_entry_t) + sizeof(gpu_bpm_filter_qry_info_t);
  const size_t bytesCandidate       = sizeof(gpu_bpm_filter_cand_info_t);
  const size_t bytesResult          = sizeof(gpu_bpm_filter_alg_entry_t);
  const size_t bytesReorderResult   = sizeof(gpu_bpm_filter_alg_entry_t);
  const size_t bytesBinningProcess  = sizeof(uint32_t) * GPU_BPM_FILTER_PENDING_NUM_TASK_LIST;
  const size_t bytesScheduling      = sizeof(uint32_t) * GPU_SCHEDULER_NUM_REORDER_BUFFERS;
  const size_t bytesCutOffProcess   = sizeof(gpu_bpm_filter_cand_error_t);
  // Calculate the necessary bytes for each BPM filter operation
  return((bytesPerQuery/(float)candidatesPerQuery) + bytesCandidate + bytesScheduling
      + bytesCutOffProcess + bytesResult + bytesReorderResult + bytesBinningProcess);
}

uint32_t gpu_bpm_filter_candidates_for_binning_padding()
{
  uint32_t idBucket;
  uint32_t bucketPaddingCandidates = 0;
  /* Worst number of dummy candidates added for padding the binning */
  for(idBucket = 1; idBucket < GPU_BPM_FILTER_NUM_BUCKETS_FOR_BINNING-1; ++idBucket)
    bucketPaddingCandidates += (GPU_WARP_SIZE / idBucket);
  /* Increased bytes per buffer taking account the padding*/
  return(bucketPaddingCandidates);
}

void gpu_bpm_filter_reallocate_host_buffer_layout(gpu_buffer_t* const mBuff)
{
  void* rawAlloc = mBuff->h_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.fbpm.queries.h_queries = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.queries.h_queries + mBuff->data.fbpm.maxPEQEntries);
  mBuff->data.fbpm.queries.h_qinfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.queries.h_qinfo + mBuff->data.fbpm.maxQueries);
  mBuff->data.fbpm.candidates.h_candidates = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.candidates.h_candidates + mBuff->data.fbpm.maxCandidates);
  mBuff->data.fbpm.reorderBuffer.threadMapScheduler.h_reorderBuffer = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.reorderBuffer.threadMapScheduler.h_reorderBuffer + mBuff->data.fbpm.maxReorderBuffer);
  mBuff->data.fbpm.reorderBuffer.h_initPosPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.reorderBuffer.h_initPosPerBucket + mBuff->data.fbpm.maxBuckets);
  mBuff->data.fbpm.reorderBuffer.h_initWarpPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.reorderBuffer.h_initWarpPerBucket + mBuff->data.fbpm.maxBuckets);
  //Adjust the host buffer layout (output)
  mBuff->data.fbpm.alignments.h_reorderAlignments = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.alignments.h_reorderAlignments + mBuff->data.fbpm.maxReorderBuffer);
  mBuff->data.fbpm.alignments.h_alignments = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.alignments.h_alignments + mBuff->data.fbpm.maxAlignments);
  //Adjust the host buffer layout (intermediate-temporal)
  mBuff->data.fbpm.cutoff.h_error = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.cutoff.h_error + mBuff->data.fbpm.maxCutOffEntries);
  mBuff->data.fbpm.cutoff.h_pendingTasks_tmpSpace[0] = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.cutoff.h_pendingTasks_tmpSpace[0] + mBuff->data.fbpm.maxPendingTasks);
  mBuff->data.fbpm.cutoff.h_pendingTasks_tmpSpace[1] = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.cutoff.h_pendingTasks_tmpSpace[1] + mBuff->data.fbpm.maxPendingTasks);
  mBuff->data.fbpm.reorderBuffer.taskMapScheduler.h_reorderBuffer = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.reorderBuffer.taskMapScheduler.h_reorderBuffer + mBuff->data.fbpm.maxCandidates);
}

void gpu_bpm_filter_reallocate_device_buffer_layout(gpu_buffer_t* const mBuff)
{
  void* rawAlloc = mBuff->d_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.fbpm.queries.d_queries = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.queries.d_queries + mBuff->data.fbpm.maxPEQEntries);
  mBuff->data.fbpm.queries.d_qinfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.queries.d_qinfo + mBuff->data.fbpm.maxQueries);
  mBuff->data.fbpm.candidates.d_candidates = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.candidates.d_candidates + mBuff->data.fbpm.maxCandidates);
  mBuff->data.fbpm.reorderBuffer.threadMapScheduler.d_reorderBuffer = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.reorderBuffer.threadMapScheduler.d_reorderBuffer + mBuff->data.fbpm.maxReorderBuffer);
  mBuff->data.fbpm.reorderBuffer.d_initPosPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.reorderBuffer.d_initPosPerBucket + mBuff->data.fbpm.maxBuckets);
  mBuff->data.fbpm.reorderBuffer.d_initWarpPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.reorderBuffer.d_initWarpPerBucket + mBuff->data.fbpm.maxBuckets);
  //Adjust the host buffer layout (output)
  mBuff->data.fbpm.alignments.d_reorderAlignments = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.alignments.d_reorderAlignments + mBuff->data.fbpm.maxReorderBuffer);
  mBuff->data.fbpm.alignments.d_alignments = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.fbpm.alignments.d_alignments + mBuff->data.fbpm.maxAlignments);
}

void gpu_bpm_filter_init_buffer_(void* const bpmBuffer, const uint32_t averageQuerySize, const uint32_t candidatesPerQuery)
{
  gpu_buffer_t* const mBuff                   = (gpu_buffer_t *) bpmBuffer;
  const double        sizeBuff                = mBuff->sizeBuffer * 0.95;
  const size_t        averageNumPEQEntries    = GPU_DIV_CEIL(averageQuerySize, GPU_BPM_FILTER_PEQ_ENTRY_LENGTH);
  const uint32_t      numInputs               = (uint32_t)(sizeBuff / gpu_bpm_filter_size_per_candidate(averageQuerySize, candidatesPerQuery));
  const uint32_t      maxCandidates           = numInputs - gpu_bpm_filter_candidates_for_binning_padding();
  const uint32_t      bucketPaddingCandidates = gpu_bpm_filter_candidates_for_binning_padding();
  //set the type of the buffer
  mBuff->typeBuffer = GPU_BPM_FILTER;
  //Set real size of the input
  mBuff->data.fbpm.maxCandidates    = maxCandidates;
  mBuff->data.fbpm.maxAlignments    = maxCandidates;
  mBuff->data.fbpm.maxPEQEntries    = (maxCandidates / candidatesPerQuery) * averageNumPEQEntries;
  mBuff->data.fbpm.maxQueries       = (maxCandidates / candidatesPerQuery);
  mBuff->data.fbpm.maxReorderBuffer = maxCandidates + bucketPaddingCandidates;
  mBuff->data.fbpm.maxCutOffEntries = maxCandidates;
  mBuff->data.fbpm.maxPendingTasks  = maxCandidates;
  mBuff->data.fbpm.maxBuckets       = GPU_BPM_FILTER_NUM_BUCKETS_FOR_BINNING;
  // Set the corresponding buffer layout
  gpu_bpm_filter_reallocate_host_buffer_layout(mBuff);
  gpu_bpm_filter_reallocate_device_buffer_layout(mBuff);
}

void gpu_bpm_filter_init_and_realloc_buffer_(void* const bpmBuffer, const uint32_t totalPEQEntries, const uint32_t totalCandidates, const uint32_t totalQueries)
{
  // Buffer re-initialization
  gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  const uint32_t averarageNumPEQEntries = totalPEQEntries / totalQueries;
  const uint32_t averageQuerySize       = (totalPEQEntries * GPU_BPM_FILTER_PEQ_ENTRY_LENGTH) / totalQueries;
  const uint32_t candidatesPerQuery     = totalCandidates / totalQueries;
  // Remap the buffer layout with new information trying to fit better
  gpu_bpm_filter_init_buffer_(bpmBuffer, averageQuerySize, candidatesPerQuery);
  // Checking if we need to reallocate a bigger buffer
  if( (totalPEQEntries > gpu_bpm_filter_buffer_get_max_peq_entries_(bpmBuffer)) ||
      (totalCandidates > gpu_bpm_filter_buffer_get_max_candidates_(bpmBuffer))  ||
      (totalQueries    > gpu_bpm_filter_buffer_get_max_queries_(bpmBuffer))){
    // Resize the GPU buffer to fit the required input
    const uint32_t  idSupDevice             = mBuff->idSupportedDevice;
    const float     resizeFactor            = 2.0;
    const size_t    bytesPerBPMBuffer       = totalCandidates * gpu_bpm_filter_size_per_candidate(averarageNumPEQEntries,candidatesPerQuery);
    //Recalculate the minimum buffer size
    //printf("RESIZE[BPM_FILTER] %d %d \n",  mBuff->sizeBuffer, bytesPerBPMBuffer * resizeFactor);
    mBuff->sizeBuffer = bytesPerBPMBuffer * resizeFactor;
    //FREE HOST AND DEVICE BUFFER
    GPU_ERROR(gpu_buffer_free(mBuff));
    //Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
    //ALLOCATE HOST AND DEVICE BUFFER
    CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
    CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));
    // Re-map the buffer layout with the new size
    gpu_bpm_filter_init_buffer_(bpmBuffer, averageQuerySize, candidatesPerQuery);
  }
}

/************************************************************
Functions to send & process a BPM buffer to GPU
************************************************************/


gpu_error_t gpu_bpm_filter_reorder_process_cutoff(gpu_buffer_t* const mBuff)
{
  const gpu_bpm_filter_queries_buffer_t* const    qry     = &mBuff->data.fbpm.queries;
  const gpu_bpm_filter_candidates_buffer_t* const cand    = &mBuff->data.fbpm.candidates;
  gpu_bpm_filter_alignments_buffer_t* const 	  res     = &mBuff->data.fbpm.alignments;
  gpu_scheduler_buffer_t* const             	  rebuff  = &mBuff->data.fbpm.reorderBuffer;
  gpu_bpm_filter_cutoff_buffer_t* const           cutoff  = &mBuff->data.fbpm.cutoff;
  const uint32_t                                  idKey   = mBuff->data.fbpm.cutoff.minTileDepth;
  // Initializing local iterators
  const uint32_t numBuckets = GPU_BPM_FILTER_NUM_BUCKETS_FOR_BINNING;
  uint32_t idBucket, idBuff, idPendingTask, idNextPendingTask = 0;
  uint32_t numThreadsPerQuery, numQueriesPerWarp;
  uint32_t tmpBuckets[numBuckets], numCandidatesPerBucket[numBuckets], numWarpsPerBucket[numBuckets];
  uint32_t numElementsTaskMapped = 0, numElementsThreadMapped = 0;
  // Initializing buckets (32 buckets => max 4096 bases)
  rebuff->numWarps = 0;
  // Initializing the work scheduler
  rebuff->taskMapScheduler.elementsPerBuffer   = 0;
  rebuff->threadMapScheduler.elementsPerBuffer = 0;
  // Initializing the thread scheduler histograms
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
    numCandidatesPerBucket[idBucket]      = 0;
    numWarpsPerBucket[idBucket]           = 0;
    tmpBuckets[idBucket]                  = 0;
  }
  // Fill buckets with elements per bucket (counting process)
  for(idPendingTask = 0; idPendingTask < cutoff->numCurrPendingTasks; idPendingTask++){
	uint32_t idCandidate = cutoff->h_currPendingTasks[idPendingTask];
	// Discard the candidate
	if((res->h_alignments[idCandidate].score  != GPU_BPM_FILTER_SCORE_INF) &&
	   (res->h_alignments[idCandidate].column != GPU_BPM_FILTER_SCORE_INF)){
	  const uint32_t idQuery = cand->h_candidates[idCandidate].query;
	  if(qry->h_qinfo[idQuery].idTile < idKey){
	    // Check if the task it was filtered
		idBucket = (qry->h_qinfo[idQuery].tileSize - 1) / GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD;
		idBucket = GPU_MIN(idBucket, rebuff->numBuckets - 1);
		numCandidatesPerBucket[idBucket]++;
	  }else{
	    //Deferring task due to dependencies for the next kernel launch
	    cutoff->h_nextPendingTasks[idNextPendingTask] = idCandidate;
	    idNextPendingTask++;
	  }
	}
  }
  // Number of warps per bucket
  for(idBucket = 0; idBucket < rebuff->numBuckets - 1; idBucket++){
    numThreadsPerQuery = idBucket + 1;
    numQueriesPerWarp = GPU_WARP_SIZE / numThreadsPerQuery;
    numWarpsPerBucket[idBucket] = GPU_DIV_CEIL(numCandidatesPerBucket[idBucket], numQueriesPerWarp);
    rebuff->h_initPosPerBucket[idBucket] = numElementsThreadMapped;
    numElementsThreadMapped += numWarpsPerBucket[idBucket] * numQueriesPerWarp;
  }
  // Fill the start position warps for each bucket
  for(idBucket = 1; idBucket < rebuff->numBuckets; idBucket++)
    rebuff->h_initWarpPerBucket[idBucket] = rebuff->h_initWarpPerBucket[idBucket-1] + numWarpsPerBucket[idBucket-1];
  // Allocate buffer (candidates)
  for(idBuff = 0; idBuff < numElementsThreadMapped; idBuff++)
    rebuff->threadMapScheduler.h_reorderBuffer[idBuff] = GPU_UINT32_ONES;
  // Set the number of real results in the reorder buffer
  res->numReorderedAlignments = numElementsThreadMapped;
  // Reorder by size the candidates
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++)
    tmpBuckets[idBucket] = rebuff->h_initPosPerBucket[idBucket];
  // Filling the thread and task schedulers with the corresponding task
  for(idPendingTask = 0; idPendingTask < cutoff->numCurrPendingTasks; idPendingTask++){
	uint32_t idCandidate = cutoff->h_currPendingTasks[idPendingTask];
    if((qry->h_qinfo[cand->h_candidates[idCandidate].query].idTile < idKey) &&
       (res->h_alignments[idCandidate].score != GPU_BPM_FILTER_SCORE_INF)   &&
       (res->h_alignments[idCandidate].column != GPU_BPM_FILTER_SCORE_INF)){
      idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].tileSize - 1) / GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD;
      // Discards the tiles larger than the pre-processed bin sizes
      if (idBucket < (rebuff->numBuckets - 1)){
        const uint32_t idReorderedCandidate = tmpBuckets[idBucket];
        rebuff->threadMapScheduler.h_reorderBuffer[idReorderedCandidate] = idCandidate;
        tmpBuckets[idBucket]++;
        rebuff->taskMapScheduler.h_reorderBuffer[numElementsTaskMapped]  = idReorderedCandidate;
        numElementsTaskMapped++;
      }
    }
  }
  // Fill the paddings with replicated candidates (always the last candidate)
  for(idBuff = 0; idBuff < numElementsThreadMapped; idBuff++){
    if(rebuff->threadMapScheduler.h_reorderBuffer[idBuff] == GPU_UINT32_ONES){
    	rebuff->threadMapScheduler.h_reorderBuffer[idBuff] = rebuff->threadMapScheduler.h_reorderBuffer[idBuff-1];
    }
  }
  // Calculate the number of warps necessaries in the GPU
  for(idBucket = 0; idBucket < (rebuff->numBuckets - 1); idBucket++)
    rebuff->numWarps += numWarpsPerBucket[idBucket];
  // Updating the task schedulers
  rebuff->taskMapScheduler.elementsPerBuffer   = numElementsTaskMapped;
  rebuff->threadMapScheduler.elementsPerBuffer = numElementsThreadMapped;
  cutoff->numNextPendingTasks                  = idNextPendingTask;
  // Swap temporary buffers
  uint32_t* tmpPending = cutoff->h_currPendingTasks;
  cutoff->numCurrPendingTasks = cutoff->numNextPendingTasks;
  cutoff->h_currPendingTasks = cutoff->h_nextPendingTasks;
  cutoff->h_nextPendingTasks = tmpPending;
  // Succeed
  return (SUCCESS);
}


gpu_error_t gpu_bpm_filter_reorder_process_full(gpu_buffer_t* const mBuff)
{
  //Wrapping data structures
  const gpu_bpm_filter_queries_buffer_t* const    qry            = &mBuff->data.fbpm.queries;
  const gpu_bpm_filter_candidates_buffer_t* const cand           = &mBuff->data.fbpm.candidates;
  gpu_bpm_filter_alignments_buffer_t* const 	  res            = &mBuff->data.fbpm.alignments;
  gpu_scheduler_buffer_t* const             	  rebuff         = &mBuff->data.fbpm.reorderBuffer;
  //Initializing local data structures
  const uint32_t numBuckets = GPU_BPM_FILTER_NUM_BUCKETS_FOR_BINNING;
  uint32_t  idBucket, idCandidate, idBuff;
  uint32_t  numThreadsPerQuery, numQueriesPerWarp;
  uint32_t  tmpBuckets[numBuckets], numCandidatesPerBucket[numBuckets], numWarpsPerBucket[numBuckets];
  uint32_t  elementsPerBuffer = 0;
  uint32_t* reorderBuffer = rebuff->threadMapScheduler.h_reorderBuffer;
  // Initializing buckets (32 buckets => max 4096 bases)
  rebuff->numWarps = 0;
  // Initializing the work scheduler
  rebuff->taskMapScheduler.elementsPerBuffer   = 0;
  rebuff->threadMapScheduler.elementsPerBuffer = 0;
  // Initializing the thread scheduler histograms
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
    numCandidatesPerBucket[idBucket]      = 0;
    numWarpsPerBucket[idBucket]           = 0;
    tmpBuckets[idBucket]                  = 0;
  }
  // Fill buckets with elements per bucket
  for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
    idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].tileSize - 1) / GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD;
    idBucket = (idBucket < (rebuff->numBuckets - 1)) ? idBucket : (rebuff->numBuckets - 1);
    numCandidatesPerBucket[idBucket]++;
  }
  // Number of warps per bucket
  for(idBucket = 0; idBucket < rebuff->numBuckets - 1; idBucket++){
    numThreadsPerQuery = idBucket + 1;
    numQueriesPerWarp = GPU_WARP_SIZE / numThreadsPerQuery;
    numWarpsPerBucket[idBucket] = GPU_DIV_CEIL(numCandidatesPerBucket[idBucket], numQueriesPerWarp);
    rebuff->h_initPosPerBucket[idBucket] = elementsPerBuffer;
    elementsPerBuffer += numWarpsPerBucket[idBucket] * numQueriesPerWarp;
  }
  // Fill the start position warps for each bucket
  for(idBucket = 1; idBucket < rebuff->numBuckets; idBucket++)
    rebuff->h_initWarpPerBucket[idBucket] = rebuff->h_initWarpPerBucket[idBucket-1] + numWarpsPerBucket[idBucket-1];
  // Allocate buffer (candidates)
  for(idBuff = 0; idBuff < elementsPerBuffer; idBuff++)
    reorderBuffer[idBuff] = GPU_UINT32_ONES;
  // Set the number of real results in the reorder buffer
  res->numReorderedAlignments = elementsPerBuffer;
  // Reorder by size the candidates
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++)
    tmpBuckets[idBucket] = rebuff->h_initPosPerBucket[idBucket];
  for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
    idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].tileSize - 1) / GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD;
    if (idBucket < (rebuff->numBuckets - 1)){
      reorderBuffer[tmpBuckets[idBucket]] = idCandidate;
      tmpBuckets[idBucket]++;
    }
  }
  // Fill the paddings with replicated candidates (always the last candidate)
  for(idBuff = 0; idBuff < elementsPerBuffer; idBuff++)
    if(reorderBuffer[idBuff] == GPU_UINT32_ONES) reorderBuffer[idBuff] = reorderBuffer[idBuff-1];
  // Calculate the number of warps necessaries in the GPU
  for(idBucket = 0; idBucket < (rebuff->numBuckets - 1); idBucket++)
    rebuff->numWarps += numWarpsPerBucket[idBucket];
  // Updating the task schedulers
  rebuff->threadMapScheduler.elementsPerBuffer = elementsPerBuffer;
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_reordering_buffer_init(gpu_scheduler_buffer_t* const rebuff)
{
  // Declaration of data iterators
  uint32_t idBucket                            = 0;
  // Re-initialize the reorderBuffer (to reuse the buffer)
  rebuff->numBuckets                           = GPU_BPM_FILTER_NUM_BUCKETS_FOR_BINNING;
  rebuff->threadMapScheduler.elementsPerBuffer = 0;
  rebuff->taskMapScheduler.elementsPerBuffer   = 0;
  rebuff->numWarps                             = 0;
  // Initialize buckets (32 buckets => max 4096 bases)
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
	rebuff->h_initPosPerBucket[idBucket]  = 0;
	rebuff->h_initWarpPerBucket[idBucket] = 0;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_reorder_process_light(gpu_buffer_t* const mBuff)
{
  const uint32_t                                  binSize = mBuff->data.fbpm.queryBinSize;
  const gpu_bpm_filter_candidates_buffer_t* const cand    = &mBuff->data.fbpm.candidates;
  gpu_scheduler_buffer_t* const                   rebuff  = &mBuff->data.fbpm.reorderBuffer;
  uint32_t idBucket;
  // Calculate the number of warps necessaries in the GPU
  rebuff->numWarps = GPU_DIV_CEIL(binSize * cand->numCandidates, GPU_WARP_SIZE);
  // Fill the start warp_id & candidate for each bucket
  for(idBucket = binSize; idBucket < rebuff->numBuckets; idBucket++){
    rebuff->h_initPosPerBucket[idBucket]  = cand->numCandidates;
    rebuff->h_initWarpPerBucket[idBucket] = rebuff->numWarps;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_reordering_buffer(gpu_buffer_t* const mBuff)
{
  const bool                       activeBinning  = mBuff->data.fbpm.queryBinning;
  const bool                       activeCutoff   = mBuff->data.fbpm.activeCutOff;
  gpu_scheduler_buffer_t* const    rebuff         = &mBuff->data.fbpm.reorderBuffer;
  // Re-initialize the reorderBuffer (to reuse the buffer)
  GPU_ERROR(gpu_bpm_filter_reordering_buffer_init(rebuff));
  // Avoiding the internal rearranging data when is not necessary
  if(activeCutoff){
    GPU_ERROR(gpu_bpm_filter_reorder_process_cutoff(mBuff));
  }else if(activeBinning){
    GPU_ERROR(gpu_bpm_filter_reorder_process_full(mBuff));
  }else{
    GPU_ERROR(gpu_bpm_filter_reorder_process_light(mBuff));
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff)
{
  const gpu_bpm_filter_queries_buffer_t* const     qry      = &mBuff->data.fbpm.queries;
  const gpu_bpm_filter_candidates_buffer_t* const  cand     = &mBuff->data.fbpm.candidates;
  const gpu_scheduler_buffer_t* const              rebuff   = &mBuff->data.fbpm.reorderBuffer;
  const gpu_bpm_filter_alignments_buffer_t* const  res      = &mBuff->data.fbpm.alignments;
  const cudaStream_t                               idStream = mBuff->listStreams[mBuff->idStream];
  size_t                                           cpySize  = 0;
  float                                            bufferUtilization;
  // Input data: Size of Queries and Candidates
  cpySize += qry->totalQueriesEntries * sizeof(gpu_bpm_filter_qry_entry_t);
  cpySize += qry->numQueries * sizeof(gpu_bpm_filter_qry_info_t);
  cpySize += cand->numCandidates * sizeof(gpu_bpm_filter_cand_info_t);
  // Intermediate data: Size of taskMapScheduler
  cpySize += cand->numCandidates * sizeof(uint32_t);
  // Intermediate data: Size of threadMapScheduler
  cpySize += cand->numCandidates * sizeof(uint32_t);
  // Intermediate data: Histograms for the GPU thread schedulers
  cpySize += rebuff->numBuckets * sizeof(uint32_t);
  cpySize += rebuff->numBuckets * sizeof(uint32_t);
  // Intermediate data: Accumulators for tile error rates
  cpySize += cand->numCandidates * sizeof(gpu_bpm_filter_cand_error_t);
  // Output data: Final Alignments and Rearranged alignments
  cpySize += res->numAlignments * sizeof(gpu_bpm_filter_alg_entry_t);
  cpySize += res->numAlignments * sizeof(gpu_bpm_filter_alg_entry_t);
  bufferUtilization = (double)cpySize / (double)mBuff->sizeBuffer;
  // Compacting transferences with high buffer occupation
  if(bufferUtilization > 0.15){
    if(mBuff->data.fbpm.activeCutOff){
      cpySize  = ((void *) (cand->d_candidates + cand->numCandidates)) - ((void *) qry->d_queries);
    }else{
      cpySize  = ((void *) (rebuff->d_initWarpPerBucket + rebuff->numBuckets)) - ((void *) qry->d_queries);
    }
    CUDA_ERROR(cudaMemcpyAsync(qry->d_queries, qry->h_queries, cpySize, cudaMemcpyHostToDevice, idStream));
  }else{
    // Transfer Binary Queries to GPU
    cpySize = qry->totalQueriesEntries * sizeof(gpu_bpm_filter_qry_entry_t);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_queries, qry->h_queries, cpySize, cudaMemcpyHostToDevice, idStream));
    // Transfer to GPU the information associated with Binary Queries
    cpySize = qry->numQueries * sizeof(gpu_bpm_filter_qry_info_t);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_qinfo, qry->h_qinfo, cpySize, cudaMemcpyHostToDevice, idStream));
    // Transfer Candidates to GPU
    cpySize = cand->numCandidates * sizeof(gpu_bpm_filter_cand_info_t);
    CUDA_ERROR(cudaMemcpyAsync(cand->d_candidates, cand->h_candidates, cpySize, cudaMemcpyHostToDevice, idStream));
    if(!mBuff->data.fbpm.activeCutOff){
      // Transfer reordered buffer to GPU
      cpySize = rebuff->threadMapScheduler.elementsPerBuffer * sizeof(uint32_t);
      CUDA_ERROR(cudaMemcpyAsync(rebuff->threadMapScheduler.d_reorderBuffer, rebuff->threadMapScheduler.h_reorderBuffer, cpySize, cudaMemcpyHostToDevice, idStream));
      // Transfer bucket information to GPU
      cpySize = rebuff->numBuckets * sizeof(uint32_t);
      CUDA_ERROR(cudaMemcpyAsync(rebuff->d_initPosPerBucket, rebuff->h_initPosPerBucket, cpySize, cudaMemcpyHostToDevice, idStream));
      CUDA_ERROR(cudaMemcpyAsync(rebuff->d_initWarpPerBucket, rebuff->h_initWarpPerBucket, cpySize, cudaMemcpyHostToDevice, idStream));
    }
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff)
{
  const cudaStream_t                              idStream =  mBuff->listStreams[mBuff->idStream];
  const gpu_bpm_filter_alignments_buffer_t* const res      = &mBuff->data.fbpm.alignments;
  // Avoiding transferences of the intermediate results (corner case: binning input work regularization)
  if(mBuff->data.fbpm.queryBinning){
    const size_t cpySize = res->numReorderedAlignments * sizeof(gpu_bpm_filter_alg_entry_t);
    CUDA_ERROR(cudaMemcpyAsync(res->h_reorderAlignments, res->d_reorderAlignments, cpySize, cudaMemcpyDeviceToHost, idStream));
  }else{
    const size_t cpySize = res->numAlignments * sizeof(gpu_bpm_filter_alg_entry_t);
    CUDA_ERROR(cudaMemcpyAsync(res->h_alignments, res->d_alignments, cpySize, cudaMemcpyDeviceToHost, idStream));
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_intermediate_data_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff)
{
  const gpu_scheduler_buffer_t * const rebuff   = &mBuff->data.fbpm.reorderBuffer;
  const cudaStream_t                  idStream  =  mBuff->listStreams[mBuff->idStream];
  // Intermediate data transference
  const size_t cpySize  = ((void *) (rebuff->d_initWarpPerBucket + rebuff->numBuckets)) - ((void *) rebuff->threadMapScheduler.d_reorderBuffer);
  CUDA_ERROR(cudaMemcpyAsync(rebuff->threadMapScheduler.d_reorderBuffer, rebuff->threadMapScheduler.h_reorderBuffer, cpySize, cudaMemcpyHostToDevice, idStream));
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_intermediate_data_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff)
{
  cudaStream_t                	     idStream  =  mBuff->listStreams[mBuff->idStream];
  gpu_bpm_filter_alignments_buffer_t *res      = &mBuff->data.fbpm.alignments;
  // Intermediate data transference
  const size_t cpySize = res->numReorderedAlignments * sizeof(gpu_bpm_filter_alg_entry_t);
  CUDA_ERROR(cudaMemcpyAsync(res->h_reorderAlignments, res->d_reorderAlignments, cpySize, cudaMemcpyDeviceToHost, idStream));
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_update_tile_depth(gpu_buffer_t* const mBuff)
{
  gpu_bpm_filter_cutoff_buffer_t* const            cutoff     = &mBuff->data.fbpm.cutoff;
  const gpu_bpm_filter_queries_buffer_t* const     qry        = &mBuff->data.fbpm.queries;
  const gpu_bpm_filter_candidates_buffer_t* const  cand       = &mBuff->data.fbpm.candidates;
  const uint32_t                                   numBuckets = GPU_BPM_FILTER_CUTOFF_MAX_HIST_SIZE;
  uint32_t                                         minIdKey   = GPU_UINT32_ONES;
  uint32_t                                         maxIdKey   = 0;
  // Indexes for the iterators
  uint32_t idCandidate = 0, idBucket = 0;
  // Initialise the buckets
  for(idBucket = 0; idBucket < numBuckets; idBucket++){
    cutoff->tileDepthHist[idBucket] = 0;
  }
  // Fill buckets with elements per bucket
  for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
	  const uint32_t idQuery = cand->h_candidates[idCandidate].query;
    const uint32_t idTile = qry->h_qinfo[idQuery].idTile;
    const uint32_t idBucket = GPU_MIN(idTile,numBuckets);
    cutoff->tileDepthHist[idBucket]++;
    minIdKey = GPU_MIN(idBucket,minIdKey);
    maxIdKey = GPU_MAX(idBucket,maxIdKey);
  }
  if(cand->numCandidates) cutoff->pendingTasks = true;
    else cutoff->pendingTasks = false;
  cutoff->minTileDepth = minIdKey;
  cutoff->maxTileDepth = maxIdKey + 1;
  // Succeed
  return (SUCCESS);
}


gpu_error_t gpu_bpm_filter_init_work(gpu_buffer_t* const mBuff)
{
  const gpu_bpm_filter_candidates_buffer_t* const  cand   = &mBuff->data.fbpm.candidates;
  gpu_bpm_filter_cutoff_buffer_t* const            cutoff = &mBuff->data.fbpm.cutoff;
  uint32_t idTask = 0;
  // Initializing the work search space
  gpu_bpm_filter_update_tile_depth(mBuff);
  // Setting number of pending tasks
  cutoff->numCurrPendingTasks = cand->numCandidates;
  // Initializing memory space for pending lists
  cutoff->h_currPendingTasks = cutoff->h_pendingTasks_tmpSpace[0];
  cutoff->h_nextPendingTasks = cutoff->h_pendingTasks_tmpSpace[1];
  // Assigning the tasks to pending
  for(idTask = 0; idTask < cutoff->numCurrPendingTasks; idTask++){
    cutoff->h_currPendingTasks[idTask] = idTask;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_schedule_work(gpu_buffer_t* const mBuff)
{
  gpu_bpm_filter_cutoff_buffer_t* const cutoff = &mBuff->data.fbpm.cutoff;
  const int32_t  minIdKey          = cutoff->minTileDepth;
  const int32_t  maxIdKey          = cutoff->maxTileDepth;
  // Setting the minimun amount of work to offload (tradeoff between synch/launch overhead & overall work)
  const uint32_t minTilesPerLaunch = 16 * GPU_BPM_FILTER_MIN_TILES_PER_LAUNCH;
  // Indexes for the iterators
  int32_t  minIdBucket = minIdKey, maxIdBucket = maxIdKey;
  uint32_t totalTiles = 0;
  // Search for the minimal amount of work
  while((totalTiles <= minTilesPerLaunch) && (minIdBucket < maxIdKey)){
    totalTiles += cutoff->tileDepthHist[minIdBucket];
    minIdBucket++;
  }
  // Update the pending work
  if(totalTiles != 0) cutoff->pendingTasks = true;
    else cutoff->pendingTasks = false;
  // Update Keys
  cutoff->minTileDepth = minIdBucket;
  cutoff->maxTileDepth = maxIdBucket;
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_init_alignments(gpu_buffer_t* const mBuff)
{
  // Avoiding transferences of the intermediate results (binning input work regularization)
  gpu_bpm_filter_alignments_buffer_t *res    = &mBuff->data.fbpm.alignments;
  const gpu_bpm_filter_alg_entry_t initAlignment = {0,0};
  uint32_t idRes;
  // Reverting the original input organization
  for(idRes = 0; idRes < res->numAlignments; idRes++)
    res->h_alignments[idRes] = initAlignment;
  // Succeed
  return (SUCCESS);
}

void gpu_bpm_filter_send_buffer_(void* const bpmBuffer, const uint32_t numPEQEntries, const uint32_t numQueries,
                          	  	 const uint32_t numCandidates, const uint32_t maxQuerySize, const uint32_t queryBinSize)
{
  gpu_buffer_t* const mBuff                          = (gpu_buffer_t *) bpmBuffer;
  const uint32_t      idSupDevice                    = mBuff->idSupportedDevice;
  // Set real size of the internal data
  mBuff->data.fbpm.queryBinning                      = !GPU_ISPOW2(queryBinSize);
  mBuff->data.fbpm.maxQuerySize                      = maxQuerySize;
  mBuff->data.fbpm.activeCutOff						 = (GPU_BPM_FILTER_CUTOFF_ACTIVE && (maxQuerySize > GPU_BPM_FILTER_CUTOFF_MAX_TILE_LENGHT));
  mBuff->data.fbpm.queryBinSize                      = mBuff->data.fbpm.queryBinning ? 0 : queryBinSize;
  mBuff->data.fbpm.queries.totalQueriesEntries       = numPEQEntries;
  mBuff->data.fbpm.queries.numQueries                = numQueries;
  mBuff->data.fbpm.candidates.numCandidates          = numCandidates;
  mBuff->data.fbpm.alignments.numAlignments          = numCandidates;
  // ReorderAlignments elements are allocated just for divergent size queries
  mBuff->data.fbpm.alignments.numReorderedAlignments = 0;
  // Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
  // Inspect if all queries have 1 or more tiles and initialise cutoff
  if(mBuff->data.fbpm.activeCutOff){
	  // Turn active the cutoff if user specifies
      GPU_ERROR(gpu_bpm_filter_update_tile_depth(mBuff));
	  // Turn active the cutoff if user specifies
	  mBuff->data.fbpm.activeCutOff	= GPU_BPM_FILTER_CUTOFF_ACTIVE;
	  // Cut-off techniques requires the bining strategy
	  mBuff->data.fbpm.queryBinning	= GPU_BPM_FILTER_CUTOFF_ACTIVE;
  }else{
	  // Disabling cutoff keys, we want to process all the work in parallel
	  mBuff->data.fbpm.cutoff.minTileDepth = GPU_BPM_FILTER_CUTOFF_DISABLED_KEY;
	  mBuff->data.fbpm.cutoff.maxTileDepth = GPU_BPM_FILTER_CUTOFF_DISABLED_KEY;
	  // Turning cutoff disabled, processing will not require the extra overhead
	  mBuff->data.fbpm.activeCutOff	= false;
  }
  // Process the filtering alignments (checking cutoff strategy)
  if(mBuff->data.fbpm.activeCutOff){
    // Clean & initialize the output buffer
    GPU_ERROR(gpu_bpm_filter_init_alignments(mBuff));
    // CPU->GPU Transfers & Process Kernel in Asynchronous way
	GPU_ERROR(gpu_bpm_filter_transfer_CPU_to_GPU(mBuff));
	// Initializing the queue of pending tasks
	GPU_ERROR(gpu_bpm_filter_init_work(mBuff));
	// Enabling and updating the cutoff keys
	GPU_ERROR(gpu_bpm_filter_schedule_work(mBuff));
	// Iterating to reduce the amount of tiles
	while(mBuff->data.fbpm.cutoff.pendingTasks){
	  // Generating the amount of necessary work
	  GPU_ERROR(gpu_bpm_filter_reordering_buffer(mBuff));
	  // Included support for future GPUs with PTX ASM code (JIT compiling)
	  GPU_ERROR(gpu_bpm_filter_intermediate_data_transfer_CPU_to_GPU(mBuff));
	  GPU_ERROR(gpu_bpm_filter_process_buffer(mBuff));
	  GPU_ERROR(gpu_bpm_filter_intermediate_data_transfer_GPU_to_CPU(mBuff));
	  // Host-Device synchronization
	  GPU_ERROR(gpu_bpm_filter_device_synch(mBuff));
	  // Post-processing (re-arrange output scores and cutoff the alignment work)
	  GPU_ERROR(gpu_bpm_filter_reordering_alignments_cutoff(mBuff));
      // Updating the cutoff keys
	  GPU_ERROR(gpu_bpm_filter_schedule_work(mBuff));
	}
  }else{
	// CPU->GPU Transfers & Process Kernel in Asynchronous way
	GPU_ERROR(gpu_bpm_filter_reordering_buffer(mBuff));
	GPU_ERROR(gpu_bpm_filter_transfer_CPU_to_GPU(mBuff));
	// Included support for future GPUs with PTX ASM code (JIT compiling)
	GPU_ERROR(gpu_bpm_filter_process_buffer(mBuff));
	// GPU->CPU Transfers
	GPU_ERROR(gpu_bpm_filter_transfer_GPU_to_CPU(mBuff));
  }
}

/************************************************************
Functions to receive & process a BPM buffer from GPU
************************************************************/

gpu_error_t gpu_bpm_filter_reordering_alignments(gpu_buffer_t* const mBuff)
{
  // Avoiding transferences of the intermediate results (binning input work regularization)
  if(mBuff->data.fbpm.queryBinning){
    gpu_scheduler_buffer_t             *rebuff = &mBuff->data.fbpm.reorderBuffer;
    gpu_bpm_filter_alignments_buffer_t *res    = &mBuff->data.fbpm.alignments;
    uint32_t idRes;
    // Reverting the original input organization
    for(idRes = 0; idRes < res->numReorderedAlignments; idRes++)
      res->h_alignments[rebuff->threadMapScheduler.h_reorderBuffer[idRes]] = res->h_reorderAlignments[idRes];
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_invalidate_all_chained_tiles(gpu_buffer_t* const mBuff, const uint32_t idCandidate)
{
  // Getting data structures
  const gpu_bpm_filter_queries_buffer_t* const     qry      = &mBuff->data.fbpm.queries;
  const gpu_bpm_filter_candidates_buffer_t* const  cand     = &mBuff->data.fbpm.candidates;
  const gpu_bpm_filter_alignments_buffer_t* const  res      = &mBuff->data.fbpm.alignments;
  gpu_bpm_filter_cutoff_buffer_t* const            cutoff   = &mBuff->data.fbpm.cutoff;
  // Getting data limits
  const gpu_bpm_filter_alg_entry_t invalidatedTile = {GPU_BPM_FILTER_SCORE_INF,GPU_BPM_FILTER_SCORE_INF};
  const uint32_t numCandidates = cand->numCandidates;
  // Setting iterators
  uint32_t prevIdCandidate = idCandidate;
  uint32_t currIdCandidate = prevIdCandidate + 1;
  // Checking if there are more chained tiles to invalidate
  uint32_t prevIdQuery = cand->h_candidates[prevIdCandidate].query;
  uint32_t prevIdTile = qry->h_qinfo[prevIdQuery].idTile;
  // Invalidating all the chained tiles from the same candidate
  while((currIdCandidate < numCandidates) && (prevIdTile < qry->h_qinfo[cand->h_candidates[currIdCandidate].query].idTile)){
    // Getting the tile identification
	const uint32_t currIdQuery = cand->h_candidates[currIdCandidate].query;
	const uint32_t currIdTile  = qry->h_qinfo[currIdQuery].idTile;
	// Invalidate the next tile
	res->h_alignments[currIdCandidate] = invalidatedTile;
	// Update the amount of scheduled work
	cutoff->tileDepthHist[currIdTile]--;
	// Update iterator information
	prevIdTile = currIdTile;
	currIdCandidate++;
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_reordering_alignments_cutoff(gpu_buffer_t* const mBuff)
{
  const gpu_bpm_filter_queries_buffer_t* const     qry      = &mBuff->data.fbpm.queries;
  const gpu_bpm_filter_candidates_buffer_t* const  cand     = &mBuff->data.fbpm.candidates;
  const gpu_scheduler_buffer_t* const              rebuff   = &mBuff->data.fbpm.reorderBuffer;
  const gpu_bpm_filter_alignments_buffer_t* const  res      = &mBuff->data.fbpm.alignments;
  gpu_bpm_filter_cutoff_buffer_t* const            cutoff   = &mBuff->data.fbpm.cutoff;
  //Setting constant special values for filtered tiles
  const gpu_bpm_filter_alg_entry_t invalidatedTile = {GPU_BPM_FILTER_SCORE_INF,GPU_BPM_FILTER_SCORE_INF};
  uint32_t idRes = 0;
  // Reverting the original permutation of the output
  for(idRes = 0; idRes < rebuff->taskMapScheduler.elementsPerBuffer; idRes++){
    const uint32_t idTask         = rebuff->taskMapScheduler.h_reorderBuffer[idRes];
	  const uint32_t idCandidate    = rebuff->threadMapScheduler.h_reorderBuffer[idTask];
	  const uint32_t idQuery        = cand->h_candidates[idCandidate].query;
	const uint32_t idTile         = qry->h_qinfo[idQuery].idTile;
	const uint32_t maxErrorTile   = qry->h_qinfo[idQuery].tileMaxError;
	const uint32_t maxErrorChain  = qry->h_qinfo[idQuery].chainMaxError;
	const uint32_t candErrorTile  = res->h_reorderAlignments[idTask].score;
	// Update and accumulate the chained score error along all candidate tiles
	uint32_t candErrorAcc         = candErrorTile;
	if(idTile != 0){
	  // Avoiding corner case for the first tile
	  candErrorAcc += cutoff->h_error[idCandidate - 1];
	}
	//Check if this candidate has been previously invalidated
	if((res->h_alignments[idCandidate].score != GPU_BPM_FILTER_SCORE_INF) &&
	   (res->h_alignments[idCandidate].column != GPU_BPM_FILTER_SCORE_INF)){
	  // Update the current tile-chain accumulated error
	  cutoff->h_error[idCandidate] = candErrorAcc;
	  // Update the final values
	  if((candErrorTile > maxErrorTile) || (candErrorAcc > maxErrorChain)){
	    // Invalidating the current tile
	    res->h_alignments[idCandidate] = invalidatedTile;
	    // Invalidating the next chained tiles belonging the same candidate & query
	    GPU_ERROR(gpu_bpm_filter_invalidate_all_chained_tiles(mBuff, idCandidate));
	  }else{
		// Saving the current filtered BPM result
	    res->h_alignments[idCandidate] = res->h_reorderAlignments[idTask];
	  }
	  //Remove work from the scheduler
	    cutoff->tileDepthHist[idTile]--;
	  }
  }
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_bpm_filter_discard_alignments(gpu_buffer_t* const mBuff)
{
  // Avoiding transferences of the intermediate results (binning input work regularization)
  if(mBuff->data.fbpm.queryBinning){
    gpu_scheduler_buffer_t*             rebuff = &mBuff->data.fbpm.reorderBuffer;
    gpu_bpm_filter_alignments_buffer_t* res    = &mBuff->data.fbpm.alignments;
    uint32_t idRes;
    // Reverting the original input organization
    for(idRes = 0; idRes < res->numReorderedAlignments; idRes++)
      res->h_alignments[rebuff->threadMapScheduler.h_reorderBuffer[idRes]] = res->h_reorderAlignments[idRes];
  }
  // Succeed
  return (SUCCESS);
}

void gpu_bpm_filter_receive_buffer_(void* const bpmBuffer)
{
  gpu_buffer_t* const mBuff       = (gpu_buffer_t *) bpmBuffer;
  const uint32_t      idSupDevice = mBuff->idSupportedDevice;
  const cudaStream_t  idStream    = mBuff->listStreams[mBuff->idStream];
  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
  if(!mBuff->data.fbpm.activeCutOff){
    //Synchronize Stream (the thread wait for the commands done in the stream)
    CUDA_ERROR(cudaStreamSynchronize(idStream));
    //Reorder the final results
    GPU_ERROR(gpu_bpm_filter_reordering_alignments(mBuff));
  }
}

gpu_error_t gpu_bpm_filter_device_synch(gpu_buffer_t* const mBuff)
{
  const uint32_t      idSupDevice = mBuff->idSupportedDevice;
  const cudaStream_t  idStream    = mBuff->listStreams[mBuff->idStream];
  // Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
  // Synchronize Stream (the thread wait for the commands done in the stream)
  CUDA_ERROR(cudaStreamSynchronize(idStream));
  // Succeed
  return (SUCCESS);
}

#endif /* GPU_BPM_PRIMITIVES_FILTER_C_ */
