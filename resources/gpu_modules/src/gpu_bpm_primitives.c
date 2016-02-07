/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Host scheduler for BPM on GPU
 */

#ifndef GPU_BPM_PRIMITIVES_C_
#define GPU_BPM_PRIMITIVES_C_

#include "../include/gpu_bpm_primitives.h"

/************************************************************
Functions to get the GPU BPM buffers
************************************************************/

uint32_t gpu_bpm_buffer_get_max_peq_entries_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.bpm.maxPEQEntries);
}

uint32_t gpu_bpm_buffer_get_max_candidates_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.bpm.maxCandidates);
}

uint32_t gpu_bpm_buffer_get_max_queries_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.bpm.maxQueries);
}

gpu_bpm_qry_entry_t* gpu_bpm_buffer_get_peq_entries_(const void* const bpmBuffer){
  const gpu_buffer_t * const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.bpm.queries.h_queries);
}

gpu_bpm_cand_info_t* gpu_bpm_buffer_get_candidates_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.bpm.candidates.h_candidates);
}

gpu_bpm_qry_info_t* gpu_bpm_buffer_get_peq_info_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.bpm.queries.h_qinfo);
}

gpu_bpm_alg_entry_t* gpu_bpm_buffer_get_alignments_(const void* const bpmBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  return(mBuff->data.bpm.alignments.h_alignments);
}


/************************************************************
Functions to init all the BPM resources
************************************************************/

float gpu_bpm_size_per_candidate(const uint32_t averageQuerySize, const uint32_t candidatesPerQuery)
{
  const size_t averageNumPEQEntries = GPU_DIV_CEIL(averageQuerySize, GPU_BPM_PEQ_ENTRY_LENGTH);
  const size_t bytesPerQuery        = averageNumPEQEntries * sizeof(gpu_bpm_qry_entry_t) + sizeof(gpu_bpm_qry_info_t);
  const size_t bytesCandidate       = sizeof(gpu_bpm_cand_info_t);
  const size_t bytesResult          = sizeof(gpu_bpm_alg_entry_t);
  const size_t bytesReorderResult   = sizeof(gpu_bpm_alg_entry_t);
  const size_t bytesBinningProcess  = sizeof(uint32_t);

  return((bytesPerQuery/(float)candidatesPerQuery) + bytesCandidate
      + bytesResult + bytesReorderResult + bytesBinningProcess);
}

uint32_t gpu_bpm_candidates_for_binning_padding()
{
  uint32_t idBucket;
  uint32_t bucketPaddingCandidates = 0;

  /* Worst number of dummy candidates added for padding the binning */
  for(idBucket = 1; idBucket < GPU_BPM_NUM_BUCKETS_FOR_BINNING-1; ++idBucket)
    bucketPaddingCandidates += (GPU_WARP_SIZE / idBucket);

  /* Increased bytes per buffer taking account the padding*/
  return(bucketPaddingCandidates);
}

void gpu_bpm_reallocate_host_buffer_layout(gpu_buffer_t* mBuff)
{
  void* rawAlloc = mBuff->h_rawData;

  //Adjust the host buffer layout (input)
  mBuff->data.bpm.queries.h_queries = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.queries.h_queries + mBuff->data.bpm.maxPEQEntries);
  mBuff->data.bpm.queries.h_qinfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.queries.h_qinfo + mBuff->data.bpm.maxQueries);
  mBuff->data.bpm.candidates.h_candidates = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.candidates.h_candidates + mBuff->data.bpm.maxCandidates);
  mBuff->data.bpm.reorderBuffer.h_reorderBuffer = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.reorderBuffer.h_reorderBuffer + mBuff->data.bpm.maxReorderBuffer);
  mBuff->data.bpm.reorderBuffer.h_initPosPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.reorderBuffer.h_initPosPerBucket + mBuff->data.bpm.maxBuckets);
  mBuff->data.bpm.reorderBuffer.h_initWarpPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.reorderBuffer.h_initWarpPerBucket + mBuff->data.bpm.maxBuckets);

  //Adjust the host buffer layout (output)
  mBuff->data.bpm.alignments.h_reorderAlignments = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.alignments.h_reorderAlignments + mBuff->data.bpm.maxReorderBuffer);
  mBuff->data.bpm.alignments.h_alignments = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.alignments.h_alignments + mBuff->data.bpm.maxAlignments);
}

void gpu_bpm_reallocate_device_buffer_layout(gpu_buffer_t* mBuff)
{
  void* rawAlloc = mBuff->d_rawData;

  //Adjust the host buffer layout (input)
  mBuff->data.bpm.queries.d_queries = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.queries.d_queries + mBuff->data.bpm.maxPEQEntries);
  mBuff->data.bpm.queries.d_qinfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.queries.d_qinfo + mBuff->data.bpm.maxQueries);
  mBuff->data.bpm.candidates.d_candidates = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.candidates.d_candidates + mBuff->data.bpm.maxCandidates);
  mBuff->data.bpm.reorderBuffer.d_reorderBuffer = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.reorderBuffer.d_reorderBuffer + mBuff->data.bpm.maxReorderBuffer);
  mBuff->data.bpm.reorderBuffer.d_initPosPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.reorderBuffer.d_initPosPerBucket + mBuff->data.bpm.maxBuckets);
  mBuff->data.bpm.reorderBuffer.d_initWarpPerBucket = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.reorderBuffer.d_initWarpPerBucket + mBuff->data.bpm.maxBuckets);

  //Adjust the host buffer layout (output)
  mBuff->data.bpm.alignments.d_reorderAlignments = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.alignments.d_reorderAlignments + mBuff->data.bpm.maxReorderBuffer);
  mBuff->data.bpm.alignments.d_alignments = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.bpm.alignments.d_alignments + mBuff->data.bpm.maxAlignments);
}

void gpu_bpm_init_buffer_(void* const bpmBuffer, const uint32_t averageQuerySize, const uint32_t candidatesPerQuery)
{
  gpu_buffer_t* const mBuff                   = (gpu_buffer_t *) bpmBuffer;
  const double        sizeBuff                = mBuff->sizeBuffer * 0.95;
  const size_t        averageNumPEQEntries    = GPU_DIV_CEIL(averageQuerySize, GPU_BPM_PEQ_ENTRY_LENGTH);
  const uint32_t      numInputs               = (uint32_t)(sizeBuff / gpu_bpm_size_per_candidate(averageQuerySize, candidatesPerQuery));
  const uint32_t      maxCandidates           = numInputs - gpu_bpm_candidates_for_binning_padding();
  const uint32_t      bucketPaddingCandidates = gpu_bpm_candidates_for_binning_padding();

  //set the type of the buffer
  mBuff->typeBuffer = GPU_BPM;

  //Set real size of the input
  mBuff->data.bpm.maxCandidates    = maxCandidates;
  mBuff->data.bpm.maxAlignments    = maxCandidates;
  mBuff->data.bpm.maxPEQEntries    = (maxCandidates / candidatesPerQuery) * averageNumPEQEntries;
  mBuff->data.bpm.maxQueries       = (maxCandidates / candidatesPerQuery);
  mBuff->data.bpm.maxReorderBuffer = maxCandidates + bucketPaddingCandidates;
  mBuff->data.bpm.maxBuckets       = GPU_BPM_NUM_BUCKETS_FOR_BINNING;

  gpu_bpm_reallocate_host_buffer_layout(mBuff);
  gpu_bpm_reallocate_device_buffer_layout(mBuff);
}

void gpu_bpm_init_and_realloc_buffer_(void *bpmBuffer, const uint32_t totalPEQEntries, const uint32_t totalCandidates, const uint32_t totalQueries)
{
  gpu_buffer_t* const mBuff = (gpu_buffer_t *) bpmBuffer;
  const uint32_t averarageNumPEQEntries = totalPEQEntries / totalQueries;
  const uint32_t averageQuerySize       = (totalPEQEntries * GPU_BPM_PEQ_ENTRY_LENGTH) / totalQueries;
  const uint32_t candidatesPerQuery     = totalCandidates / totalQueries;

  gpu_bpm_init_buffer_(bpmBuffer, averageQuerySize, candidatesPerQuery);

  if( (totalPEQEntries > gpu_bpm_buffer_get_max_peq_entries_(bpmBuffer)) &&
      (totalCandidates > gpu_bpm_buffer_get_max_candidates_(bpmBuffer))  &&
      (totalQueries    > gpu_bpm_buffer_get_max_queries_(bpmBuffer))){
    // Resize the GPU buffer to fit the required input
    const uint32_t  idSupDevice             = mBuff->idSupportedDevice;
    const float     resizeFactor            = 2.0;
    const size_t    bytesPerBPMBuffer       = totalCandidates * gpu_bpm_size_per_candidate(averarageNumPEQEntries,candidatesPerQuery);

    //Recalculate the minimum buffer size
    mBuff->sizeBuffer = bytesPerBPMBuffer * resizeFactor;

    //FREE HOST AND DEVICE BUFFER
    GPU_ERROR(gpu_free_buffer(mBuff));

    //Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

    //ALLOCATE HOST AND DEVICE BUFFER
    CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
    CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));

    gpu_bpm_init_buffer_(bpmBuffer, averageQuerySize, candidatesPerQuery);
  }
}

/************************************************************
Functions to send & process a BPM buffer to GPU
************************************************************/

gpu_error_t gpu_bpm_reorder_process(const gpu_bpm_queries_buffer_t* const qry, const gpu_bpm_candidates_buffer_t* const cand,
                                    gpu_bpm_reorder_buffer_t* const rebuff, gpu_bpm_alignments_buffer_t* const res)
{
  const uint32_t numBuckets = GPU_BPM_NUM_BUCKETS_FOR_BINNING;
  uint32_t idBucket, idCandidate, idBuff;
  uint32_t numThreadsPerQuery, numQueriesPerWarp;
  uint32_t tmpBuckets[numBuckets], numCandidatesPerBucket[numBuckets], numWarpsPerBucket[numBuckets];

  //Init buckets (32 buckets => max 4096 bases)
  rebuff->numWarps = 0;
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
    numCandidatesPerBucket[idBucket]      = 0;
    numWarpsPerBucket[idBucket]           = 0;
    tmpBuckets[idBucket]                  = 0;
  }

  //Fill buckets with elements per bucket
  for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
    idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].size - 1) / GPU_BPM_PEQ_LENGTH_PER_CUDA_THREAD;
    idBucket = (idBucket < (rebuff->numBuckets - 1)) ? idBucket : (rebuff->numBuckets - 1);
    numCandidatesPerBucket[idBucket]++;
  }

  //Number of warps per bucket
  rebuff->candidatesPerBuffer = 0;
  for(idBucket = 0; idBucket < rebuff->numBuckets - 1; idBucket++){
    numThreadsPerQuery = idBucket + 1;
    numQueriesPerWarp = GPU_WARP_SIZE / numThreadsPerQuery;
    numWarpsPerBucket[idBucket] = GPU_DIV_CEIL(numCandidatesPerBucket[idBucket], numQueriesPerWarp);
    rebuff->h_initPosPerBucket[idBucket] = rebuff->candidatesPerBuffer;
    rebuff->candidatesPerBuffer += numWarpsPerBucket[idBucket] * numQueriesPerWarp;
  }

  //Fill the start position warps for each bucket
  for(idBucket = 1; idBucket < rebuff->numBuckets; idBucket++)
    rebuff->h_initWarpPerBucket[idBucket] = rebuff->h_initWarpPerBucket[idBucket-1] + numWarpsPerBucket[idBucket-1];

  //Allocate buffer (candidates)
  for(idBuff = 0; idBuff < rebuff->candidatesPerBuffer; idBuff++)
    rebuff->h_reorderBuffer[idBuff] = GPU_UINT32_ONES;

  //Set the number of real results in the reorder buffer
  res->numReorderedAlignments = rebuff->candidatesPerBuffer;

  //Reorder by size the candidates
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++)
    tmpBuckets[idBucket] = rebuff->h_initPosPerBucket[idBucket];
  for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
    idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].size - 1) / GPU_BPM_PEQ_LENGTH_PER_CUDA_THREAD;
    if (idBucket < (rebuff->numBuckets - 1)){
      rebuff->h_reorderBuffer[tmpBuckets[idBucket]] = idCandidate;
      tmpBuckets[idBucket]++;
    }
  }

  //Fill the paddings with replicated candidates (always the last candidate)
  for(idBuff = 0; idBuff < rebuff->candidatesPerBuffer; idBuff++)
    if(rebuff->h_reorderBuffer[idBuff] == GPU_UINT32_ONES) rebuff->h_reorderBuffer[idBuff] = rebuff->h_reorderBuffer[idBuff-1];

  //Calculate the number of warps necessaries in the GPU
  for(idBucket = 0; idBucket < (rebuff->numBuckets - 1); idBucket++)
    rebuff->numWarps += numWarpsPerBucket[idBucket];

  return (SUCCESS);
}

gpu_error_t gpu_bpm_reordering_buffer(gpu_buffer_t *mBuff)
{
  const uint32_t                           binning  = mBuff->data.bpm.queryBinning;
  const gpu_bpm_queries_buffer_t* const    qry      = &mBuff->data.bpm.queries;
  const gpu_bpm_candidates_buffer_t* const cand     = &mBuff->data.bpm.candidates;
        gpu_bpm_alignments_buffer_t* const res      = &mBuff->data.bpm.alignments;
        gpu_bpm_reorder_buffer_t* const    rebuff   = &mBuff->data.bpm.reorderBuffer;
        uint32_t                           idBucket;

  //Re-initialize the reorderBuffer (to reuse the buffer)
  rebuff->numBuckets          = GPU_BPM_NUM_BUCKETS_FOR_BINNING;
  rebuff->candidatesPerBuffer = 0;
  rebuff->numWarps            = 0;

  //Initialize buckets (32 buckets => max 4096 bases)
  for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
    rebuff->h_initPosPerBucket[idBucket]  = 0;
    rebuff->h_initWarpPerBucket[idBucket] = 0;
  }

  if(binning){
    //Calculate the number of warps necessaries in the GPU
    rebuff->numWarps = GPU_DIV_CEIL(binning * cand->numCandidates, GPU_WARP_SIZE);
    //Fill the start warp_id & candidate for each bucket
    for(idBucket = binning; idBucket < rebuff->numBuckets; idBucket++){
      rebuff->h_initPosPerBucket[idBucket]  = cand->numCandidates;
      rebuff->h_initWarpPerBucket[idBucket] = rebuff->numWarps;
    }
  }else{
    GPU_ERROR(gpu_bpm_reorder_process(qry, cand, rebuff, res));
  }

  return (SUCCESS);
}

gpu_error_t gpu_bpm_transfer_CPU_to_GPU(gpu_buffer_t *mBuff)
{
  gpu_bpm_queries_buffer_t    *qry      = &mBuff->data.bpm.queries;
  gpu_bpm_candidates_buffer_t *cand     = &mBuff->data.bpm.candidates;
  gpu_bpm_reorder_buffer_t    *rebuff   = &mBuff->data.bpm.reorderBuffer;
  gpu_bpm_alignments_buffer_t *res      = &mBuff->data.bpm.alignments;
  cudaStream_t                idStream  = mBuff->idStream;
  size_t                      cpySize   = 0;
  float                       bufferUtilization;

  cpySize += qry->totalQueriesEntries * sizeof(gpu_bpm_qry_entry_t);
  cpySize += qry->numQueries * sizeof(gpu_bpm_qry_info_t);
  cpySize += cand->numCandidates * sizeof(gpu_bpm_cand_info_t);
  cpySize += rebuff->candidatesPerBuffer * sizeof(uint32_t);
  cpySize += rebuff->numBuckets * sizeof(uint32_t);
  cpySize += rebuff->numBuckets * sizeof(uint32_t);
  cpySize += res->numReorderedAlignments * sizeof(gpu_bpm_alg_entry_t);
  cpySize += res->numAlignments * sizeof(gpu_bpm_alg_entry_t);
  bufferUtilization = (double)cpySize / (double)mBuff->sizeBuffer;

  if(bufferUtilization > 0.25){
    cpySize  = ((void *) (rebuff->d_initWarpPerBucket + rebuff->numBuckets)) - ((void *) qry->d_queries);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_queries, qry->h_queries, cpySize, cudaMemcpyHostToDevice, idStream));
  }else{
    //Transfer Binary Queries to GPU
    cpySize = qry->totalQueriesEntries * sizeof(gpu_bpm_qry_entry_t);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_queries, qry->h_queries, cpySize, cudaMemcpyHostToDevice, idStream));

    //Transfer to GPU the information associated with Binary Queries
    cpySize = qry->numQueries * sizeof(gpu_bpm_qry_info_t);
    CUDA_ERROR(cudaMemcpyAsync(qry->d_qinfo, qry->h_qinfo, cpySize, cudaMemcpyHostToDevice, idStream));

    //Transfer Candidates to GPU
    cpySize = cand->numCandidates * sizeof(gpu_bpm_cand_info_t);
    CUDA_ERROR(cudaMemcpyAsync(cand->d_candidates, cand->h_candidates, cpySize, cudaMemcpyHostToDevice, idStream));

    //Transfer reordered buffer to GPU
    cpySize = rebuff->candidatesPerBuffer * sizeof(uint32_t);
    CUDA_ERROR(cudaMemcpyAsync(rebuff->d_reorderBuffer, rebuff->h_reorderBuffer, cpySize, cudaMemcpyHostToDevice, idStream));

    //Transfer bucket information to GPU
    cpySize = rebuff->numBuckets * sizeof(uint32_t);
    CUDA_ERROR(cudaMemcpyAsync(rebuff->d_initPosPerBucket, rebuff->h_initPosPerBucket, cpySize, cudaMemcpyHostToDevice, idStream));
    CUDA_ERROR(cudaMemcpyAsync(rebuff->d_initWarpPerBucket, rebuff->h_initWarpPerBucket, cpySize, cudaMemcpyHostToDevice, idStream));
  }

  return (SUCCESS);
}

gpu_error_t gpu_bpm_transfer_GPU_to_CPU(gpu_buffer_t *mBuff)
{
  cudaStream_t                idStream  =  mBuff->idStream;
  gpu_bpm_alignments_buffer_t *res      = &mBuff->data.bpm.alignments;
  size_t                      cpySize;

  if(mBuff->data.bpm.queryBinning){
    cpySize = res->numAlignments * sizeof(gpu_bpm_alg_entry_t);
    CUDA_ERROR(cudaMemcpyAsync(res->h_alignments, res->d_alignments, cpySize, cudaMemcpyDeviceToHost, idStream));
  }else{
    cpySize = res->numReorderedAlignments * sizeof(gpu_bpm_alg_entry_t);
    CUDA_ERROR(cudaMemcpyAsync(res->h_reorderAlignments, res->d_reorderAlignments, cpySize, cudaMemcpyDeviceToHost, idStream));
  }

  return (SUCCESS);
}

void gpu_bpm_send_buffer_(void* const bpmBuffer, const uint32_t numPEQEntries, const uint32_t numQueries,
                          const uint32_t numCandidates, const uint32_t queryBinning)
{
  gpu_buffer_t* const mBuff     = (gpu_buffer_t *) bpmBuffer;
  const uint32_t    idSupDevice = mBuff->idSupportedDevice;

  //Set real size of the things
  mBuff->data.bpm.queryBinning                = queryBinning;
  mBuff->data.bpm.queries.totalQueriesEntries = numPEQEntries;
  mBuff->data.bpm.queries.numQueries          = numQueries;
  mBuff->data.bpm.candidates.numCandidates    = numCandidates;
  mBuff->data.bpm.alignments.numAlignments    = numCandidates;

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  //CPU->GPU Transfers & Process Kernel in Asynchronous way
  GPU_ERROR(gpu_bpm_reordering_buffer(mBuff));
  GPU_ERROR(gpu_bpm_transfer_CPU_to_GPU(mBuff));

  /* INCLUDED SUPPORT for future GPUs with PTX ASM code (JIT compiling) */
  GPU_ERROR(gpu_bpm_process_buffer(mBuff));

  //GPU->CPU Transfers
  GPU_ERROR(gpu_bpm_transfer_GPU_to_CPU(mBuff));
}

/************************************************************
Functions to receive & process a BPM buffer from GPU
************************************************************/

gpu_error_t gpu_bpm_reordering_alignments(gpu_buffer_t *mBuff)
{
  if(mBuff->data.bpm.queryBinning == 0){
    gpu_bpm_reorder_buffer_t    *rebuff = &mBuff->data.bpm.reorderBuffer;
    gpu_bpm_alignments_buffer_t *res    = &mBuff->data.bpm.alignments;
    uint32_t idRes;

    for(idRes = 0; idRes < res->numReorderedAlignments; idRes++)
      res->h_alignments[rebuff->h_reorderBuffer[idRes]] = res->h_reorderAlignments[idRes];
  }
  return (SUCCESS);
}

void gpu_bpm_receive_buffer_(void* const bpmBuffer)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) bpmBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  //Synchronize Stream (the thread wait for the commands done in the stream)
  CUDA_ERROR(cudaStreamSynchronize(mBuff->idStream));
  //Reorder the final results
  GPU_ERROR(gpu_bpm_reordering_alignments(mBuff));
}

#endif /* GPU_BPM_PRIMITIVES_C_ */
