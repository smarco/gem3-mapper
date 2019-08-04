/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_PRIMITIVES_ASEARCH_C_
#define GPU_FMI_PRIMITIVES_ASEARCH_C_

#include "../include/gpu_fmi_primitives.h"
#include "../include/gpu_sa_primitives.h"

#ifdef GPU_FMI_DEBUG
//Data classified by num of regions
uint32_t histogram_queries[15]            = {0};
uint32_t histogram_query_candidates[15]   = {0};
//Data classified by region size
uint32_t histogram_regions[110]           = {0};
uint32_t histogram_region_candidates[110] = {0};
uint32_t histogram_coverage[110]          = {0};
#endif

/************************************************************
Functions to get the GPU FMI buffers
************************************************************/

gpu_fmi_search_query_t* gpu_fmi_asearch_buffer_get_queries_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.asearch.queries.h_queries);
}

gpu_fmi_search_query_info_t* gpu_fmi_asearch_buffer_get_queries_info_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.asearch.queries.h_queryInfo);
}

gpu_fmi_search_region_t* gpu_fmi_asearch_buffer_get_regions_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.asearch.queries.h_regions);
}

gpu_sa_search_inter_t* gpu_fmi_asearch_buffer_get_regions_intervals_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.asearch.regions.h_intervals);
}

gpu_fmi_search_region_info_t* gpu_fmi_asearch_buffer_get_regions_offsets_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.asearch.regions.h_regionsOffsets);
}

/************************************************************
Functions to debug internals
************************************************************/

#ifdef GPU_FMI_DEBUG
void gpu_buffer_fmi_asearch_process_histogram(const gpu_buffer_t* const mBuff){
  uint32_t idQuery = 0, numBases = 0, numQueries = mBuff->data.asearch.queries.numQueries;
  for(idQuery = 0; idQuery < numQueries; ++idQuery){
    //Getting query stats
    uint32_t idRegion = 0, numRegions = mBuff->data.asearch.queries.h_regions[idQuery].num_regions;
    uint32_t offset   = mBuff->data.asearch.queries.h_regions[idQuery].init_offset;
    for(idRegion = 0, numBases = 0; idRegion < numRegions; ++idRegion){
      //Getting region stats
      uint32_t hi  = mBuff->data.asearch.regions.h_intervals[offset + idRegion].hi;
      uint32_t low = mBuff->data.asearch.regions.h_intervals[offset + idRegion].low;
      uint32_t init_offset   = mBuff->data.asearch.regions.h_regionsOffsets[offset + idRegion].init_offset;
      uint32_t end_offset    = mBuff->data.asearch.regions.h_regionsOffsets[offset + idRegion].end_offset;
      //Setting region stats
      uint32_t numCandidates = hi - low;
      uint32_t regionSize    = end_offset - init_offset;
      numBases += regionSize;
      // Setting the num candidates per read
      histogram_query_candidates[numRegions] += numCandidates;
      // Setting the num candidates per region
      histogram_region_candidates[regionSize] += numCandidates;
      // Setting the region size
      histogram_regions[regionSize]++;
    }
    // Setting the num reads
    histogram_queries[numRegions]++;
    // Setting all regions coverage over the read
    histogram_coverage[numBases]++;
  }
}
#endif  

/************************************************************
Functions to get the maximum elements of the buffers
************************************************************/

uint32_t gpu_fmi_asearch_buffer_get_max_queries_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.asearch.numMaxQueries);
}

uint32_t gpu_fmi_asearch_buffer_get_max_regions_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.asearch.numMaxRegions);
}

uint32_t gpu_fmi_asearch_buffer_get_max_bases_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.asearch.numMaxBases);
}

/************************************************************
Functions to initialize the buffers (E. SEARCH)
************************************************************/

size_t gpu_fmi_asearch_size_per_query(const uint32_t averageQuerySize, const uint32_t averageRegionsPerQuery)
{
  //Memory size dedicated to each query
  const size_t bytesPerQueryRAW  = averageQuerySize * sizeof(gpu_fmi_search_query_t);
  const size_t bytesPerQueryInfo = sizeof(gpu_fmi_search_query_info_t) + sizeof(gpu_fmi_search_region_t);
  const size_t bytesPerQuery     = bytesPerQueryRAW + bytesPerQueryInfo;
  //Memory size dedicated to query regions (maxRegionsRatio means % of regions per query)
  const size_t bytesPerRegion    = sizeof(gpu_sa_search_inter_t) + sizeof(gpu_fmi_search_region_info_t);
  //Return maximum memory size required per each query
  return((averageRegionsPerQuery * bytesPerRegion) + bytesPerQuery);
}

void gpu_fmi_asearch_reallocate_host_buffer_layout(gpu_buffer_t* mBuff)
{
  const void* rawAlloc = mBuff->h_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.asearch.queries.h_queries = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.queries.h_queries + mBuff->data.asearch.numMaxBases);
  mBuff->data.asearch.queries.h_queryInfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.queries.h_queryInfo + mBuff->data.asearch.numMaxQueries);
  //Adjust the host buffer layout (output)
  mBuff->data.asearch.queries.h_regions = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.queries.h_regions + mBuff->data.asearch.numMaxQueries);
  mBuff->data.asearch.regions.h_intervals = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.regions.h_intervals + mBuff->data.asearch.numMaxRegions);
  mBuff->data.asearch.regions.h_regionsOffsets = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.regions.h_regionsOffsets + mBuff->data.asearch.numMaxRegions);
}

void gpu_fmi_asearch_reallocate_device_buffer_layout(gpu_buffer_t* mBuff)
{
  const void* rawAlloc = mBuff->d_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.asearch.queries.d_queries = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.queries.d_queries + mBuff->data.asearch.numMaxBases);
  mBuff->data.asearch.queries.d_queryInfo = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.queries.d_queryInfo + mBuff->data.asearch.numMaxQueries);
  //Adjust the host buffer layout (output)
  mBuff->data.asearch.queries.d_regions = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.queries.d_regions + mBuff->data.asearch.numMaxQueries);
  mBuff->data.asearch.regions.d_intervals = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.regions.d_intervals + mBuff->data.asearch.numMaxRegions);
  mBuff->data.asearch.regions.d_regionsOffsets = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.asearch.regions.d_regionsOffsets + mBuff->data.asearch.numMaxRegions);
}

void gpu_fmi_asearch_init_buffer_(void* const fmiBuffer, const uint32_t averageQuerySize, const uint32_t maxRegionsFactor)
{
  gpu_buffer_t* const mBuff                  = (gpu_buffer_t *) fmiBuffer;
  const size_t        sizeBuff               = mBuff->sizeBuffer * 0.90;
  const uint32_t      averageRegionsPerQuery = GPU_MAX(GPU_DIV_CEIL(averageQuerySize, maxRegionsFactor), GPU_FMI_MIN_REGIONS);
  const size_t        bytesPerQuery          = gpu_fmi_asearch_size_per_query(averageQuerySize, averageRegionsPerQuery);
  const uint32_t      numQueries             = sizeBuff / bytesPerQuery;
  //set the type of the buffer
  mBuff->typeBuffer = GPU_FMI_ADAPT_SEARCH;
  // Set real size of the input
  mBuff->data.asearch.numMaxQueries         = numQueries;
  mBuff->data.asearch.numMaxBases           = numQueries * averageQuerySize;
  mBuff->data.asearch.numMaxRegions         = numQueries * averageRegionsPerQuery;
  // Internal data information
  mBuff->data.asearch.maxRegionsFactor      = maxRegionsFactor;
  // Set the corresponding buffer layout
  gpu_fmi_asearch_reallocate_host_buffer_layout(mBuff);
  gpu_fmi_asearch_reallocate_device_buffer_layout(mBuff);
}

void gpu_fmi_asearch_init_and_realloc_buffer_(void* const fmiBuffer, const uint32_t maxRegionsFactor, const uint32_t totalBases,
                                              const uint32_t totalQueries, const uint32_t totalRegions)
{
  // Buffer reinitialization
  gpu_buffer_t* const mBuff                  = (gpu_buffer_t *) fmiBuffer;
  const uint32_t      averageQuerySize       = GPU_DIV_CEIL(totalBases, totalQueries);
  const uint32_t      averageRegionsPerQuery = GPU_MAX(GPU_DIV_CEIL(averageQuerySize, maxRegionsFactor), GPU_FMI_MIN_REGIONS);
  // Remap the buffer layout with new information trying to fit better
  gpu_fmi_asearch_init_buffer_(fmiBuffer, averageQuerySize, maxRegionsFactor);
  // Checking if we need to reallocate a bigger buffer
  if( (totalBases   > gpu_fmi_asearch_buffer_get_max_bases_(fmiBuffer))   ||
      (totalQueries > gpu_fmi_asearch_buffer_get_max_queries_(fmiBuffer)) ||
      (totalRegions > gpu_fmi_asearch_buffer_get_max_regions_(fmiBuffer))){
    // Resize the GPU buffer to fit the required input
    const uint32_t      idSupDevice             = mBuff->idSupportedDevice;
    const float         resizeFactor            = 2.0;
    const size_t        bytesPerSearchBuffer    = totalQueries * gpu_fmi_asearch_size_per_query(averageQuerySize, averageRegionsPerQuery);
    //Recalculate the minimum buffer size
    //printf("RESIZE[FMI_SEARCH] %d %d \n",  mBuff->sizeBuffer, bytesPerSearchBuffer * resizeFactor);
    mBuff->sizeBuffer = bytesPerSearchBuffer * resizeFactor;
    //FREE HOST AND DEVICE BUFFER
    GPU_ERROR(gpu_buffer_free(mBuff));
    //Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
    //ALLOCATE HOST AND DEVICE BUFFER
    CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
    CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));
    // Remap the buffer layout with the new size
    gpu_fmi_asearch_init_buffer_(fmiBuffer, averageQuerySize, maxRegionsFactor);
  }
}

/************************************************************
Functions to transfer data HOST <-> DEVICE (E. SEARCH)
************************************************************/

gpu_error_t gpu_fmi_asearch_transfer_CPU_to_GPU(gpu_buffer_t *mBuff)
{
  const gpu_fmi_asearch_queries_buffer_t* qryBuff  = &mBuff->data.asearch.queries;
  const gpu_fmi_asearch_regions_buffer_t* regBuff  = &mBuff->data.asearch.regions;
  const cudaStream_t                      idStream =  mBuff->listStreams[mBuff->idStream];
  size_t                                  cpySize  =  0;
  float                                   bufferUtilization;
  // Defining buffer offsets
  cpySize += qryBuff->numBases   * sizeof(gpu_fmi_search_query_t);
  cpySize += qryBuff->numQueries * sizeof(gpu_fmi_search_query_info_t);
  cpySize += qryBuff->numQueries * sizeof(gpu_fmi_search_region_t);
  cpySize += regBuff->numRegions * sizeof(gpu_sa_search_inter_t);
  cpySize += regBuff->numRegions * sizeof(gpu_fmi_search_region_info_t);
  bufferUtilization = (double)cpySize / (double)mBuff->sizeBuffer;
  // Compacting tranferences with high buffer occupation
  if(bufferUtilization > 0.15){
    cpySize  = ((void *) (qryBuff->d_regions + qryBuff->numQueries)) - ((void *) qryBuff->d_queries);
    CUDA_ERROR(cudaMemcpyAsync(qryBuff->d_queries, qryBuff->h_queries, cpySize, cudaMemcpyHostToDevice, idStream));
  }else{
    //Transfer Queries to GPU
    cpySize = qryBuff->numBases * sizeof(gpu_fmi_search_query_t);
    CUDA_ERROR(cudaMemcpyAsync(qryBuff->d_queries, qryBuff->h_queries, cpySize, cudaMemcpyHostToDevice, idStream));
    //Transfer to GPU the information associated with Queries
    cpySize = qryBuff->numQueries * sizeof(gpu_fmi_search_query_info_t);
    CUDA_ERROR(cudaMemcpyAsync(qryBuff->d_queryInfo, qryBuff->h_queryInfo, cpySize, cudaMemcpyHostToDevice, idStream));
    //Transfer Candidates to GPU
    cpySize = qryBuff->numQueries * sizeof(gpu_fmi_search_region_t);
    CUDA_ERROR(cudaMemcpyAsync(qryBuff->d_regions, qryBuff->h_regions, cpySize, cudaMemcpyHostToDevice, idStream));
  }
  // Suceed
  return (SUCCESS);
}

gpu_error_t gpu_fmi_asearch_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff)
{
  const gpu_fmi_asearch_queries_buffer_t* qryBuff   = &mBuff->data.asearch.queries;
  const gpu_fmi_asearch_regions_buffer_t* regBuff   = &mBuff->data.asearch.regions;
  const cudaStream_t                      idStream  =  mBuff->listStreams[mBuff->idStream];
        size_t                            cpySize   =  0;
        float                             bufferUtilization;
  // Defining buffer offsets
  cpySize += qryBuff->numBases   * sizeof(gpu_fmi_search_query_t);
  cpySize += qryBuff->numQueries * sizeof(gpu_fmi_search_query_info_t);
  cpySize += qryBuff->numQueries * sizeof(gpu_fmi_search_region_t);
  cpySize += regBuff->numRegions * sizeof(gpu_sa_search_inter_t);
  cpySize += regBuff->numRegions * sizeof(gpu_fmi_search_region_info_t);
  bufferUtilization = (double)cpySize / (double)mBuff->sizeBuffer;
  // Compacting tranferences with high buffer occupation
  if(bufferUtilization > 0.15){
    cpySize  = ((void *) (regBuff->h_regionsOffsets + regBuff->numRegions)) - ((void *) qryBuff->h_regions);
    CUDA_ERROR(cudaMemcpyAsync(qryBuff->h_regions, qryBuff->d_regions, cpySize, cudaMemcpyDeviceToHost, idStream));
  }else{
    // Transfer Candidates to GPU
    cpySize = qryBuff->numQueries * sizeof(gpu_fmi_search_region_t);
    CUDA_ERROR(cudaMemcpyAsync(qryBuff->h_regions, qryBuff->d_regions, cpySize, cudaMemcpyDeviceToHost, idStream));
    // Transfer Candidates to GPU
    cpySize = regBuff->numRegions * sizeof(gpu_sa_search_inter_t);
    CUDA_ERROR(cudaMemcpyAsync(regBuff->h_intervals, regBuff->d_intervals, cpySize, cudaMemcpyDeviceToHost, idStream));
    // Transfer Candidates to GPU
    cpySize = regBuff->numRegions * sizeof(gpu_fmi_search_region_info_t);
    CUDA_ERROR(cudaMemcpyAsync(regBuff->h_regionsOffsets, regBuff->d_regionsOffsets, cpySize, cudaMemcpyDeviceToHost, idStream));
  }
  // Suceed
  return (SUCCESS);
}

void gpu_fmi_asearch_send_buffer_(void* const fmiBuffer, const uint32_t numQueries, const uint32_t numBases, const uint32_t numRegions,
                                  const uint32_t occMinThreshold, const uint32_t extraSteps, const uint32_t alphabetSize)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;
  //Set real size of the input
  mBuff->data.asearch.extraSteps         = extraSteps;
  mBuff->data.asearch.occShrinkFactor    = GPU_BASE2LOG(alphabetSize);
  mBuff->data.asearch.occMinThreshold    = occMinThreshold;
  mBuff->data.asearch.queries.numQueries = numQueries;
  mBuff->data.asearch.queries.numBases   = numBases;
  mBuff->data.asearch.regions.numRegions = numRegions;
  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
  GPU_ERROR(gpu_fmi_asearch_transfer_CPU_to_GPU(mBuff));
  GPU_ERROR(gpu_fmi_asearch_process_buffer(mBuff));
  GPU_ERROR(gpu_fmi_asearch_transfer_GPU_to_CPU(mBuff));
}

void gpu_fmi_asearch_receive_buffer_(const void* const fmiBuffer)
{
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  const cudaStream_t  idStream    =  mBuff->listStreams[mBuff->idStream];
  //Synchronize Stream (the thread wait for the commands done in the stream)
  CUDA_ERROR(cudaStreamSynchronize(idStream));
  #ifdef GPU_FMI_DEBUG
    gpu_buffer_fmi_asearch_process_histogram(mBuff);
  #endif
}

#ifdef GPU_FMI_DEBUG
void gpu_fmi_asearch_print_histograms()
{
  uint32_t i = 0;
  printf("histogram_queries: Setting the num reads \n");
  for(i = 0; i < 15; i++)
    printf("%u\t%u\n", i, histogram_queries[i]);
  printf("\n\n");

  printf("histogram_query_candidates: Setting num candidates per read \n");
  for(i = 0; i < 15; i++)
    printf("%u\t%u\n", i, histogram_query_candidates[i]);
  printf("\n\n");

  printf("histogram_regions: Setting num region per size \n");
  for(i = 0; i < 110; i++)
    printf("%u\t%u\n", i, histogram_regions[i]);
  printf("\n\n");

  printf("histogram_region_candidates: Setting num candidates per size \n");
  for(i = 0; i < 110; i++)
    printf("%u\t%u\n", i, histogram_region_candidates[i]);
  printf("\n\n");

  printf("histogram_coverage: Setting all regions coverage over the read \n");
  for(i = 0; i < 110; i++)
    printf("%u\t%u\n", i, histogram_coverage[i]);
  printf("\n\n");
}
#endif

#endif /* GPU_FMI_PRIMITIVES_ASEARCH_C_ */

