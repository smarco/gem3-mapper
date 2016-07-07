/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include "gpu_index_modules.h"
#include "gpu_commons.h"
#include "gpu_fmi_structure.h"
#include "gpu_sa_index.h"

#ifndef GPU_FMI_PRIMITIVES_ASEARCH_H_
#define GPU_FMI_PRIMITIVES_ASEARCH_H_

/********************************
Common constants for Device & Host
*********************************/

/* Defines related to FMI adaptative BACKWARD-SEARCH primitives */
#define GPU_FMI_ENTRIES_PER_QUERY      2
#define GPU_FMI_BASE_QUERY_LENGTH      8 //Bits per base query
#define GPU_FMI_BASES_PER_QUERY_ENTRY  (GPU_UINT64_LENGTH / GPU_FMI_BASE_QUERY_LENGTH) //Bases per internal entry

#define GPU_FMI_THREADS_PER_QUERY (GPU_FMI_THREADS_PER_ENTRY * GPU_FMI_ENTRIES_PER_SEED)
#define GPU_FMI_D                 2
#define GPU_FMI_STEPS             4
#define GPU_FMI_OCC_THRESHOLD     20  //Min number of regions per seed
#define GPU_FMI_RATIO_REGIONS     10  //Extract 10 seeds (covering a 9% query error)
#define GPU_FMI_MIN_REGIONS       1

/*****************************
Internal Objects (Adaptative  Search)
*****************************/

typedef struct {
  uint32_t                     numBases;
  uint32_t                     numQueries;
  gpu_fmi_search_query_t       *d_queries;
  gpu_fmi_search_query_t       *h_queries;
  gpu_fmi_search_query_info_t  *d_queryInfo;
  gpu_fmi_search_query_info_t  *h_queryInfo;
  gpu_fmi_search_region_t      *h_regions;
  gpu_fmi_search_region_t      *d_regions;
} gpu_fmi_asearch_queries_buffer_t;

typedef struct {
  uint32_t                      numRegions;
  gpu_sa_search_inter_t         *h_intervals;
  gpu_sa_search_inter_t         *d_intervals;
  gpu_fmi_search_region_info_t  *h_regionsOffsets;
  gpu_fmi_search_region_info_t  *d_regionsOffsets;
} gpu_fmi_asearch_regions_buffer_t;


/*****************************
Internal Objects (General)
*****************************/
typedef struct {
  uint32_t                          numMaxBases;
  uint32_t                          numMaxQueries;
  uint32_t                          numMaxRegions;
  uint32_t                          maxRegionsFactor;
  uint32_t                          occMinThreshold;
  uint32_t                          extraSteps;
  uint32_t                          occShrinkFactor;
  gpu_fmi_asearch_queries_buffer_t  queries;
  gpu_fmi_asearch_regions_buffer_t  regions;
} gpu_fmi_asearch_buffer_t;

#include "gpu_buffer.h"

/* Functions to init the buffers (E. SEARCH) */
size_t      gpu_fmi_asearch_size_per_query(const uint32_t averageQuerySize, const uint32_t averageRegionsPerQuery);
void        gpu_fmi_asearch_reallocate_host_buffer_layout(gpu_buffer_t* mBuff);
void        gpu_fmi_asearch_reallocate_device_buffer_layout(gpu_buffer_t* mBuff);

/* Functions to transfer data HOST <-> DEVICE (E. SEARCH) */
gpu_error_t gpu_fmi_asearch_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff);
gpu_error_t gpu_fmi_asearch_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff);

/* DEVICE Kernels */
gpu_error_t gpu_fmi_asearch_process_buffer(gpu_buffer_t* const mBuff);


#endif /* GPU_FMI_PRIMITIVES_ASEARCH_H_ */

