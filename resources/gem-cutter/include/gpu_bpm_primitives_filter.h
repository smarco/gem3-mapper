/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_BPM_PRIMITIVES_FILTER_H_
#define GPU_BPM_PRIMITIVES_FILTER_H_

#include "gpu_commons.h"
#include "gpu_scheduler.h"

/********************************
Common constants for Device & Host
*********************************/

/* Defines to distribute and balance the work in GPU*/
#define GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD  128
#define GPU_BPM_FILTER_NUM_BUCKETS_FOR_BINNING     (GPU_WARP_SIZE + 1)

#define GPU_BPM_FILTER_CANDIDATES_BUFFER_PADDING   10
#define GPU_BPM_FILTER_MIN_ELEMENTS                150  // MIN elements per buffer (related to the SM -2048th-)

#define GPU_BPM_FILTER_SCORE_INF                   (GPU_UINT32_ONES)
#define GPU_BPM_FILTER_CUTOFF_DISABLED_KEY         (GPU_UINT32_ONES)
#define GPU_BPM_FILTER_CUTOFF_ACTIVE     	         true

#define GPU_BPM_FILTER_PENDING_NUM_TASK_LIST       2

#define GPU_BPM_FILTER_CUTOFF_MAX_TILE_LENGHT      512
#define GPU_BPM_FILTER_CUTOFF_MAX_READ_LENGHT      65536
#define GPU_BPM_FILTER_CUTOFF_MAX_HIST_SIZE        (GPU_BPM_FILTER_CUTOFF_MAX_READ_LENGHT / GPU_BPM_FILTER_CUTOFF_MAX_TILE_LENGHT)
#define GPU_BPM_FILTER_CUTOFF_MAX_TILE_BIN         (GPU_BPM_FILTER_CUTOFF_MAX_TILE_LENGHT / GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD)

#define GPU_BPM_FILTER_MIN_THREADS_PER_LAUNCH      128 // related to the SM -2048th-
#define GPU_BPM_FILTER_THREADS_PER_TILE            (GPU_BPM_FILTER_CUTOFF_MAX_TILE_LENGHT / GPU_BPM_FILTER_PEQ_LENGTH_PER_CUDA_THREAD)
#define GPU_BPM_FILTER_MIN_TILES_PER_LAUNCH        (GPU_BPM_FILTER_MIN_THREADS_PER_LAUNCH / GPU_BPM_FILTER_THREADS_PER_TILE)

/*****************************
Internal Object
*****************************/

typedef struct {
  uint32_t            		  numAlignments;
  uint32_t            		  numReorderedAlignments;
  gpu_bpm_filter_alg_entry_t *h_alignments;
  gpu_bpm_filter_alg_entry_t *d_alignments;
  gpu_bpm_filter_alg_entry_t *h_reorderAlignments;
  gpu_bpm_filter_alg_entry_t *d_reorderAlignments;
} gpu_bpm_filter_alignments_buffer_t;


typedef struct {
  uint32_t             		   numCandidates;
  uint32_t             		   numCandidatesEntries;
  gpu_bpm_filter_cand_info_t  *h_candidatesInfo;
  gpu_bpm_filter_cand_info_t  *d_candidatesInfo;
  gpu_bpm_filter_cand_info_t  *h_candidates;
  gpu_bpm_filter_cand_info_t  *d_candidates;
} gpu_bpm_filter_candidates_buffer_t;


typedef struct {
  uint32_t              	    totalQueriesEntries;
  uint32_t              	    numQueries;
  gpu_bpm_filter_qry_entry_t   *h_queries;
  gpu_bpm_filter_qry_entry_t   *d_queries;
  gpu_bpm_filter_qry_info_t    *h_qinfo;
  gpu_bpm_filter_qry_info_t    *d_qinfo;
} gpu_bpm_filter_queries_buffer_t;

/*************************************
Specific types for the Devices (GPUs)
**************************************/

typedef struct {
  uint4 bitmap[GPU_BPM_FILTER_PEQ_ALPHABET_SIZE];
} gpu_bpm_filter_device_qry_entry_t;


/*************************************
Specific temporal types for the Host (CPU)
**************************************/

typedef uint32_t gpu_bpm_filter_cand_error_t;

typedef struct {
  uint32_t                     numErrorEntries;
  bool                         pendingTasks;
  uint32_t                     numCurrPendingTasks;
  uint32_t                     numNextPendingTasks;
   int32_t                     maxTileDepth;
   int32_t                     minTileDepth;
  uint32_t                     tileDepthHist[GPU_BPM_FILTER_CUTOFF_MAX_HIST_SIZE];
  gpu_bpm_filter_cand_error_t  *h_error;
  uint32_t                     *h_currPendingTasks;
  uint32_t                     *h_nextPendingTasks;
  uint32_t                     *h_pendingTasks_tmpSpace[GPU_BPM_FILTER_PENDING_NUM_TASK_LIST];
} gpu_bpm_filter_cutoff_buffer_t;

/*****************************
General Object
*****************************/

typedef struct {
  uint32_t                    		   maxPEQEntries;
  uint32_t                    		   maxCandidates;
  uint32_t                    		   maxQueries;
  uint32_t                    		   maxReorderBuffer;
  uint32_t                    		   maxAlignments;
  uint32_t                    		   maxBuckets;
  uint32_t                    		   maxCutOffEntries;
  uint32_t                             maxPendingTasks;
  uint32_t                    		   queryBinSize;
  uint32_t                             maxQuerySize;
  bool                        		   queryBinning;
  bool                                 activeCutOff;
  gpu_bpm_filter_queries_buffer_t      queries;
  gpu_bpm_filter_candidates_buffer_t   candidates;
  gpu_scheduler_buffer_t               reorderBuffer;
  gpu_bpm_filter_alignments_buffer_t   alignments;
  gpu_bpm_filter_cutoff_buffer_t       cutoff;
} gpu_bpm_filter_buffer_t;

#include "gpu_buffer.h"

/* Functions to initialize all the BPM resources */
float       gpu_bpm_filter_size_per_candidate(const uint32_t averageQuerySize, const uint32_t candidatesPerQuery);
uint32_t    gpu_bpm_filter_candidates_for_binning_padding();
void        gpu_bpm_filter_reallocate_host_buffer_layout(gpu_buffer_t* const mBuff);
void        gpu_bpm_filter_reallocate_device_buffer_layout(gpu_buffer_t* const mBuff);
/* Functions to send & process a BPM buffer to GPU */
gpu_error_t gpu_bpm_filter_reordering_buffer(gpu_buffer_t* const mBuff);
gpu_error_t gpu_bpm_filter_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff);
gpu_error_t gpu_bpm_filter_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff);
/* Functions to send & process a BPM internal buffer to GPU */
gpu_error_t gpu_bpm_filter_intermediate_data_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff);
gpu_error_t gpu_bpm_filter_intermediate_data_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff);
/* Functions to receive & process a BPM buffer from GPU */
gpu_error_t gpu_bpm_filter_reordering_alignments(gpu_buffer_t* const mBuff);
gpu_error_t gpu_bpm_filter_reorder_process(const gpu_bpm_filter_queries_buffer_t* const qry, const gpu_bpm_filter_candidates_buffer_t* const cand,
                                    	   gpu_scheduler_buffer_t* const rebuff, gpu_bpm_filter_alignments_buffer_t* const res, const uint32_t idKey);
/* DEVICE Kernels */
gpu_error_t gpu_bpm_filter_process_buffer(gpu_buffer_t* const mBuff);
/* Cutoff primitives */
gpu_error_t gpu_bpm_filter_device_synch(gpu_buffer_t* const mBuff);
gpu_error_t gpu_bpm_filter_reordering_alignments_cutoff(gpu_buffer_t* const mBuff);

#endif /* GPU_BPM_PRIMITIVES_FILTER_H_ */



