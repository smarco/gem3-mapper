/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_BPM_ALIGN_PRIMITIVES_H_
#define GPU_BPM_ALIGN_PRIMITIVES_H_

#include "gpu_commons.h"
#include "gpu_scheduler.h"

/********************************
Common constants for Device & Host
*********************************/

/* Defines to distribute and balance the work in GPU*/
#define GPU_BPM_ALIGN_PEQ_LENGTH_PER_CUDA_THREAD  128
#define GPU_BPM_ALIGN_NUM_BUCKETS_FOR_BINNING     (GPU_WARP_SIZE + 1)

#define GPU_BPM_ALIGN_CANDIDATES_BUFFER_PADDING   10
#define GPU_BPM_ALIGN_MIN_ELEMENTS                150  // MIN elements per buffer (related to the SM -2048th-)
#define GPU_BPM_ALIGN_MAX_SIZE_CANDIDATE          750

/*****************************
Internal Objects
*****************************/

typedef struct {
  uint32_t            	       numCigars;
  uint32_t            	       numCigarEntries;
  gpu_bpm_align_cigar_entry_t *h_cigars;
  gpu_bpm_align_cigar_entry_t *d_cigars;
  gpu_bpm_align_cigar_info_t  *h_cigarsInfo;
  gpu_bpm_align_cigar_info_t  *d_cigarsInfo;
} gpu_bpm_align_cigars_buffer_t;

typedef struct {
  uint32_t                    numCandidates;
  gpu_bpm_align_cand_info_t  *h_candidatesInfo;
  gpu_bpm_align_cand_info_t  *d_candidatesInfo;
} gpu_bpm_align_candidates_buffer_t;

typedef struct {
  uint32_t              	     totalQueriesPEQs;
  uint32_t              	     totalQueriesBases;
  uint32_t              	     numQueries;
  gpu_bpm_align_qry_entry_t   *h_queries;
  gpu_bpm_align_qry_entry_t   *d_queries;
  gpu_bpm_align_peq_entry_t   *h_peq;
  gpu_bpm_align_peq_entry_t   *d_peq;
  gpu_bpm_align_qry_info_t    *h_qinfo;
  gpu_bpm_align_qry_info_t    *d_qinfo;
} gpu_bpm_align_queries_buffer_t;

/*************************************
Specific types for the Devices (GPUs)
**************************************/

typedef struct {
  uint4 bitmap[GPU_BPM_ALIGN_PEQ_ALPHABET_SIZE];
} gpu_bpm_align_device_qry_entry_t;

typedef struct {
  char vCoord;
  char hCoord;
  char cigarEvent;
  char matchLenght;
} gpu_bpm_align_device_cigar_entry_t;

/*****************************
General Object
*****************************/

typedef struct {
  uint32_t                    		  maxPEQEntries;
  uint32_t                    		  maxQueryBases;
  uint32_t                    		  maxCandidates;
  uint32_t                    		  maxQueries;
  uint32_t                    		  maxCandidateSize;
  uint32_t                    		  maxCigars;
  uint32_t                    		  maxCigarEntries;
  uint32_t                    		  maxReorderBuffer;
  uint32_t                    		  maxBuckets;
  uint32_t                    		  queryBinSize;
  bool                        		  queryBinning;
  gpu_bpm_align_queries_buffer_t    queries;
  gpu_bpm_align_candidates_buffer_t candidates;
  gpu_scheduler_buffer_t	          reorderBuffer;
  gpu_bpm_align_cigars_buffer_t 	  cigars;
} gpu_bpm_align_buffer_t;

#include "gpu_buffer.h"

/* Functions to initialize all the BPM resources */
float       gpu_bpm_align_size_per_candidate(const uint32_t averageQuerySize, const uint32_t candidatesPerQuery);
uint32_t    gpu_bpm_align_candidates_for_binning_padding();
void        gpu_bpm_align_reallocate_host_buffer_layout(gpu_buffer_t* const mBuff);
void        gpu_bpm_align_reallocate_device_buffer_layout(gpu_buffer_t* const mBuff);
/* Functions to send & process a BPM buffer to GPU */
gpu_error_t gpu_bpm_align_reordering_buffer(gpu_buffer_t* const mBuff);
gpu_error_t gpu_bpm_align_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff);
gpu_error_t gpu_bpm_align_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff);
/* Functions to receive & process a BPM buffer from GPU */
gpu_error_t gpu_bpm_align_reordering_alignments(gpu_buffer_t* const mBuff);
gpu_error_t gpu_bpm_align_reorder_process(const gpu_bpm_align_queries_buffer_t* const qry, const gpu_bpm_align_candidates_buffer_t* const cand,
                                          gpu_scheduler_buffer_t* const rebuff, gpu_bpm_align_cigars_buffer_t* const res);
/* DEVICE Kernels */
gpu_error_t gpu_bpm_align_process_buffer(gpu_buffer_t *mBuff);

#endif /* GPU_BPM_ALIGN_PRIMITIVES_H_ */



