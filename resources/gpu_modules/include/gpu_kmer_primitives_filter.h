/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_KMER_PRIMITIVES_H_
#define GPU_KMER_PRIMITIVES_H_

#include "gpu_commons.h"

/********************************
Common constants for Device & Host
*********************************/

#define GPU_KMER_FILTER_MIN_ELEMENTS       2048  // MIN elements per buffer (related to the SM -2048th-)

/*****************************
Internal Objects
*****************************/

typedef struct {
  uint32_t             numBases;
  uint32_t             numQueries;
  gpu_kmer_filter_qry_entry_t *d_queries;
  gpu_kmer_filter_qry_entry_t *h_queries;
  gpu_kmer_filter_qry_info_t  *d_queryInfo;
  gpu_kmer_filter_qry_info_t  *h_queryInfo;
} gpu_kmer_filter_queries_buffer_t;

typedef struct {
  uint32_t             numAlignments;
  gpu_kmer_filter_alg_entry_t *h_alignments;
  gpu_kmer_filter_alg_entry_t *d_alignments;
} gpu_kmer_filter_alignments_buffer_t;

typedef struct {
  uint32_t             numCandidates;
  gpu_kmer_filter_cand_info_t *h_candidates;
  gpu_kmer_filter_cand_info_t *d_candidates;
} gpu_kmer_filter_candidates_buffer_t;


/*****************************
General Object
*****************************/

typedef struct {
  uint32_t                     maxBases;
  uint32_t                     maxCandidates;
  uint32_t                     maxQueries;
  uint32_t                     maxAlignments;
  uint32_t                     maxError;
  gpu_kmer_filter_queries_buffer_t    queries;
  gpu_kmer_filter_candidates_buffer_t candidates;
  gpu_kmer_filter_alignments_buffer_t alignments;
} gpu_kmer_filter_buffer_t;

#include "gpu_buffer.h"

/* Functions to initialize all the KMER resources */
float       gpu_kmer_filter_size_per_candidate(const uint32_t averageQuerySize, const uint32_t candidatesPerQuery);
void        gpu_kmer_filter_reallocate_host_buffer_layout(gpu_buffer_t* mBuff);
void        gpu_kmer_filter_reallocate_device_buffer_layout(gpu_buffer_t* mBuff);
/* Functions to send & process a KMER buffer to GPU */
gpu_error_t gpu_kmer_filter_transfer_CPU_to_GPU(gpu_buffer_t *mBuff);
gpu_error_t gpu_kmer_filter_transfer_GPU_to_CPU(gpu_buffer_t *mBuff);
/* DEVICE Kernels */
gpu_error_t gpu_kmer_filter_process_buffer(gpu_buffer_t *mBuff);

#endif /* GPU_KMER_PRIMITIVES_H_ */



