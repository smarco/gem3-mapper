/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Common headers and data structures for BPM on GPU library
 */

#include "gpu_commons.h"

#ifndef GPU_BPM_H_
#define GPU_BPM_H_

/********************************
Common constants for Device & Host
*********************************/

/* Defines to distribute and balance the work in GPU*/
#define	GPU_BPM_PEQ_LENGTH_PER_CUDA_THREAD		128
#define	GPU_BPM_NUM_BUCKETS_FOR_BINNING			(GPU_WARP_SIZE + 1)

#define GPU_BPM_CANDIDATES_BUFFER_PADDING		10
#define GPU_BPM_MIN_ELEMENTS		 			2048	// MIN elements per buffer (related to the SM -2048th-)


/*****************************
Internal Objects
*****************************/

typedef struct {
	uint32_t			numAlignments;
 	uint32_t			numReorderedAlignments;
 	gpu_bpm_alg_entry_t	*h_alignments;
 	gpu_bpm_alg_entry_t	*d_alignments;
 	gpu_bpm_alg_entry_t	*h_reorderAlignments;
 	gpu_bpm_alg_entry_t	*d_reorderAlignments;
} gpu_bpm_alignments_buffer_t;

typedef struct {
	uint32_t	numBuckets;
	uint32_t	candidatesPerBuffer;
	uint32_t	numWarps;
	uint32_t   	*h_reorderBuffer;
	uint32_t   	*d_reorderBuffer;
	uint32_t 	*h_initPosPerBucket;
	uint32_t 	*h_initWarpPerBucket;
	uint32_t 	*d_initPosPerBucket;
	uint32_t 	*d_initWarpPerBucket;
} gpu_bpm_reorder_buffer_t;

typedef struct {
	uint32_t 			numCandidates;
	gpu_bpm_cand_info_t *h_candidates;
	gpu_bpm_cand_info_t *d_candidates;
} gpu_bpm_candidates_buffer_t;


/*************************************
Specific types for the Devices (GPUs)
**************************************/

typedef struct {
	uint4 bitmap[GPU_BPM_PEQ_ALPHABET_SIZE];
} gpu_bpm_device_qry_entry_t;


/*************************************
Specific types for the Host (CPU)
**************************************/

typedef struct {
	uint32_t 				totalQueriesEntries;
	uint32_t 				numQueries;
	gpu_bpm_qry_entry_t 	*h_queries;
	gpu_bpm_qry_entry_t 	*d_queries;
	gpu_bpm_qry_info_t 		*h_qinfo;
	gpu_bpm_qry_info_t 		*d_qinfo;
} gpu_bpm_queries_buffer_t;


/*****************************
General Object
*****************************/

typedef struct {
	uint32_t 					maxPEQEntries;
	uint32_t					maxCandidates;
	uint32_t					maxQueries;
	uint32_t					maxReorderBuffer;
	uint32_t					maxAlignments;
	uint32_t					maxBuckets;
	gpu_bpm_queries_buffer_t 	queries;
	gpu_bpm_candidates_buffer_t candidates;
	gpu_bpm_reorder_buffer_t 	reorderBuffer;
	gpu_bpm_alignments_buffer_t alignments;
} gpu_bpm_buffer_t;

#include "gpu_buffers.h"

/* Functions to initialize all the BPM resources */
float 		gpu_bpm_size_per_candidate(const uint32_t averageNumPEQEntries, const uint32_t candidatesPerQuery);
uint32_t 	gpu_bpm_candidates_for_binning_padding();
void 		gpu_bpm_reallocate_host_buffer_layout(gpu_buffer_t* mBuff);
void 		gpu_bpm_reallocate_device_buffer_layout(gpu_buffer_t* mBuff);
/* Functions to send & process a BPM buffer to GPU */
gpu_error_t gpu_bpm_reordering_buffer(gpu_buffer_t *mBuff);
gpu_error_t gpu_bpm_transfer_CPU_to_GPU(gpu_buffer_t *mBuff);
gpu_error_t gpu_bpm_transfer_GPU_to_CPU(gpu_buffer_t *mBuff);
/* Functions to receive & process a BPM buffer from GPU */
gpu_error_t gpu_bpm_reordering_alignments(gpu_buffer_t *mBuff);
/* DEVICE Kernels */
gpu_error_t gpu_bpm_process_buffer(gpu_buffer_t *mBuff);

#endif /* GPU_BPM_H_ */



