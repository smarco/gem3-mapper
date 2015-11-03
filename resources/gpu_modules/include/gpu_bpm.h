/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Common headers and data structures for BPM on GPU library
 */

#ifndef GPU_BPM_H_
#define GPU_BPM_H_

#include "gpu_commons.h"
#include "gpu_scheduler.h"

/********************************
Common constants for Device & Host
*********************************/

/* Defines to distribute and balance the work in GPU*/
#define	GPU_BPM_PEQ_LENGTH_PER_CUDA_THREAD		128
#define	GPU_BPM_NUM_BUCKETS_FOR_BINNING			(GPU_WARP_SIZE + 1)

#define GPU_BPM_CANDIDATES_BUFFER_PADDING		10
#define GPU_BPM_MIN_ELEMENTS		 			2048	// Min elements per buffer (related to the SM -2048th-)


/*****************************
Internal Objects
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

#endif /* GPU_BPM_H_ */



