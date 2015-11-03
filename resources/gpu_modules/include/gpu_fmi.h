/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Common headers and data structures for BPM on GPU library
 */

#ifndef GPU_FMI_H_
#define GPU_FMI_H_

#include "gpu_commons.h"
#include "gpu_scheduler.h"
#include "gpu_index.h"

/********************************
Common constants for Device & Host
*********************************/

/* Defines related to FMI BACKWARD-SEARCH primitives */
#define GPU_FMI_ENTRIES_PER_SEED			2
#define GPU_FMI_SEED_ENTRY_LENGTH			128															// 128 bits
#define GPU_FMI_SEED_FIELD_SIZE				8															// 8 bits
#define GPU_FMI_SEED_MAX_CHARS				60															// 60 bases
#define GPU_FMI_SEED_CHAR_LENGTH			2 		    												// 2 bits

#define GPU_FMI_SEED_THREADS_PER_ENTRY		(GPU_FMI_THREADS_PER_ENTRY * GPU_FMI_ENTRIES_PER_SEED)
#define GPU_FMI_SEED_ENTRIES_PER_WARP		(GPU_WARP_SIZE / GPU_FMI_SEED_THREADS_PER_ENTRY)
#define GPU_FMI_SEED_BASES_PER_ENTRY		(GPU_UINT64_LENGTH / GPU_FMI_SEED_CHAR_LENGTH)
#define GPU_FMI_SEARCH_SEEDS_BUFFER_PADDING 10
#define GPU_FMI_SEARCH_MIN_ELEMENTS		 	(2048 / GPU_FMI_SEED_THREADS_PER_ENTRY)						// Min elements per buffer (related to the SM -2048th-)

/* Defines related to FMI DECODE primitives */
#define GPU_FMI_ENTRIES_PER_DECODE		 	2 											       			// Necessary FMI entries for unaligned FMI threads
#define GPU_FMI_THREADS_PER_COUNTERS	 	2 											       			// Necessary threads to load all the counters

#define GPU_FMI_DECODE_THREADS_PER_ENTRY 	(GPU_FMI_THREADS_PER_ENTRY * GPU_FMI_ENTRIES_PER_DECODE)  	// FMI_THREADS_PER_ENTRY + 1 (rounded to pow of 2)
#define GPU_FMI_DECODE_THREADS_PER_LOAD	 	(GPU_FMI_THREADS_PER_ENTRY + GPU_FMI_THREADS_PER_COUNTERS)  // Threads to load an entry for the decode primitive
#define GPU_FMI_DECODE_POS_BUFFER_PADDING 	10
#define GPU_FMI_DECODE_MIN_ELEMENTS		 	(2048 / GPU_FMI_DECODE_THREADS_PER_ENTRY) 					// Min elements per buffer (related to the SM)

/*****************************
Internal Objects (General)
*****************************/

typedef struct {
	uint32_t						 numMaxSeeds;
	uint32_t						 numMaxIntervals;
	gpu_fmi_search_seeds_buffer_t 	 seeds;
	gpu_fmi_search_sa_inter_buffer_t saIntervals;
} gpu_fmi_search_buffer_t;

typedef struct {
	uint32_t						 numMaxInitPositions;
	uint32_t						 numMaxEndPositions;
	uint32_t						 samplingRate;
	gpu_fmi_decode_init_pos_buffer_t initPositions;
	gpu_fmi_decode_end_pos_buffer_t  endPositions;
} gpu_fmi_decode_buffer_t;

/*****************************
Internal Objects (Search)
*****************************/

typedef struct {
	uint32_t					numSeeds;
	gpu_fmi_search_seed_t		*d_seeds;
	gpu_fmi_search_seed_t		*h_seeds;
} gpu_fmi_search_seeds_buffer_t;

typedef struct {
	uint32_t					numIntervals;
	gpu_fmi_search_sa_inter_t 	*h_intervals;
	gpu_fmi_search_sa_inter_t	*d_intervals;
} gpu_fmi_search_sa_inter_buffer_t;


/*****************************
Internal Objects (Decode)
*****************************/

typedef struct {
	uint32_t					numDecodings;
	gpu_fmi_decode_init_pos_t	*h_initBWTPos;
	gpu_fmi_decode_init_pos_t	*d_initBWTPos;
} gpu_fmi_decode_init_pos_buffer_t;

typedef struct {
	uint32_t					numDecodings;
	gpu_fmi_decode_end_pos_t	*h_endBWTPos;
	gpu_fmi_decode_end_pos_t	*d_endBWTPos;
} gpu_fmi_decode_end_pos_buffer_t;


GPU_INLINE size_t 	gpu_fmi_search_input_size();
GPU_INLINE void 	gpu_fmi_search_reallocate_host_buffer_layout(gpu_buffer_t* mBuff);
GPU_INLINE void 	gpu_fmi_search_reallocate_device_buffer_layout(gpu_buffer_t* mBuff);
GPU_INLINE void 	gpu_fmi_search_init_buffer_(void* fmiBuffer);


#endif /* GPU_FMI_H_ */



