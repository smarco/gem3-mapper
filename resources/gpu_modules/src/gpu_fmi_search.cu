/*
 * PROJECT: Thread-cooperative FM-index on GPU
 * FILE: genIndex.h
 * DATE: 1/9/2015
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: FM-Index DNA backward-search customized for GEM Mapper
 */

#include "../include/gpu_fmi_core.h"


void __global__ gpu_fmi_search_kernel(const gpu_fmi_device_entry_t *fmi, const uint64_t bwtSize,
									  const uint32_t numSeeds, const ulonglong2 *seeds, ulonglong2 *resIntervals)
{
	const uint32_t globalThreadIdx     = blockIdx.x * GPU_MAX_THREADS_PER_SM + threadIdx.x;
	const uint32_t localWarpThreadIdx  = globalThreadIdx     % GPU_WARP_SIZE;
	const uint32_t localEntryIdx       = localWarpThreadIdx  / GPU_FMI_THREADS_PER_ENTRY;
	const uint32_t localEntryThreadIdx = localWarpThreadIdx  % GPU_FMI_THREADS_PER_ENTRY;
	const uint32_t idSeed 	 		   = globalThreadIdx     / GPU_FMI_SEED_THREADS_PER_ENTRY;

	if ( (threadIdx.x < GPU_MAX_THREADS_PER_SM) && (globalThreadIdx < (numSeeds * GPU_FMI_SEED_THREADS_PER_ENTRY)) ){

		const uint32_t   localIdSeed = idSeed % GPU_FMI_SEED_ENTRIES_PER_WARP;
		const ulonglong2 seed 	     = seeds[idSeed];
		const uint32_t   seedSize    = seed.y >> (GPU_UINT64_LENGTH - GPU_FMI_SEED_FIELD_SIZE);
			  uint64_t   currentSeed = seed.x;

			  uint64_t sharedInterval, interval = (localEntryIdx % GPU_FMI_ENTRIES_PER_SEED) ? bwtSize : 0;
			  uint32_t idStep = 0, foundSeed = 0;

		__shared__ gpu_fmi_exch_bmp_mem_t   exchBMP[GPU_FMI_ENTRIES_PER_BLOCK];
				   gpu_fmi_exch_bmp_mem_t * seedExchBMP = &exchBMP[threadIdx.x / GPU_FMI_THREADS_PER_ENTRY];

		while((idStep < seedSize) && (foundSeed == 0)){			
			const uint64_t entryIdx    	  =  interval / GPU_FMI_ENTRY_SIZE;
			const uint32_t bitmapPosition =  interval % GPU_FMI_ENTRY_SIZE;

			// Gathering the base of the seed
			currentSeed = (idStep == GPU_FMI_SEED_BASES_PER_ENTRY) ? seed.y : currentSeed;
			const uint32_t bit0 =  currentSeed & 0x1L;
			const uint32_t bit1 = (currentSeed & 0x2L) >> 1;
			currentSeed >>= GPU_FMI_SEED_CHAR_LENGTH;
 
			// Loading FM-index entry in thread cooperative way
			const uint32_t missedEntry   = (entryIdx % GPU_FMI_ALTERNATE_COUNTERS != bit1) ? 1 : 0;
			const uint64_t entryIdxFixed = (localEntryThreadIdx == 0) ? entryIdx + missedEntry : entryIdx;
			uint4 loadEntry              = fmi[entryIdxFixed].v[localEntryThreadIdx];

			// Compute LF-Mapping (th0 of each group contain the result)
 			interval = LF_Mapping(loadEntry, seedExchBMP, missedEntry, localEntryThreadIdx, bitmapPosition, bit1, bit0);

			// Shared results
			const uint32_t lane = (localEntryThreadIdx == 0) ? localWarpThreadIdx + GPU_FMI_THREADS_PER_ENTRY : localEntryIdx * GPU_FMI_THREADS_PER_ENTRY;
			sharedInterval      = shfl_64(interval, lane);

			// Early exit condition
			foundSeed           = (__ballot(interval == sharedInterval) >> (localIdSeed * GPU_FMI_SEED_THREADS_PER_ENTRY)) & GPU_UINT32_MASK_ONE_LOW;

			// Update interval for bitmap threads
			if (localEntryThreadIdx) interval = sharedInterval;

			// Increment for the next Backward-Search
			idStep++;
		}
		// Save intervals
		if((localWarpThreadIdx % GPU_FMI_SEED_THREADS_PER_ENTRY) == 0)
			resIntervals[idSeed] = make_ulonglong2(interval, sharedInterval);
	}
}

extern "C"
gpu_error_t gpu_fmi_search_launch_kernel(const gpu_fmi_device_entry_t *d_fmi, const uint64_t bwtSize, const uint32_t numSeeds, const ulonglong2 *d_seeds, ulonglong2 *d_intervals)
{
	const uint32_t threads = 128;
	const uint32_t blocks  = GPU_DIV_CEIL(numSeeds * GPU_FMI_SEED_THREADS_PER_ENTRY, threads);
	const uint32_t nreps   = 10;

	float elapsed_time_ms = 0.0f;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

		for(uint32_t iteration = 0; iteration < nreps; ++iteration)
			gpu_fmi_search_kernel<<<blocks,threads>>>(d_fmi, bwtSize, numSeeds, d_seeds, d_intervals);

	cudaEventRecord(stop, 0);
	cudaThreadSynchronize();
	cudaEventElapsedTime(&elapsed_time_ms, start, stop);
	elapsed_time_ms /= nreps;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	printf("\t Time Kernel GPU:  %8.2f ms\n", elapsed_time_ms);

	return(SUCCESS);
}

extern "C"
gpu_error_t gpu_fmi_search_process_buffer(gpu_buffer_t *mBuff)
{
	gpu_index_buffer_t 			  	 *index    	  =  mBuff->index;
	gpu_fmi_search_seeds_buffer_t 	 *seeds    	  = &mBuff->data.search.seeds;
	gpu_fmi_search_sa_inter_buffer_t *saIntervals = &mBuff->data.search.saIntervals;
	uint32_t 					     numSeeds	  =  mBuff->data.search.seeds.numSeeds;
	cudaStream_t 				     idStream	  =  mBuff->idStream;
	uint32_t					     idSupDev	  =  mBuff->idSupportedDevice;

	uint32_t threadsPerBlock = GPU_MAX_THREADS_PER_BLOCK;
	uint32_t numThreads = numSeeds * GPU_FMI_SEED_THREADS_PER_ENTRY;
	uint32_t blocksPerGrid = GPU_DIV_CEIL(numThreads, threadsPerBlock);

	//printf("KEPLER 2ndGen: LAUNCH KERNEL -- Bloques: %d - Th_block %d\n", blocksPerGrid, threadsPerBlock);
	gpu_fmi_search_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((gpu_fmi_device_entry_t*) index->d_fmi[idSupDev], index->bwtSize,
																		   seeds->numSeeds, (ulonglong2*) seeds->d_seeds,
																		   (ulonglong2*) saIntervals->d_intervals);
	return(SUCCESS);
}



