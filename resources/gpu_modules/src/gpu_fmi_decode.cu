/*
 * PROJECT: Thread-cooperative FM-index on GPU
 * FILE: genIndex.h
 * DATE: 1/9/2015
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: FM-Index DNA backward-search customized for GEM Mapper
 */

#include "../include/gpu_fmi_core.h"

void __global__ gpu_fmi_decoding_positions_kernel(const gpu_fmi_device_entry_t *fmi, const uint64_t bwtSize, const uint32_t numDecodings,
											  	  const uint64_t *d_initBWTPos, ulonglong2 *d_endBWTPos, const uint32_t samplingRate)
{
	const uint32_t globalThreadIdx      = blockIdx.x * GPU_MAX_THREADS_PER_SM + threadIdx.x;
	const uint32_t localWarpThreadIdx   = globalThreadIdx % GPU_WARP_SIZE;
	const uint32_t idDecoding 	 	    = globalThreadIdx / GPU_FMI_DECODE_THREADS_PER_ENTRY;

	if((threadIdx.x < GPU_MAX_THREADS_PER_SM) && (idDecoding < numDecodings)){

		const uint32_t decodeEntryIdx        = localWarpThreadIdx  / GPU_FMI_DECODE_THREADS_PER_ENTRY;
		const uint32_t decodeEntryThreadIdx  = localWarpThreadIdx  % GPU_FMI_DECODE_THREADS_PER_ENTRY;
		const uint32_t fmiEntryThreadIdx     = localWarpThreadIdx  % GPU_FMI_THREADS_PER_ENTRY;
		const uint32_t loadFMIEntryThreadIdx = generateLoadThreadIdx(decodeEntryThreadIdx);

			  uint4    loadEntry;
			  uint64_t interval = d_initBWTPos[idDecoding], idStep = 0;
			  uint32_t foundBaseN = 0;

		__shared__ gpu_fmi_exch_bmp_mem_t   exchBMP[GPU_FMI_ENTRIES_PER_BLOCK];
				   gpu_fmi_exch_bmp_mem_t * decExchBMP = &exchBMP[threadIdx.x / GPU_FMI_THREADS_PER_ENTRY];

		while((interval % samplingRate) && (foundBaseN == 0)){
			const uint64_t entryIdx    	  =  interval / GPU_FMI_ENTRY_SIZE;
			const uint32_t bitmapPosition =  interval % GPU_FMI_ENTRY_SIZE;
			      uint32_t bit1, bit0, bit2;

			// Loading FM-index entry in thread cooperative way
			if(( 0 < decodeEntryThreadIdx) && (decodeEntryThreadIdx < GPU_FMI_DECODE_THREADS_PER_LOAD))
				loadEntry = fmi[entryIdx].v[loadFMIEntryThreadIdx];

			// Gathering the base and sharing it with the rest of the threads
			gather_base_from_BWT(loadEntry, decExchBMP, bitmapPosition, fmiEntryThreadIdx, decodeEntryThreadIdx, decodeEntryIdx, &bit1, &bit0, &bit2);

			// Gathering the counters
			const uint32_t missedEntry = (entryIdx % GPU_FMI_ALTERNATE_COUNTERS != bit1) ? 1 : 0;
 			loadEntry = gather_counters_FMI(loadEntry, missedEntry, decodeEntryIdx, decodeEntryThreadIdx);

			// Compute LF-Mapping th0 of each group contain the result)
 			interval = LF_Mapping(loadEntry, decExchBMP, missedEntry, fmiEntryThreadIdx, bitmapPosition, bit1, bit0);

			// Share interval along the thread group
 			const uint32_t lane = decodeEntryIdx * GPU_FMI_DECODE_THREADS_PER_ENTRY;
			interval 		    = shfl_64(interval, lane);

			// Exit condition
			if(bit2 == 0) foundBaseN = 1;

			// Increment for the next Backward-Search
			idStep++;
		}

		ulonglong2 BWTPos = {interval, idStep};
		if(foundBaseN) BWTPos = make_ulonglong2(GPU_UINT64_MAX_VALUE, GPU_UINT64_MAX_VALUE);

		// Save intervals
		if(decodeEntryThreadIdx  == 0) d_endBWTPos[idDecoding] = BWTPos;
	}
}

extern "C"
gpu_error_t gpu_fmi_decoding_launch_kernel(gpu_fmi_device_entry_t *d_fmi, uint64_t bwtSize, uint32_t numDecodings, uint64_t *d_initBWTPos, ulonglong2 *d_endBWTPos)
{
	const uint32_t samplingRate = 4;
	const uint32_t threads = 128;
	const uint32_t blocks  = GPU_DIV_CEIL(numDecodings * GPU_FMI_DECODE_THREADS_PER_ENTRY, threads);
	const uint32_t nreps   = 10;

	printf("GPU numDecodings=%u \n", numDecodings);

	float elapsed_time_ms = 0.0f;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

		for(uint32_t iteration = 0; iteration < nreps; ++iteration)
			gpu_fmi_decoding_positions_kernel<<<blocks,threads>>>(d_fmi, bwtSize, numDecodings, d_initBWTPos, d_endBWTPos, samplingRate);

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
gpu_error_t gpu_fmi_decode_process_buffer(gpu_buffer_t *mBuff)
{
	gpu_index_buffer_t 			  	 *index    	  = mBuff->index;
	gpu_fmi_decode_buffer_t			 *decBuff	  = mBuff->data.decode;
	gpu_fmi_decode_init_pos_buffer_t *initPos     = mBuff->data.decode.initPositions;
	gpu_fmi_decode_end_pos_buffer_t  *endPos 	  = mBuff->data.decode.endPositions;
	uint32_t 					     numDecodings = mBuff->data.decode.initPositions.numDecodings;
	cudaStream_t 				     idStream	  = mBuff->idStream;
	uint32_t					     idSupDev 	  = mBuff->device->idSupportedDevice;

	const uint32_t threadsPerBlock = GPU_MAX_THREADS_PER_BLOCK;
	const uint32_t numThreads = numDecodings * GPU_FMI_DECODE_THREADS_PER_ENTRY;
	const uint32_t blocksPerGrid = GPU_DIV_CEIL(numThreads, threadsPerBlock);

	//printf("KEPLER 2ndGen: LAUNCH KERNEL -- Bloques: %d - Th_block %d\n", blocksPerGrid, threadsPerBlock);
	gpu_fmi_decoding_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((gpu_fmi_device_entry_t*) index->d_fmi[idSupDev], index->bwtSize,
																			 initPos->numDecodings, (uint64_t*) initPos->d_initBWTPos,
																			 (ulonglong2*) endPos->d_endBWTPos, decBuff->samplingRate);

	return(SUCCESS);
}



