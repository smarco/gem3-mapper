/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Host scheduler for BPM on GPU
 */

#include "../include/gpu_fmi.h"

/************************************************************
Functions to get the GPU FMI buffers
************************************************************/

GPU_INLINE gpu_fmi_entry_t* gpu_fmi_buffer_get_index_(void *fmiBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	return(mBuff->index->h_fmi);
}

GPU_INLINE gpu_fmi_search_seed_t* gpu_fmi_search_buffer_get_seeds_(void *fmiBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	return(mBuff->data.search.seeds.h_seeds);
}

GPU_INLINE gpu_fmi_search_sa_inter_t* gpu_fmi_search_buffer_get_sa_intervals_(void *fmiBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	return(mBuff->data.search.saIntervals.h_intervals);
}

GPU_INLINE gpu_fmi_decode_init_pos_t* gpu_fmi_decode_buffer_get_init_pos_(void *fmiBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	return(mBuff->data.decode.initPositions.h_initBWTPos);
}

GPU_INLINE gpu_fmi_decode_end_pos_t* gpu_fmi_decode_buffer_get_end_pos_(void *fmiBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	return(mBuff->data.decode.endPositions.h_endBWTPos);
}


/************************************************************
Functions to get the maximum elements of the buffers
************************************************************/

GPU_INLINE uint32_t gpu_fmi_search_buffer_get_max_seeds_(void *fmiBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	return(mBuff->data.search.numMaxSeeds);
}

GPU_INLINE uint32_t gpu_fmi_decode_buffer_get_max_positions_(void *fmiBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	return(mBuff->data.decode.numMaxInitPositions);
}

/************************************************************
Functions to init the buffers (E. SEARCH)
************************************************************/

GPU_INLINE size_t gpu_fmi_search_input_size()
{
	return(sizeof(gpu_fmi_search_seed_t) + sizeof(gpu_fmi_search_sa_inter_t));
}

GPU_INLINE void gpu_fmi_search_reallocate_host_buffer_layout(gpu_buffer_t* mBuff)
{
	void* rawAlloc = mBuff->h_rawData;
	//Adjust the host buffer layout (input)
	mBuff->data.search.seeds.h_seeds = GPU_ALIGN_TO(rawAlloc,16);
	rawAlloc = (void *) (mBuff->data.search.seeds.h_seeds + mBuff->data.search.numMaxSeeds);
	//Adjust the host buffer layout (output)
	mBuff->data.search.saIntervals.h_intervals = GPU_ALIGN_TO(rawAlloc,16);
	rawAlloc = (void *) (mBuff->data.search.saIntervals.h_intervals + mBuff->data.search.numMaxIntervals);
}

GPU_INLINE void gpu_fmi_search_reallocate_device_buffer_layout(gpu_buffer_t* mBuff)
{
	void* rawAlloc = mBuff->d_rawData;
	//Adjust the host buffer layout (input)
	mBuff->data.search.seeds.d_seeds = GPU_ALIGN_TO(rawAlloc,16);
	rawAlloc = (void *) (mBuff->data.search.seeds.d_seeds + mBuff->data.search.numMaxSeeds);
	//Adjust the host buffer layout (output)
	mBuff->data.search.saIntervals.d_intervals = GPU_ALIGN_TO(rawAlloc,16);
	rawAlloc = (void *) (mBuff->data.search.saIntervals.d_intervals + mBuff->data.search.numMaxIntervals);
}

GPU_INLINE void gpu_fmi_search_init_buffer_(void* fmiBuffer)
{
	gpu_buffer_t *mBuff      = (gpu_buffer_t *) fmiBuffer;
	void* 		    	 h_rawAlloc = mBuff->h_rawData;
	void* 		  		 d_rawAlloc = mBuff->d_rawData;
	size_t	  	  		 sizeBuff   = mBuff->sizeBuffer;
	uint32_t	  		 numInputs  = (sizeBuff / gpu_fmi_search_input_size()) - GPU_FMI_SEARCH_SEEDS_BUFFER_PADDING;

	//Set real size of the input
	mBuff->data.search.numMaxSeeds     = numInputs;
	mBuff->data.search.numMaxIntervals = numInputs;
	gpu_fmi_search_reallocate_host_buffer_layout(mBuff);
	gpu_fmi_search_reallocate_device_buffer_layout(mBuff);
}


/************************************************************
Functions to init the buffers (DECODE)
************************************************************/

GPU_INLINE size_t gpu_fmi_decode_input_size()
{
	return(sizeof(gpu_fmi_decode_init_pos_t) + sizeof(gpu_fmi_decode_end_pos_t));
}

GPU_INLINE void gpu_fmi_decode_reallocate_host_buffer_layout(gpu_buffer_t* mBuff)
{
	void* rawAlloc = mBuff->h_rawData;
	//Adjust the host buffer layout (input)
	mBuff->data.decode.initPositions.h_initBWTPos = GPU_ALIGN_TO(rawAlloc,16);
	rawAlloc = (void *) (mBuff->data.decode.initPositions.h_initBWTPos + mBuff->data.decode.numMaxInitPositions);
	//Adjust the host buffer layout (output)
	mBuff->data.decode.endPositions.h_endBWTPos   = GPU_ALIGN_TO(rawAlloc,16);
	rawAlloc = (void *) (mBuff->data.decode.endPositions.h_endBWTPos + mBuff->data.decode.numMaxEndPositions);
}

GPU_INLINE void gpu_fmi_decode_reallocate_device_buffer_layout(gpu_buffer_t* mBuff)
{
	void* rawAlloc = mBuff->d_rawData;
	//Adjust the host buffer layout (input)
	mBuff->data.decode.initPositions.d_initBWTPos = GPU_ALIGN_TO(rawAlloc,16);
	rawAlloc = (void *) (mBuff->data.decode.initPositions.d_initBWTPos + mBuff->data.decode.numMaxInitPositions);
	//Adjust the host buffer layout (output)
	mBuff->data.decode.endPositions.d_endBWTPos   = GPU_ALIGN_TO(rawAlloc,16);
	rawAlloc = (void *) (mBuff->data.decode.endPositions.d_endBWTPos + mBuff->data.decode.numMaxEndPositions);
}

GPU_INLINE void gpu_fmi_decode_init_buffer_(void* fmiBuffer)
{
	gpu_buffer_t *mBuff      = (gpu_buffer_t *) fmiBuffer;
	size_t	  	  sizeBuff   = mBuff->sizeBuffer;
	uint32_t	  numMaxPositions  = (sizeBuff / gpu_fmi_decode_input_size()) - GPU_FMI_DECODE_POS_BUFFER_PADDING;

	//Set real size of the input
	mBuff->data.decode.numMaxInitPositions = numMaxPositions;
	mBuff->data.decode.numMaxEndPositions = numMaxPositions;
	gpu_fmi_decode_reallocate_host_buffer_layout(mBuff);
	gpu_fmi_decode_reallocate_device_buffer_layout(mBuff);
}


/************************************************************
Functions to transfer data HOST <-> DEVICE (E. SEARCH)
************************************************************/

GPU_INLINE gpu_error_t gpu_fmi_search_transfer_CPU_to_GPU(gpu_buffer_t *mBuff)
{
	const gpu_fmi_search_seeds_buffer_t* seedBuff = &mBuff->data.search.seeds;
	const cudaStream_t 				     idStream =  mBuff->idStream;
	const size_t 					     cpySize  =  seedBuff->numSeeds * sizeof(gpu_fmi_search_seed_t);

	//Transfer seeds from CPU to the GPU
	CUDA_ERROR(cudaMemcpyAsync(seedBuff->d_seeds, seedBuff->h_seeds, cpySize, cudaMemcpyHostToDevice, idStream));

	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_fmi_search_transfer_GPU_to_CPU(gpu_buffer_t *mBuff)
{
	const gpu_fmi_search_sa_inter_buffer_t* interBuff = &mBuff->data.search.saIntervals;
	const cudaStream_t 				   		idStream  =  mBuff->idStream;
	const size_t 							cpySize   =  interBuff->numIntervals * sizeof(gpu_fmi_search_sa_inter_t);

	//Transfer SA intervals (occurrence results) from CPU to the GPU
	CUDA_ERROR(cudaMemcpyAsync(interBuff->h_intervals, interBuff->d_intervals, cpySize, cudaMemcpyDeviceToHost, idStream));

	return (SUCCESS);
}

GPU_INLINE void gpu_fmi_search_send_buffer_(void* fmiBuffer, const uint32_t numSeeds)
{
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	const uint32_t idSupDevice = mBuff->idSupportedDevice;

	//Set real size of the input
	mBuff->data.search.seeds.numSeeds = numSeeds;
	mBuff->data.search.saIntervals.numIntervals = numSeeds;

	//Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));
	GPU_ERROR(gpu_fmi_search_transfer_CPU_to_GPU(mBuff));
	GPU_ERROR(gpu_fmi_search_process_buffer(mBuff));
	GPU_ERROR(gpu_fmi_search_transfer_GPU_to_CPU(mBuff));
}

GPU_INLINE void gpu_fmi_search_receive_buffer_(void *fmiBuffer)
{
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	uint32_t error;

	//Synchronize Stream (the thread wait for the commands done in the stream)
	CUDA_ERROR(cudaStreamSynchronize(mBuff->idStream));
}

/************************************************************
Functions to transfer data HOST <-> DEVICE (Decode)
************************************************************/

GPU_INLINE gpu_error_t gpu_fmi_decode_transfer_CPU_to_GPU(gpu_buffer_t *mBuff)
{
	const gpu_fmi_decode_init_pos_buffer_t* initPosBuff = &mBuff->data.decode.initPositions;
	const cudaStream_t 				   	    idStream    =  mBuff->idStream;
	const size_t 							cpySize     =  initPosBuff->numDecodings * sizeof(gpu_fmi_decode_init_pos_t);

	//Transfer seeds from CPU to the GPU
	CUDA_ERROR(cudaMemcpyAsync(initPosBuff->d_initBWTPos, initPosBuff->h_initBWTPos, cpySize, cudaMemcpyHostToDevice, idStream));

	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_fmi_decode_transfer_GPU_to_CPU(gpu_buffer_t *mBuff)
{
	const gpu_fmi_search_sa_inter_buffer_t* interBuff = &mBuff->data.search.saIntervals;
	const cudaStream_t 				   		idStream  =  mBuff->idStream;
	const size_t 							cpySize   =  interBuff->numIntervals * sizeof(gpu_fmi_search_sa_inter_t);

	//Transfer SA intervals (occurrence results) from CPU to the GPU
	CUDA_ERROR(cudaMemcpyAsync(interBuff->h_intervals, interBuff->d_intervals, cpySize, cudaMemcpyDeviceToHost, idStream));

	return (SUCCESS);
}

GPU_INLINE void gpu_fmi_decode_send_buffer_(void* fmiBuffer, const uint32_t numDecodings, const uint32_t samplingRate)
{
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	const uint32_t idSupDevice = mBuff->idSupportedDevice;

	//Set real size of the input
	mBuff->data.decode.initPositions.numDecodings = numDecodings;
	mBuff->data.decode.endPositions.numDecodings  = numDecodings;
	mBuff->data.decode.samplingRate    = samplingRate;

	//Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

	GPU_ERROR(gpu_fmi_decode_transfer_CPU_to_GPU(mBuff));
	GPU_ERROR(gpu_fmi_decode_process_buffer(mBuff));
	GPU_ERROR(gpu_fmi_decode_transfer_GPU_to_CPU(mBuff));
}

GPU_INLINE void gpu_fmi_decode_receive_buffer_(void *fmiBuffer)
{
	gpu_buffer_t *mBuff = (gpu_buffer_t *) fmiBuffer;
	uint32_t error;

	//Synchronize Stream (the thread wait for the commands done in the stream)
	CUDA_ERROR(cudaStreamSynchronize(mBuff->idStream));
}
