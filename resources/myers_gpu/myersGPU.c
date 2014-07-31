#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "myers-common.h"


#define		BASES_PER_THREAD		128
#define		BASES_PER_PEQ_ENTRY		128

#define		SIZE_WARP				32
#define		NUMBUCKETS				(SIZE_WARP + 1)


#ifndef 	DEVICE
			#define 	DEVICE 		0
#endif



/************************************************************
Functions to handle errors
************************************************************/

void CudaError( cudaError_t err, const char *file,  int line ) {
   	if (err != cudaSuccess) {
      		fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err),  file, line );
       		exit(EXIT_FAILURE);
   	}
}

char *myersGetErrorString(myersError_t error){
    switch(error) {
        case E_OPENING_FILE:  		 return "Myers - Error: opening file"; break;
        case E_READING_FILE:  		 return "Myers - Error: reading file"; break;
        case E_INSUFFICIENT_MEM_GPU: return "Myers - Error: there aren't enough GPU memory space"; break;
        case E_ALLOCATE_MEM: 		 return "Myers - Error: allocating data"; break;
        case E_INCOMPATIBLE_GPU:	 return "Myers - Error: incompatible GPU (old CC version)"; break;
        default: 					 return "Myers - Unknown error";
    }
}

void MyersError( myersError_t err, const char *file,  int line ) {
   	if (err != 0) {
      		fprintf(stderr, "%s in %s at line %d\n", myersGetErrorString(err),  file, line );
       		exit(EXIT_FAILURE);
   	}
}



/************************************************************
Functions to get the Myers buffers
************************************************************/

inline uint32_t getMaxPEQEntries(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->maxPEQEntries);
}

inline uint32_t getMaxCandidates(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->maxCandidates);
}

inline uint32_t getMaxQueries(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->maxQueries);
}

inline qryEntry_t* getPEQBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->queries->h_queries);
}

inline candInfo_t* getCandidatesBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->candidates->h_candidates);
}

inline qryInfo_t* getPEQInfoBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->queries->h_qinfo);
}

inline resEntry_t* getResultsBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->results->h_results);
}



/************************************************************
Functions to init all the Myers resources
************************************************************/

myersError_t freeReferenceHost(reference_buffer_t *reference)
{
    if(reference->h_reference != NULL){
        CUDA_ERROR(cudaFreeHost(reference->h_reference));
        reference->h_reference = NULL;
    }

    return(0);
}

myersError_t readReference(const char *fn, reference_buffer_t *reference)
{
	FILE *fp = NULL;
	size_t result;

	fp = fopen(fn, "rb");
		if (fp == NULL) return (E_OPENING_FILE);

	//TODO: Cambiar a 8bytes para direccionar tamaños mayores de 4 Giga entradas
    result = fread(&reference->numEntries, sizeof(uint32_t), 1, fp);
		if (result != 1) return (E_READING_FILE);
    result = fread(&reference->size, sizeof(uint32_t), 1, fp);
		if (result != 1) return (E_READING_FILE);

    //TODO: Cambiar reserva en el host para referencias en streaming

	CUDA_ERROR(cudaHostAlloc((void**)&reference->h_reference, reference->numEntries * sizeof(uint32_t), cudaHostAllocMapped));
	result = fread(reference->h_reference, sizeof(uint32_t), reference->numEntries, fp);
		if (result != reference->numEntries) return (E_READING_FILE);

	fclose(fp);
	return (SUCCESS);
}

myersError_t initReference(const char *refFile, reference_buffer_t **reference)
{
	size_t free, total;
	reference_buffer_t *ref = (reference_buffer_t *) malloc(sizeof(reference_buffer_t));
	uint32_t error;

    ref->numEntries = 0;
    ref->size = 0;

	ref->h_reference = NULL;
	ref->d_reference = NULL;

	MYERS_ERROR(readReference(refFile, ref));

	CUDA_ERROR(cudaMemGetInfo(&free, &total));
	if ((ref->numEntries * sizeof(uint32_t)) > free) return(E_INSUFFICIENT_MEM_GPU);

	//Synchronous allocate & transfer Binary Reference to GPU
	CUDA_ERROR(cudaMalloc((void**) &ref->d_reference, ref->numEntries * sizeof(uint32_t)));
	CUDA_ERROR(cudaMemcpy(ref->d_reference, ref->h_reference, ref->numEntries * sizeof(uint32_t), cudaMemcpyHostToDevice));

	MYERS_ERROR(freeReferenceHost(ref));
	(* reference) = ref;
	return (SUCCESS);
}

myersError_t initQueries(queries_buffer_t *queries, uint32_t maxPEQEntries, uint32_t maxQueries)
{
	queries->totalQueriesEntries = 0;
	queries->numQueries = 0;
	queries->h_queries = NULL;
	queries->d_queries = NULL;
	queries->h_qinfo = NULL;
	queries->d_qinfo = NULL;

	//Allocate PEQ entries in CPU & GPU
	CUDA_ERROR(cudaHostAlloc((void**) &queries->h_queries, maxPEQEntries * sizeof(qryEntry_t), cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &queries->d_queries, maxPEQEntries * sizeof(qryEntry_t)));

	//Allocate queries info in CPU & GPU
	CUDA_ERROR(cudaHostAlloc((void**) &queries->h_qinfo, maxQueries * sizeof(qryInfo_t), cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &queries->d_qinfo, maxQueries * sizeof(qryInfo_t)));

	return (SUCCESS);
}

myersError_t initCandidates(candidates_buffer_t *candidates, uint32_t maxCandidates)
{
	candidates->numCandidates = 0;
	candidates->h_candidates = NULL;
	candidates->d_candidates = NULL;

	//Allocate candidates info in CPU & GPU
	//candidates->h_candidates = (candInfo_t *) malloc(maxCandidates * sizeof(candInfo_t));
	//	if (candidates->h_candidates == NULL) return (34);
	CUDA_ERROR(cudaHostAlloc((void**) &candidates->h_candidates, maxCandidates * sizeof(candInfo_t), cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &candidates->d_candidates, maxCandidates * sizeof(candInfo_t)));

	return (SUCCESS);
}

myersError_t initReorderBuffer(reorder_buffer_t *reorderBuffer, uint32_t maxReorderBuffer)
{
	reorderBuffer->numBuckets = NUMBUCKETS;
	reorderBuffer->candidatesPerBufffer = 0;
	reorderBuffer->numWarps = 0;

	reorderBuffer->h_reorderBuffer		= NULL;
	reorderBuffer->d_reorderBuffer		= NULL;
	reorderBuffer->h_initPosPerBucket	= NULL;
	reorderBuffer->d_initPosPerBucket	= NULL;
	reorderBuffer->h_initWarpPerBucket	= NULL;
	reorderBuffer->d_initWarpPerBucket	= NULL;

	//Allocate Binary Reference in GPU
	CUDA_ERROR(cudaHostAlloc((void**) &reorderBuffer->h_reorderBuffer, maxReorderBuffer * sizeof(uint32_t), cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &reorderBuffer->d_reorderBuffer, maxReorderBuffer * sizeof(uint32_t)));

	CUDA_ERROR(cudaHostAlloc((void**) &reorderBuffer->h_initPosPerBucket, reorderBuffer->numBuckets * sizeof(uint32_t), cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &reorderBuffer->d_initPosPerBucket, reorderBuffer->numBuckets * sizeof(uint32_t)));

	CUDA_ERROR(cudaHostAlloc((void**) &reorderBuffer->h_initWarpPerBucket, reorderBuffer->numBuckets * sizeof(uint32_t), cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &reorderBuffer->d_initWarpPerBucket, reorderBuffer->numBuckets * sizeof(uint32_t)));

	return (SUCCESS);
}

myersError_t initResults(results_buffer_t *results, uint32_t maxReorderBuffer, uint32_t maxCandidates)
{
	results->numResults = 0;
	results->numReorderedResults = 0;

	results->h_results			= NULL;
	results->h_reorderResults	= NULL;
	results->d_reorderResults	= NULL;

	//Allocate candidates info in CPU
	results->h_results = (resEntry_t *) malloc(maxCandidates * sizeof(resEntry_t));
		if (results->h_results == NULL) return (E_ALLOCATE_MEM);

	//Allocate Binary Reference in GPU
	CUDA_ERROR(cudaHostAlloc((void**) &results->h_reorderResults, maxReorderBuffer * sizeof(resEntry_t), cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &results->d_reorderResults, maxReorderBuffer * sizeof(resEntry_t)));

	return (SUCCESS);
}

myersError_t initBuffer(buffer_t *buffer, uint32_t idBuffer, uint32_t idDevice, uint32_t numBuffers,
			   uint32_t maxCandidates, uint32_t maxQueries, uint32_t maxPEQEntries,
			   uint32_t bucketPaddingCandidates, reference_buffer_t *reference)
{
    struct cudaDeviceProp deviceProp;

	buffer->queries			=	NULL;
	buffer->candidates		=	NULL;
	buffer->reorderBuffer	=	NULL;
	buffer->results 		=	NULL;

	buffer->numBuffers 		= 	numBuffers;
	buffer->idBuffer 		=	idBuffer;
	buffer->maxPEQEntries	=	maxPEQEntries;
	buffer->maxCandidates	=	maxCandidates;
	buffer->maxQueries		=	maxQueries;
	buffer->maxReorderBuffer=	maxCandidates + bucketPaddingCandidates;

	buffer->reference 		=	reference;

	buffer->queries			=	(queries_buffer_t *)    malloc(sizeof(queries_buffer_t));
		if (buffer->queries == NULL) return (E_ALLOCATE_MEM);
	buffer->candidates		=	(candidates_buffer_t *) malloc(sizeof(candidates_buffer_t));
		if (buffer->candidates == NULL) return (E_ALLOCATE_MEM);
	buffer->reorderBuffer	=	(reorder_buffer_t *)    malloc(sizeof(reorder_buffer_t));
		if (buffer->reorderBuffer == NULL) return (E_ALLOCATE_MEM);
	buffer->results 		=	(results_buffer_t *)    malloc(sizeof(results_buffer_t));
		if (buffer->results == NULL) return (E_ALLOCATE_MEM);

	buffer->idDevice = idDevice;
    CUDA_ERROR(cudaGetDeviceProperties(&deviceProp, buffer->idDevice));
    //printf("\nDevice %d has compute capability %d.%d.\n", DEVICE, deviceProp.major, deviceProp.minor);
    buffer->majorCC = deviceProp.major;
    buffer->minorCC = deviceProp.minor;

	//Create the cuda stream per each buffer
	CUDA_ERROR(cudaStreamCreate(&buffer->idStream));

	MYERS_ERROR(initQueries(buffer->queries, buffer->maxPEQEntries, buffer->maxQueries));
	MYERS_ERROR(initCandidates(buffer->candidates, buffer->maxCandidates));
	MYERS_ERROR(initReorderBuffer(buffer->reorderBuffer, buffer->maxReorderBuffer));
	MYERS_ERROR(initResults(buffer->results, buffer->maxReorderBuffer, buffer->maxCandidates));

	return(SUCCESS);
}


void initMyers(void ***myersBuffer, uint32_t numBuffers, const char *referenceFileName,
			  int32_t averageQuerySize, int32_t candidatesPerQuery)
{
	buffer_t **buffer = NULL;
	reference_buffer_t *reference;

	size_t freeMem, totalMem, bytesPerBuffer;
	uint32_t idBuffer, idBucket, bucketPaddingCandidates;
	uint32_t maxCandidates, maxQueries, maxPEQEntries;
	uint32_t numBuffersPerDevice, idDevice, idGlobalBuffer, idLocalBuffer;
	int numDevices;

	buffer = (buffer_t **) malloc(numBuffers * sizeof(buffer_t *));
		if (buffer == NULL) MYERS_ERROR(E_ALLOCATE_MEM);

	//Worst case space
	if(averageQuerySize <= 0) averageQuerySize = 1024;
	if(candidatesPerQuery <= 0) candidatesPerQuery = 1;

	CUDA_ERROR(cudaGetDeviceCount(&numDevices));
	//printf("numDevices: %d\n", numDevices);
	numBuffersPerDevice = numBuffers / numDevices;

	for(idDevice = 0; idDevice < numDevices; idDevice++){
		//printf("Identidicador de GPU %d \n", idDevice);
		CUDA_ERROR(cudaSetDevice(idDevice));
		MYERS_ERROR(initReference(referenceFileName, &reference));

		uint32_t averageNumPEQEntries = (averageQuerySize / BASES_PER_PEQ_ENTRY) +
										((averageQuerySize % BASES_PER_PEQ_ENTRY) ? 1 : 0);

		float bytesPerQuery = (((averageNumPEQEntries * sizeof(qryEntry_t))
								+ sizeof(qryInfo_t)) / (float) candidatesPerQuery)
						 		+ sizeof(candInfo_t) + sizeof(resEntry_t) + sizeof(uint32_t);

		bucketPaddingCandidates = 0;
		for(idBucket = 1; idBucket < NUMBUCKETS-1; idBucket++)
			bucketPaddingCandidates += (SIZE_WARP / idBucket);

		CUDA_ERROR(cudaMemGetInfo(&freeMem, &totalMem));

		bytesPerBuffer = (size_t) ((freeMem - (freeMem * 0.05)) / numBuffersPerDevice)
									- (bucketPaddingCandidates * sizeof(uint32_t))
									- (bucketPaddingCandidates * sizeof(resEntry_t));
		maxCandidates = (uint32_t) (bytesPerBuffer / bytesPerQuery);
		maxQueries 	  = maxCandidates / candidatesPerQuery;
		maxPEQEntries = maxQueries * averageNumPEQEntries;

		if ((maxCandidates < 1) || (maxQueries < 1) || ((bytesPerBuffer * numBuffersPerDevice) > freeMem))
			MYERS_ERROR(E_INSUFFICIENT_MEM_GPU);

		for(idLocalBuffer = 0; idLocalBuffer < numBuffersPerDevice; idLocalBuffer++){
			idGlobalBuffer = idDevice * numBuffersPerDevice + idLocalBuffer;
			buffer[idGlobalBuffer] = (buffer_t *) malloc(sizeof(buffer_t));
			MYERS_ERROR(initBuffer(buffer[idGlobalBuffer], idBuffer, idDevice, numBuffersPerDevice, maxCandidates,
					   				maxQueries, maxPEQEntries, bucketPaddingCandidates, reference));
		}
	}

	(* myersBuffer) = (void **) buffer;
}


/************************************************************
Funtions to send & process a Myers buffer to GPU
************************************************************/

myersError_t reorderingBuffer(buffer_t *mBuff)
{

	queries_buffer_t  	*qry  		= mBuff->queries;
	candidates_buffer_t	*cand 		= mBuff->candidates;
	reorder_buffer_t  	*rebuff 	= mBuff->reorderBuffer;
	results_buffer_t  	*res  		= mBuff->results;

	uint32_t idBucket, idCandidate, idBuff;
	uint32_t numThreadsPerQuery;
	uint32_t numQueriesPerWarp;
	uint32_t tmpBuckets[rebuff->numBuckets];
	uint32_t numCandidatesPerBucket[rebuff->numBuckets];
	uint32_t numWarpsPerBucket[rebuff->numBuckets];

	//Init buckets (32 buckets => max 4096 bases)
	for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
		numCandidatesPerBucket[idBucket] = 0;
		numWarpsPerBucket[idBucket] = 0;
		tmpBuckets[idBucket] = 0;
	}

	//Fill buckets with elements per bucket
	for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
		idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].size - 1) / BASES_PER_THREAD;
		idBucket = (idBucket < (rebuff->numBuckets - 1)) ? idBucket : (rebuff->numBuckets - 1);
		numCandidatesPerBucket[idBucket]++;
	}

	//Number of warps per bucket
	rebuff->candidatesPerBufffer = 0;
	for(idBucket = 0; idBucket < rebuff->numBuckets - 1; idBucket++){
		numThreadsPerQuery = idBucket + 1;
		numQueriesPerWarp = SIZE_WARP / numThreadsPerQuery;
		numWarpsPerBucket[idBucket] = (numCandidatesPerBucket[idBucket] / numQueriesPerWarp) +
											((numCandidatesPerBucket[idBucket] % numQueriesPerWarp) ? 1 : 0);
		rebuff->h_initPosPerBucket[idBucket] = rebuff->candidatesPerBufffer;
		rebuff->candidatesPerBufffer += numWarpsPerBucket[idBucket] * numQueriesPerWarp;
	}

	//Fill init position warps per bucket
	for(idBucket = 1; idBucket < rebuff->numBuckets; idBucket++)
		rebuff->h_initWarpPerBucket[idBucket] = rebuff->h_initWarpPerBucket[idBucket-1] + numWarpsPerBucket[idBucket-1];

	//Allocate buffer (candidates)
	for(idBuff = 0; idBuff < rebuff->candidatesPerBufffer; idBuff++)
		rebuff->h_reorderBuffer[idBuff] = MAX_VALUE;

	//Set the number of real results in the reorder buffer
	res->numReorderedResults = rebuff->candidatesPerBufffer;

	//Reorder by size the candidates
	for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++)
		tmpBuckets[idBucket] = rebuff->h_initPosPerBucket[idBucket];
	for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
		idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].size - 1) / BASES_PER_THREAD;
		if (idBucket < (rebuff->numBuckets - 1)){
			rebuff->h_reorderBuffer[tmpBuckets[idBucket]] = idCandidate;
			tmpBuckets[idBucket]++;
		}
	}

	//Rellenar los huecos con elementos duplicados (se duplica el último)
	for(idBuff = 0; idBuff < rebuff->candidatesPerBufffer; idBuff++)
		if(rebuff->h_reorderBuffer[idBuff] == MAX_VALUE) rebuff->h_reorderBuffer[idBuff] = rebuff->h_reorderBuffer[idBuff-1];

	//Calculate the number of warps necessaries in the GPU
	for(idBucket = 0; idBucket < (rebuff->numBuckets - 1); idBucket++)
		rebuff->numWarps += numWarpsPerBucket[idBucket];

	return (SUCCESS);
}

myersError_t transferCPUtoGPU(buffer_t *mBuff)
{

	queries_buffer_t  	*qry  		= mBuff->queries;
	candidates_buffer_t	*cand 		= mBuff->candidates;
	reorder_buffer_t  	*rebuff 	= mBuff->reorderBuffer;
	results_buffer_t  	*res  		= mBuff->results;
	cudaStream_t 		idStream 	= mBuff->idStream;

	//Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->idDevice));

	//Transfer Binary Queries to GPU
	CUDA_ERROR(cudaMemcpyAsync(qry->d_queries, qry->h_queries, qry->totalQueriesEntries * sizeof(qryEntry_t), cudaMemcpyHostToDevice, idStream));

	//Transfer to GPU the information associated with Binary Queries
	CUDA_ERROR(cudaMemcpyAsync(qry->d_qinfo, qry->h_qinfo, qry->numQueries * sizeof(qryInfo_t), cudaMemcpyHostToDevice, idStream));

	//Transfer Candidates to GPU
	CUDA_ERROR(cudaMemcpyAsync(cand->d_candidates, cand->h_candidates, cand->numCandidates * sizeof(candInfo_t), cudaMemcpyHostToDevice, idStream));

	//Transfer reordered buffer to GPU
	CUDA_ERROR(cudaMemcpyAsync(rebuff->d_reorderBuffer, rebuff->h_reorderBuffer, rebuff->candidatesPerBufffer * sizeof(uint32_t), cudaMemcpyHostToDevice, idStream));
	//Transfer bucket information to GPU
	CUDA_ERROR(cudaMemcpyAsync(rebuff->d_initPosPerBucket, rebuff->h_initPosPerBucket, rebuff->numBuckets * sizeof(uint32_t), cudaMemcpyHostToDevice, idStream));
	CUDA_ERROR(cudaMemcpyAsync(rebuff->d_initWarpPerBucket, rebuff->h_initWarpPerBucket, rebuff->numBuckets * sizeof(uint32_t), cudaMemcpyHostToDevice, idStream));

	return (SUCCESS);
}

myersError_t transferGPUtoCPU(buffer_t *mBuff)
{
	cudaStream_t 		idStream 	= mBuff->idStream;
	results_buffer_t  	*res  		= mBuff->results;

	//Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->idDevice));
	CUDA_ERROR(cudaMemcpyAsync(res->h_reorderResults, res->d_reorderResults, res->numReorderedResults * sizeof(resEntry_t), cudaMemcpyDeviceToHost, idStream));

	return (SUCCESS);
}

void sendMyersBuffer(void *myersBuffer, float distance, uint32_t numPEQEntries, uint32_t numQueries, uint32_t numCandidates)
{
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	uint32_t error;

	//Set real the size of the things
	mBuff->queries->distance 			= distance;
	mBuff->queries->totalQueriesEntries = numPEQEntries;
	mBuff->queries->numQueries 			= numQueries;
	mBuff->candidates->numCandidates 	= numCandidates;
	mBuff->results->numResults			= numCandidates;

	queries_buffer_t  	*qry  		= mBuff->queries;
	candidates_buffer_t	*cand 		= mBuff->candidates;
	reorder_buffer_t  	*rebuff 	= mBuff->reorderBuffer;
	results_buffer_t  	*res  		= mBuff->results;
	cudaStream_t 		idStream 	= mBuff->idStream;

	//Process Kernel in Asynchonous way
	MYERS_ERROR(reorderingBuffer(mBuff));
	MYERS_ERROR(transferCPUtoGPU(mBuff));

	if( mBuff->majorCC <  2) MYERS_ERROR(E_INCOMPATIBLE_GPU);
	if( mBuff->majorCC == 2) MYERS_ERROR(processMyersBufferOnFermi(mBuff));
	if((mBuff->majorCC == 3) && (mBuff->minorCC == 0)) MYERS_ERROR(processMyersBufferOnKepler1stGen(mBuff));
	if((mBuff->majorCC >= 3) && (mBuff->minorCC >= 5)) MYERS_ERROR(processMyersBufferOnKepler2ndGen(mBuff));
	if((mBuff->majorCC >= 5) && (mBuff->minorCC >= 0)) MYERS_ERROR(processMyersBufferOnMaxwell1stGen(mBuff));
	MYERS_ERROR(transferGPUtoCPU(mBuff));
}

/************************************************************
Funtions to recive & process a Myers buffer from GPU
************************************************************/

myersError_t reorderingResults(buffer_t *mBuff)
{
	reorder_buffer_t  	*rebuff = mBuff->reorderBuffer;
	results_buffer_t  	*res  	= mBuff->results;

	uint32_t idRes;

	for(idRes = 0; idRes < res->numReorderedResults; idRes++){
		res->h_results[rebuff->h_reorderBuffer[idRes]].column = res->h_reorderResults[idRes].column;
		res->h_results[rebuff->h_reorderBuffer[idRes]].score = res->h_reorderResults[idRes].score;
	}

	return (SUCCESS);
}

void receiveMyersBuffer(void *myersBuffer)
{
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	uint32_t error;

	//Synchronize Stream (the thread wait for the commands done in the stream)
	CUDA_ERROR(cudaStreamSynchronize(mBuff->idStream));
	//Reordening the final results
	MYERS_ERROR(reorderingResults(mBuff));
}



/************************************************************
Funtions to free all the Myers resources
************************************************************/

myersError_t freeReference(reference_buffer_t *reference)
{
    if(reference->d_reference != NULL){
        CUDA_ERROR(cudaFree(reference->d_reference));
        reference->d_reference = NULL;
    }

    /* TODO: Liberar referencia almacenada en forma streaming */


    if(reference != NULL){
        free(reference);
        reference = NULL;
    }

    return(SUCCESS);
}

myersError_t freeQueries(queries_buffer_t *queries)
{
    if(queries->h_queries != NULL){
        CUDA_ERROR(cudaFreeHost(queries->h_queries));
        queries->h_queries = NULL;
    }

    if(queries->d_queries != NULL){
        CUDA_ERROR(cudaFree(queries->d_queries));
        queries->d_queries = NULL;
    }

    if(queries->h_qinfo != NULL){
        CUDA_ERROR(cudaFreeHost(queries->h_qinfo));
        queries->h_qinfo = NULL;
    }

    if(queries->d_qinfo != NULL){
        CUDA_ERROR(cudaFree(queries->d_qinfo));
        queries->d_qinfo = NULL;
    }

    if(queries != NULL){
        free(queries);
        queries = NULL;
    }

    return(SUCCESS);
}

myersError_t freeResults(results_buffer_t *results)
{

    if(results->h_results != NULL){
        free(results->h_results);
        results->h_results = NULL;
    }

    if(results->h_reorderResults != NULL){
        CUDA_ERROR(cudaFreeHost(results->h_reorderResults));
        results->h_reorderResults = NULL;
    }

    if(results->d_reorderResults != NULL){
        CUDA_ERROR(cudaFree(results->d_reorderResults));
        results->d_reorderResults = NULL;
    }

    if(results != NULL){
        free(results);
        results = NULL;
    }

    return(SUCCESS);
}

myersError_t freeReorderBuffer(reorder_buffer_t *reorderBuffer)
{
	if(reorderBuffer->h_reorderBuffer != NULL){
		CUDA_ERROR(cudaFreeHost(reorderBuffer->h_reorderBuffer));
		reorderBuffer->h_reorderBuffer = NULL;
	}

	if(reorderBuffer->h_initPosPerBucket != NULL){
		CUDA_ERROR(cudaFreeHost(reorderBuffer->h_initPosPerBucket));
		reorderBuffer->h_initPosPerBucket = NULL;
	}

	if(reorderBuffer->h_initWarpPerBucket != NULL){
		CUDA_ERROR(cudaFreeHost(reorderBuffer->h_initWarpPerBucket));
		reorderBuffer->h_initWarpPerBucket = NULL;
	}

	if(reorderBuffer != NULL){
        free(reorderBuffer);
        reorderBuffer = NULL;
    }

	return(SUCCESS);
}

myersError_t freeCandidates(candidates_buffer_t *candidates)
{
	if(candidates->h_candidates != NULL){
        CUDA_ERROR(cudaFreeHost(candidates->h_candidates));
        candidates->h_candidates = NULL;
    }

    if(candidates->d_candidates != NULL){
        CUDA_ERROR(cudaFree(candidates->d_candidates));
        candidates->d_candidates = NULL;
    }

	if(candidates != NULL){
        free(candidates);
        candidates = NULL;
    }

	return(SUCCESS);
}

void endMyers(void ***myersBuffer)
{
	buffer_t **mBuff = (buffer_t **) (* myersBuffer);
	uint32_t numBuffers = mBuff[0]->numBuffers;
	uint32_t idBuffer;

	//Synchronize all the jobs with the GPU
	CUDA_ERROR(cudaDeviceSynchronize());

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		MYERS_ERROR(freeReference(mBuff[idBuffer]->reference));
		MYERS_ERROR(freeQueries(mBuff[idBuffer]->queries));
		MYERS_ERROR(freeCandidates(mBuff[idBuffer]->candidates));
		MYERS_ERROR(freeReorderBuffer(mBuff[idBuffer]->reorderBuffer));
		MYERS_ERROR(freeResults(mBuff[idBuffer]->results));

		MYERS_ERROR(cudaStreamDestroy(mBuff[idBuffer]->idStream));
	}

	if(mBuff != NULL){
    	free(mBuff);
    	mBuff = NULL;
    }

	(* myersBuffer) = (void **) mBuff;
}
