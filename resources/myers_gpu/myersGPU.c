#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <string.h>
#include "myers-common.h"

/************************************************************
Functions to handle errors
************************************************************/

MYERS_INLINE void CudaError( cudaError_t err, const char *file,  int line ) {
   	if (err != cudaSuccess) {
      		fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err),  file, line );
       		exit(EXIT_FAILURE);
   	}
}

MYERS_INLINE char *myersGetErrorString(myersError_t error){
    switch(error) {
        case E_OPENING_FILE:  			return "Myers - Error: opening file"; break;
        case E_READING_FILE:  			return "Myers - Error: reading file"; break;
        case E_INSUFFICIENT_MEM_GPU:	return "Myers - Error: there aren't enough GPU memory space"; break;
        case E_ALLOCATE_MEM: 			return "Myers - Error: allocating data"; break;
        case E_INCOMPATIBLE_GPU:		return "Myers - Error: incompatible GPU (old CC version)"; break;
        case E_REFERENCE_CODING:		return "Myers - Error: reference coding not supported"; break;
        default: 						return "Myers - Unknown error";
    }
}

MYERS_INLINE void MyersError(myersError_t err, const char *file,  int line ) {
   	if (err != 0) {
      		fprintf(stderr, "%s in %s at line %d\n", myersGetErrorString(err),  file, line );
       		exit(EXIT_FAILURE);
   	}
}



/************************************************************
Functions to get the Myers buffers
************************************************************/

MYERS_INLINE uint32_t getMaxPEQEntries(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->maxPEQEntries);
}

MYERS_INLINE uint32_t getMaxCandidates(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->maxCandidates);
}

MYERS_INLINE uint32_t getMaxQueries(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->maxQueries);
}

MYERS_INLINE qryEntry_t* getPEQBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->queries->h_queries);
}

MYERS_INLINE candInfo_t* getCandidatesBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->candidates->h_candidates);
}

MYERS_INLINE qryInfo_t* getPEQInfoBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->queries->h_qinfo);
}

MYERS_INLINE resEntry_t* getResultsBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->results->h_results);
}

MYERS_INLINE uint32_t getIdDeviceBuffer(void *myersBuffer){
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	return(mBuff->device->idDevice);
}


/************************************************************
Functions to init all the Myers resources
************************************************************/

uint32_t getDeviceFreeMemory(uint32_t idDevice)
{
	size_t free, total;
    CUDA_ERROR(cudaSetDevice(idDevice));
	CUDA_ERROR(cudaMemGetInfo(&free, &total));

	return (CONVERT_B_TO_MB(free));
}

MYERS_INLINE myersError_t freeDevicesListHost(device_info_t ***deviceList)
{
	device_info_t **device = (* deviceList);

    if(device != NULL){
        free(device);
        device = NULL;
    }

    (* deviceList) = device;
    return(SUCCESS);
}

MYERS_INLINE myersError_t freeReferenceHost(reference_buffer_t *reference)
{
    if(reference->h_reference != NULL){
        CUDA_ERROR(cudaFreeHost(reference->h_reference));
        reference->h_reference = NULL;
    }

    return(SUCCESS);
}

MYERS_INLINE myersError_t freeReferenceDevice(reference_buffer_t *reference, device_info_t **devices)
{
	uint32_t idSupportedDevice, numSupportedDevices;

	numSupportedDevices = devices[0]->numSupportedDevices;

	//Free all the references in the devices
    for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
	    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
        if(reference->d_reference[idSupportedDevice] != NULL){
            CUDA_ERROR(cudaFree(reference->d_reference[idSupportedDevice]));
            reference->d_reference[idSupportedDevice] = NULL;
        }
    }

    //Free the reference list
    if(reference->d_reference != NULL){
        free(reference->d_reference);
        reference->d_reference = NULL;
    }

    return(SUCCESS);
}

MYERS_INLINE uint32_t charToBinASCII(unsigned char base)
{
	switch(base)
	{
    	case 'A':
    	case 'a':
    	    return(ENC_DNA_CHAR_A);
    	case 'C':
    	case 'c':
    	    return(ENC_DNA_CHAR_C << (UINT32_LENGTH - REFERENCE_CHAR_LENGTH));
    	case 'G':
    	case 'g':
    	    return(ENC_DNA_CHAR_G << (UINT32_LENGTH - REFERENCE_CHAR_LENGTH));
    	case 'T':
    	case 't':
    	    return(ENC_DNA_CHAR_T << (UINT32_LENGTH - REFERENCE_CHAR_LENGTH));
    	default :
    	    return(ENC_DNA_CHAR_N << (UINT32_LENGTH - REFERENCE_CHAR_LENGTH));
	}
}

MYERS_INLINE myersError_t transformReferenceASCII(const char *referenceASCII, reference_buffer_t *reference)
{
	uint32_t indexBase, bitmap;
	uint64_t idEntry, i, referencePosition;
	unsigned char referenceChar;

	CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint32_t), cudaHostAllocMapped));

	for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
		bitmap = 0;
		for(i = 0; i < REFERENCE_CHARS_PER_ENTRY; i++){
			referencePosition = idEntry * REFERENCE_CHARS_PER_ENTRY + i;
			if (referencePosition < reference->size) referenceChar = referenceASCII[referencePosition];
				else referenceChar = 'N'; //filling reference padding
			indexBase = charToBinASCII(referenceChar);
			bitmap = (bitmap >> REFERENCE_CHAR_LENGTH) | indexBase;
		}
		reference->h_reference[referencePosition / REFERENCE_CHARS_PER_ENTRY] = bitmap;
	}
	return(SUCCESS);
}

MYERS_INLINE myersError_t loadReferenceMFASTA(const char *fn, void *reference)
{
	reference_buffer_t *ref = (reference_buffer_t *) reference;
	FILE *fp = NULL;
	char lineFile[FILE_SIZE_LINES], *tmp_reference;
	uint64_t sizeFile = 0, position = 0;
	int32_t charsRead = 0;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (E_OPENING_FILE);

	fseek(fp, 0L, SEEK_END);
	sizeFile = ftell(fp);
	rewind(fp);

	tmp_reference = (char*) malloc(sizeFile * sizeof(char));
	if (ref == NULL) return (E_ALLOCATE_MEM);

	if ((fgets(lineFile, FILE_SIZE_LINES, fp) == NULL) || (lineFile[0] != '>'))
		return (E_READING_FILE);

	while((!feof(fp)) && (fgets(lineFile, FILE_SIZE_LINES, fp) != NULL)){
		if (lineFile[0] != '>'){
			charsRead = strlen(lineFile);
			if(charsRead) charsRead--;
			memcpy((tmp_reference + position), lineFile, charsRead);
			position +=  charsRead;
		}
	}

	ref->size = position;
	ref->numEntries = DIV_CEIL(ref->size, REFERENCE_CHARS_PER_ENTRY) + REFERENCE_END_PADDING;

	MYERS_ERROR(transformReferenceASCII(tmp_reference, ref));

	fclose(fp);
	free(tmp_reference);
	return (SUCCESS);
}

MYERS_INLINE myersError_t loadReferencePROFILE(const char *fn, reference_buffer_t *reference)
{
	FILE *fp = NULL;
	size_t result;

	fp = fopen(fn, "rb");
		if (fp == NULL) return (E_OPENING_FILE);

    result = fread(&reference->numEntries, sizeof(uint64_t), 1, fp);
		if (result != 1) return (E_READING_FILE);
    result = fread(&reference->size, sizeof(uint64_t), 1, fp);
		if (result != 1) return (E_READING_FILE);

	CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint32_t), cudaHostAllocMapped));

	result = fread(reference->h_reference, sizeof(uint32_t), reference->numEntries, fp);
		if (result != reference->numEntries) return (E_READING_FILE);

	fclose(fp);
	return (SUCCESS);
}

MYERS_INLINE myersError_t transferReferenceCPUtoGPUs(reference_buffer_t *reference, device_info_t **devices)
{
	uint32_t deviceFreeMemory, idSupportedDevice;
	uint32_t numSupportedDevices = devices[0]->numSupportedDevices;

	printf("Loading Reference in devices .... \n");

	reference->d_reference = (uint32_t **) malloc(numSupportedDevices * sizeof(uint32_t *));
	if (reference->d_reference == NULL) MYERS_ERROR(E_ALLOCATE_MEM);

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		deviceFreeMemory = getDeviceFreeMemory(devices[idSupportedDevice]->idDevice);
		if ((CONVERT_B_TO_MB(reference->numEntries * sizeof(uint32_t))) > deviceFreeMemory) return(E_INSUFFICIENT_MEM_GPU);

	    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
		//Synchronous allocate & transfer Binary Reference to GPU
		CUDA_ERROR(cudaMalloc((void**) &reference->d_reference[idSupportedDevice], reference->numEntries * sizeof(uint32_t)));
		CUDA_ERROR(cudaMemcpy(reference->d_reference[idSupportedDevice], reference->h_reference, reference->numEntries * sizeof(uint32_t), cudaMemcpyHostToDevice));
	}

	return (SUCCESS);
}

MYERS_INLINE myersError_t transformReferenceGEM(const char *referenceGEM, reference_buffer_t *reference)
{
	uint32_t indexBase, bitmap;
	uint64_t idEntry, i, referencePosition;
	unsigned char referenceChar;

	CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint32_t), cudaHostAllocMapped));

	for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
		bitmap = 0;
		for(i = 0; i < REFERENCE_CHARS_PER_ENTRY; i++){
			referencePosition = idEntry * REFERENCE_CHARS_PER_ENTRY + i;
			if (referencePosition < reference->size) referenceChar = referenceGEM[referencePosition];
				else referenceChar = 'N'; //filling reference padding
			indexBase = ((uint32_t) referenceChar) << (UINT32_LENGTH - REFERENCE_CHAR_LENGTH);
			bitmap = (bitmap >> REFERENCE_CHAR_LENGTH) | indexBase;
		}
		reference->h_reference[referencePosition / REFERENCE_CHARS_PER_ENTRY] = bitmap;
	}
	return(SUCCESS);
}

MYERS_INLINE myersError_t initReference(reference_buffer_t **reference, const char *referenceRaw, uint64_t refSize, refCoding_t refCoding)
{
	reference_buffer_t *ref = (reference_buffer_t *) malloc(sizeof(reference_buffer_t));

	ref->d_reference = NULL;
	ref->h_reference = NULL;
	ref->size = refSize;
	ref->numEntries = DIV_CEIL(ref->size, REFERENCE_CHARS_PER_ENTRY) + REFERENCE_END_PADDING;

	switch(refCoding){
    	case ASCII:						MYERS_ERROR(transformReferenceASCII(referenceRaw, ref));
    		break;
    	case GEM:						MYERS_ERROR(transformReferenceGEM(referenceRaw, ref));
    		break;
    	case MFASTA_FILE:				MYERS_ERROR(loadReferenceMFASTA(referenceRaw, ref));
    		break;
    	case PROFILE_REFERENCE_FILE:	MYERS_ERROR(loadReferencePROFILE(referenceRaw, ref));
    		break;
    	default:						MYERS_ERROR(E_REFERENCE_CODING);
	}

	(* reference) = ref;
	return (SUCCESS);
}

MYERS_INLINE float sizePerCandidate(uint32_t averageNumPEQEntries, uint32_t candidatesPerQuery)
{
	uint32_t bytesPerQuery = averageNumPEQEntries * sizeof(qryEntry_t) + sizeof(qryInfo_t);
	uint32_t bytesCandidate = sizeof(candInfo_t);
	uint32_t bytesResult = sizeof(candInfo_t);
	uint32_t bytesBiningProcess = sizeof(candInfo_t);

	return((bytesPerQuery/(float)candidatesPerQuery) + bytesCandidate
			+ bytesResult + bytesBiningProcess);
}

MYERS_INLINE uint32_t candidatesPerBinningPadding()
{
	uint32_t idBucket;
	uint32_t bucketPaddingCandidates = 0;

	/* Worst number of dummy candidates added for padding the binning */
	for(idBucket = 1; idBucket < NUM_BUCKETS_FOR_BINNING-1; ++idBucket)
		bucketPaddingCandidates += (WARP_SIZE / idBucket);

	/* Increased bytes per buffer taking account the padding*/
	return(bucketPaddingCandidates);
}

MYERS_INLINE myersError_t initQueries(queries_buffer_t *queries, uint32_t maxPEQEntries, uint32_t maxQueries)
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

MYERS_INLINE myersError_t initCandidates(candidates_buffer_t *candidates, uint32_t maxCandidates)
{
	candidates->numCandidates = 0;
	candidates->h_candidates = NULL;
	candidates->d_candidates = NULL;

	//Allocate candidates info in CPU & GPU
	CUDA_ERROR(cudaHostAlloc((void**) &candidates->h_candidates, maxCandidates * sizeof(candInfo_t), cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &candidates->d_candidates, maxCandidates * sizeof(candInfo_t)));

	return (SUCCESS);
}

MYERS_INLINE myersError_t initReorderBuffer(reorder_buffer_t *reorderBuffer, uint32_t maxReorderBuffer)
{
	reorderBuffer->numBuckets = NUM_BUCKETS_FOR_BINNING;
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

MYERS_INLINE myersError_t initResults(results_buffer_t *results, uint32_t maxReorderBuffer, uint32_t maxCandidates)
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

MYERS_INLINE myersError_t initBuffer(buffer_t *buffer, uint32_t idBuffer, device_info_t *device, uint32_t numBuffers,
			   uint32_t maxCandidates, uint32_t maxQueries, uint32_t maxPEQEntries,
			   uint32_t bucketPaddingCandidates, reference_buffer_t *reference)
{
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
    buffer->device			=	device;


	buffer->queries			=	(queries_buffer_t *)    malloc(sizeof(queries_buffer_t));
		if (buffer->queries == NULL) return (E_ALLOCATE_MEM);
	buffer->candidates		=	(candidates_buffer_t *) malloc(sizeof(candidates_buffer_t));
		if (buffer->candidates == NULL) return (E_ALLOCATE_MEM);
	buffer->reorderBuffer	=	(reorder_buffer_t *)    malloc(sizeof(reorder_buffer_t));
		if (buffer->reorderBuffer == NULL) return (E_ALLOCATE_MEM);
	buffer->results 		=	(results_buffer_t *)    malloc(sizeof(results_buffer_t));
		if (buffer->results == NULL) return (E_ALLOCATE_MEM);

	//Set in which Device we create and initialize the structures
    CUDA_ERROR(cudaSetDevice(buffer->device->idDevice));

	//Create the CUDA stream per each buffer
	CUDA_ERROR(cudaStreamCreate(&buffer->idStream));

	MYERS_ERROR(initQueries(buffer->queries, buffer->maxPEQEntries, buffer->maxQueries));
	MYERS_ERROR(initCandidates(buffer->candidates, buffer->maxCandidates));
	MYERS_ERROR(initReorderBuffer(buffer->reorderBuffer, buffer->maxReorderBuffer));
	MYERS_ERROR(initResults(buffer->results, buffer->maxReorderBuffer, buffer->maxCandidates));

	return(SUCCESS);
}

MYERS_INLINE myersError_t initDeviceBuffers(buffer_t ***myersBuffer, uint32_t numBuffers, reference_buffer_t *reference,
											device_info_t **device, uint32_t maxMbPerBuffer, uint32_t averageQuerySize,
											uint32_t candidatesPerQuery)
{
	uint32_t numSupportedDevices, idSupportedDevice, idDevice, idGlobalBuffer, numBuffersPerDevice, idLocalBuffer;
	uint32_t freeDeviceMemory, bucketPaddingCandidates, averageNumPEQEntries, maxCandidates, maxPEQEntries, maxQueries, maxMbPerDevice;
	size_t bytesPadding, bytesPerBuffer;
	float bytesPerCandidate;
	int32_t remainderBuffers;

	buffer_t **buffer = (buffer_t **) malloc(numBuffers * sizeof(buffer_t *));
		if (buffer == NULL) MYERS_ERROR(E_ALLOCATE_MEM);

	printf("Allocating device buffers .... \n");

	numSupportedDevices = device[0]->numSupportedDevices;
	remainderBuffers = numBuffers;
	idGlobalBuffer = 0;

	averageNumPEQEntries = DIV_CEIL(averageQuerySize, PEQ_ENTRY_LENGTH);
	bucketPaddingCandidates = candidatesPerBinningPadding();

	bytesPadding = bucketPaddingCandidates * (sizeof(uint32_t) + sizeof(resEntry_t));
	bytesPerCandidate = sizePerCandidate(averageNumPEQEntries, candidatesPerQuery);

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		idDevice = device[idSupportedDevice]->idDevice;
		freeDeviceMemory = getDeviceFreeMemory(idDevice);

		numBuffersPerDevice = ROUND(numBuffers * device[idSupportedDevice]->relativePerformance);
		if(idSupportedDevice == numSupportedDevices-1) numBuffersPerDevice = remainderBuffers;

		if(maxMbPerBuffer != 0) maxMbPerDevice = MIN(numBuffersPerDevice * maxMbPerBuffer, freeDeviceMemory);
			else maxMbPerDevice = freeDeviceMemory;

		bytesPerBuffer = (CONVERT_MB_TO_B(maxMbPerDevice) - (CONVERT_MB_TO_B(maxMbPerDevice) * 0.05)) / numBuffersPerDevice;

		maxCandidates = (uint32_t) ((bytesPerBuffer - bytesPadding) / bytesPerCandidate);
		maxQueries    = maxCandidates / candidatesPerQuery;
		maxPEQEntries = maxQueries * averageNumPEQEntries;

		//TODO: Consider less buffers than GPUs!!
		printf("Requested: %d - Available: %d \n", CONVERT_B_TO_MB(bytesPerBuffer * numBuffersPerDevice), freeDeviceMemory);
		if ((maxCandidates < 1) || (maxQueries < 1) || (CONVERT_B_TO_MB(bytesPerBuffer * numBuffersPerDevice) > freeDeviceMemory))
			MYERS_ERROR(E_INSUFFICIENT_MEM_GPU);

		printf("Device %d: %d buffers x %d MBytes/buffer (max %d MB) = %d MBytes\n",
				idDevice, numBuffersPerDevice, CONVERT_B_TO_MB(bytesPerBuffer), maxMbPerBuffer,
				CONVERT_B_TO_MB(bytesPerBuffer) * numBuffersPerDevice);

		for(idLocalBuffer = 0; idLocalBuffer < numBuffersPerDevice; ++idLocalBuffer){
			buffer[idGlobalBuffer] = (buffer_t *) malloc(sizeof(buffer_t));
		    CUDA_ERROR(cudaSetDevice(idDevice));
			MYERS_ERROR(initBuffer(buffer[idGlobalBuffer], idGlobalBuffer, device[idSupportedDevice], numBuffers, maxCandidates,
						  maxQueries, maxPEQEntries, bucketPaddingCandidates, reference));
			idGlobalBuffer++;
		}
		remainderBuffers -= numBuffersPerDevice;
	  }

	(* myersBuffer) = buffer;
	return (SUCCESS);
}

MYERS_INLINE devArch_t getDeviceArchitecture(uint32_t idDevice)
{
	struct cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, idDevice);

	if (devProp.major <= 1) return(ARCH_TESLA);								/* CC 1.X		*/
	if (devProp.major == 2 && devProp.minor == 0) return(ARCH_FERMI_1G); 	/* CC 2.0		*/
	if (devProp.major == 2 && devProp.minor >  0) return(ARCH_FERMI_2G); 	/* CC 2.1		*/
	if (devProp.major == 3 && devProp.minor <  5) return(ARCH_KEPLER_1G); 	/* CC 3.0, 3.2	*/
	if (devProp.major == 3 && devProp.minor >= 5) return(ARCH_KEPLER_2G); 	/* CC 3.5		*/
	if (devProp.major == 5) return(ARCH_MAXWELL); 							/* CC 5.X		*/
	return(ARCH_NEWGEN);													/* CC X.X		*/
}

MYERS_INLINE uint32_t getSMCudaCores(devArch_t architecture)
{
	switch (architecture) {
		case ARCH_TESLA:		return(8);
		case ARCH_FERMI_1G:		return(32);
		case ARCH_FERMI_2G:		return(48);
		case ARCH_KEPLER_1G:	return(192);
		case ARCH_KEPLER_2G:	return(192);
		case ARCH_MAXWELL:		return(128);
		default:				return(128);
	}
}

MYERS_INLINE uint32_t getDeviceCudaCores(uint32_t idDevice)
{
	uint32_t coresPerSM;
	devArch_t architecture;

	struct cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, idDevice);

	architecture = getDeviceArchitecture(idDevice);
	coresPerSM = getSMCudaCores(architecture);

	return(devProp.multiProcessorCount * coresPerSM);
}

MYERS_INLINE myersError_t selectSupportedDevices(device_info_t ***devices, uint32_t minimumMemorySize, devArch_t selectedArchitectures)
{
	uint32_t idDevice, idSupportedDevice, numSupportedDevices, memoryFree, totalSystemPerformance = 0;
	int32_t numDevices;
	devArch_t architecture;

	CUDA_ERROR(cudaGetDeviceCount(&numDevices));
	printf("There are %d visible devices: \n", numDevices);
	printf("Compiled with CUDA SDK %d.%d \n", CUDA_VERSION/1000, CUDA_VERSION%1000);

	device_info_t **dev = (device_info_t **) malloc(numDevices * sizeof(device_info_t *));
	if (dev == NULL) MYERS_ERROR(E_ALLOCATE_MEM);

	for(idDevice = 0, idSupportedDevice = 0; idDevice < numDevices; ++idDevice){
		struct cudaDeviceProp devProp;
		cudaGetDeviceProperties(&devProp, idDevice);
		printf("Device %d: %s", idDevice, devProp.name);

		architecture = getDeviceArchitecture(idDevice);
		if((architecture & selectedArchitectures) == 0)
			printf("\t \t Compute Capability UNSUITABLE");

		memoryFree = getDeviceFreeMemory(idDevice); /* in MB */
		if(memoryFree < minimumMemorySize)
			printf("\t \t INSUFFICIENT DEVICE MEMORY (Mem Req: %u MBytes - Mem Avail: %u MBytes)", minimumMemorySize, memoryFree);

		printf("\n");

		if((architecture & selectedArchitectures) && (memoryFree > minimumMemorySize)){
			dev[idSupportedDevice] = NULL;
			dev[idSupportedDevice] = (device_info_t *) malloc(sizeof(device_info_t));
			if(dev == NULL) MYERS_ERROR(E_ALLOCATE_MEM);

			dev[idSupportedDevice]->architecture = architecture;
			dev[idSupportedDevice]->cudaCores = getDeviceCudaCores(idDevice);
			dev[idSupportedDevice]->frequency = devProp.clockRate;
			dev[idSupportedDevice]->idSupportedDevice = idSupportedDevice;
			dev[idSupportedDevice]->idDevice = idDevice;
			dev[idSupportedDevice]->numDevices = numDevices;
			dev[idSupportedDevice]->relativePerformance = dev[idSupportedDevice]->cudaCores * dev[idSupportedDevice]->frequency;
			totalSystemPerformance += dev[idSupportedDevice]->relativePerformance;
			idSupportedDevice++;
		}
	}

	numSupportedDevices = idSupportedDevice;
	if(numSupportedDevices == 0) return(E_NO_SUPPORTED_GPUS);
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		dev[idSupportedDevice]->relativePerformance /= totalSystemPerformance;
		dev[idSupportedDevice]->numSupportedDevices = numSupportedDevices;
	}

	(* devices) = dev;
	return(SUCCESS);
}

MYERS_INLINE uint32_t minMemorySizePerDevice(size_t *minimumMemorySize, reference_buffer_t *reference,
												uint32_t numBuffers, uint32_t averageQuerySize, uint32_t candidatesPerQuery)
{
	uint32_t idDevice, idBucket, numDevices;
	uint32_t averageNumPEQEntries, bucketPaddingCandidates, bytesPadding;
	uint32_t mBytesPerReference, mBytesPerBuffer; /* in MBytes */
	float bytesPerCandidate;

	averageNumPEQEntries = DIV_CEIL(averageQuerySize, PEQ_ENTRY_LENGTH);
	mBytesPerReference = CONVERT_B_TO_MB(reference->numEntries * UINT32_SIZE);
	bytesPerCandidate = sizePerCandidate(averageNumPEQEntries, candidatesPerQuery);

	/* Increased bytes per buffer taking account the padding*/
	bytesPadding = candidatesPerBinningPadding() * (sizeof(uint32_t) + sizeof(resEntry_t));

	mBytesPerBuffer =  CONVERT_B_TO_MB(MIN_CANDIDATES_PER_BUFFER * bytesPerCandidate + bytesPadding);
	(* minimumMemorySize) = mBytesPerReference + numBuffers * mBytesPerBuffer;

	return (SUCCESS);
}

MYERS_INLINE void initMyers(void ***myersBuffer, uint32_t numBuffers, uint32_t maxMbPerBuffer,
							const char *referenceRaw, refCoding_t refCoding, const uint64_t refSize,
							uint32_t averageQuerySize, uint32_t candidatesPerQuery, devArch_t selectedArchitectures)
{
  buffer_t				**buffer = NULL;
  reference_buffer_t 	*reference = NULL;
  device_info_t			**devices = NULL;
  size_t				minimumMemorySize = 0; /*in MBytes */

  MYERS_ERROR(initReference(&reference, referenceRaw, refSize, refCoding));
  MYERS_ERROR(minMemorySizePerDevice(&minimumMemorySize, reference, numBuffers, averageQuerySize, candidatesPerQuery));
  MYERS_ERROR(selectSupportedDevices(&devices, minimumMemorySize, selectedArchitectures));

  MYERS_ERROR(transferReferenceCPUtoGPUs(reference, devices));
  MYERS_ERROR(initDeviceBuffers(&buffer, numBuffers, reference, devices, maxMbPerBuffer, averageQuerySize, candidatesPerQuery));

  MYERS_ERROR(freeReferenceHost(reference));
  MYERS_ERROR(freeDevicesListHost(&devices));

  (* myersBuffer) = (void **) buffer;
}

/************************************************************
Functions to send & process a Myers buffer to GPU
************************************************************/

MYERS_INLINE myersError_t reorderingBuffer(buffer_t *mBuff)
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

	//Re-init the reorderBuffer (to reuse the buffer)
	rebuff->numBuckets = NUM_BUCKETS_FOR_BINNING;
	rebuff->candidatesPerBufffer = 0;
	rebuff->numWarps = 0;

	//Init buckets (32 buckets => max 4096 bases)
	for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++){
		numCandidatesPerBucket[idBucket] = 0;
		numWarpsPerBucket[idBucket] = 0;
		tmpBuckets[idBucket] = 0;
	}

	//Fill buckets with elements per bucket
	for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
		idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].size - 1) / PEQ_LENGTH_PER_CUDA_THREAD;
		idBucket = (idBucket < (rebuff->numBuckets - 1)) ? idBucket : (rebuff->numBuckets - 1);
		numCandidatesPerBucket[idBucket]++;
	}

	//Number of warps per bucket
	rebuff->candidatesPerBufffer = 0;
	for(idBucket = 0; idBucket < rebuff->numBuckets - 1; idBucket++){
		numThreadsPerQuery = idBucket + 1;
		numQueriesPerWarp = WARP_SIZE / numThreadsPerQuery;
		numWarpsPerBucket[idBucket] = DIV_CEIL(numCandidatesPerBucket[idBucket], numQueriesPerWarp);
		rebuff->h_initPosPerBucket[idBucket] = rebuff->candidatesPerBufffer;
		rebuff->candidatesPerBufffer += numWarpsPerBucket[idBucket] * numQueriesPerWarp;
	}

	//Fill the start position warps for each bucket
	for(idBucket = 1; idBucket < rebuff->numBuckets; idBucket++)
		rebuff->h_initWarpPerBucket[idBucket] = rebuff->h_initWarpPerBucket[idBucket-1] + numWarpsPerBucket[idBucket-1];

	//Allocate buffer (candidates)
	for(idBuff = 0; idBuff < rebuff->candidatesPerBufffer; idBuff++)
		rebuff->h_reorderBuffer[idBuff] = UINT32_ONES;

	//Set the number of real results in the reorder buffer
	res->numReorderedResults = rebuff->candidatesPerBufffer;

	//Reorder by size the candidates
	for(idBucket = 0; idBucket < rebuff->numBuckets; idBucket++)
		tmpBuckets[idBucket] = rebuff->h_initPosPerBucket[idBucket];
	for(idCandidate = 0; idCandidate < cand->numCandidates; idCandidate++){
		idBucket = (qry->h_qinfo[cand->h_candidates[idCandidate].query].size - 1) / PEQ_LENGTH_PER_CUDA_THREAD;
		if (idBucket < (rebuff->numBuckets - 1)){
			rebuff->h_reorderBuffer[tmpBuckets[idBucket]] = idCandidate;
			tmpBuckets[idBucket]++;
		}
	}

	//Fill the paddings with replicated candidates (always the last candidate)
	for(idBuff = 0; idBuff < rebuff->candidatesPerBufffer; idBuff++)
		if(rebuff->h_reorderBuffer[idBuff] == UINT32_ONES) rebuff->h_reorderBuffer[idBuff] = rebuff->h_reorderBuffer[idBuff-1];

	//Calculate the number of warps necessaries in the GPU
	for(idBucket = 0; idBucket < (rebuff->numBuckets - 1); idBucket++)
		rebuff->numWarps += numWarpsPerBucket[idBucket];

	return (SUCCESS);
}

MYERS_INLINE myersError_t transferCPUtoGPU(buffer_t *mBuff)
{

	queries_buffer_t  	*qry  		= mBuff->queries;
	candidates_buffer_t	*cand 		= mBuff->candidates;
	reorder_buffer_t  	*rebuff 	= mBuff->reorderBuffer;
	results_buffer_t  	*res  		= mBuff->results;
	cudaStream_t 		idStream 	= mBuff->idStream;

	//Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device->idDevice));

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

MYERS_INLINE myersError_t transferGPUtoCPU(buffer_t *mBuff)
{
	cudaStream_t 		idStream 	= mBuff->idStream;
	results_buffer_t  	*res  		= mBuff->results;

	//Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device->idDevice));
	CUDA_ERROR(cudaMemcpyAsync(res->h_reorderResults, res->d_reorderResults, res->numReorderedResults * sizeof(resEntry_t), cudaMemcpyDeviceToHost, idStream));

	return (SUCCESS);
}

MYERS_INLINE void sendMyersBuffer(void *myersBuffer, uint32_t numPEQEntries, uint32_t numQueries, uint32_t numCandidates)
{
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	uint32_t error;

	//Set real size of the things
	mBuff->queries->totalQueriesEntries = numPEQEntries;
	mBuff->queries->numQueries 			= numQueries;
	mBuff->candidates->numCandidates 	= numCandidates;
	mBuff->results->numResults			= numCandidates;

	queries_buffer_t  	*qry  		= mBuff->queries;
	candidates_buffer_t	*cand 		= mBuff->candidates;
	reorder_buffer_t  	*rebuff 	= mBuff->reorderBuffer;
	results_buffer_t  	*res  		= mBuff->results;
	cudaStream_t 		idStream 	= mBuff->idStream;

	//Process Kernel in Asynchronous way
	MYERS_ERROR(reorderingBuffer(mBuff));
	MYERS_ERROR(transferCPUtoGPU(mBuff));

	#if (CUDA_VERSION >= 3000)
		if(mBuff->device->architecture & (ARCH_FERMI_1G | ARCH_FERMI_2G)) MYERS_ERROR(processMyersBufferOnFermi(mBuff));
	#endif

	#if (CUDA_VERSION >= 5000)
		if(mBuff->device->architecture & ARCH_KEPLER_1G) MYERS_ERROR(processMyersBufferOnKepler1stGen(mBuff));
		if(mBuff->device->architecture & ARCH_KEPLER_2G) MYERS_ERROR(processMyersBufferOnKepler2ndGen(mBuff));
	#endif

	#if (CUDA_VERSION >= 6000)
		if(mBuff->device->architecture & (ARCH_MAXWELL | ARCH_NEWGEN)) MYERS_ERROR(processMyersBufferOnMaxwell1stGen(mBuff));
		/* INCLUDED SUPPORT for future GPUs with PTX ASM code (JIT compiling) */
	#endif

	MYERS_ERROR(transferGPUtoCPU(mBuff));
}

/************************************************************
Functions to receive & process a Myers buffer from GPU
************************************************************/

MYERS_INLINE myersError_t reorderingResults(buffer_t *mBuff)
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

MYERS_INLINE void receiveMyersBuffer(void *myersBuffer)
{
	buffer_t *mBuff = (buffer_t *) myersBuffer;
	uint32_t error;

	//Synchronize Stream (the thread wait for the commands done in the stream)
	CUDA_ERROR(cudaStreamSynchronize(mBuff->idStream));
	//Reorder the final results
	MYERS_ERROR(reorderingResults(mBuff));
}



/************************************************************
Functions to free all the Myers resources
************************************************************/

MYERS_INLINE myersError_t freeReference(reference_buffer_t **referenceBuffer, device_info_t **devices)
{
	reference_buffer_t *reference = (* referenceBuffer);

    MYERS_ERROR(freeReferenceHost(reference));
    MYERS_ERROR(freeReferenceDevice(reference, devices));

    if(reference != NULL){
        free(reference);
        reference = NULL;
    }

    (* referenceBuffer) = reference;
    return(SUCCESS);
}

MYERS_INLINE myersError_t freeQueries(queries_buffer_t **queriesBuffer)
{
	queries_buffer_t *queries = (* queriesBuffer);

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

    (* queriesBuffer) = queries;
    return(SUCCESS);
}

MYERS_INLINE myersError_t freeResults(results_buffer_t **resultsBuffer)
{
	results_buffer_t *results = (* resultsBuffer);

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

    (* resultsBuffer) = results;
    return(SUCCESS);
}

MYERS_INLINE myersError_t freeReorderBuffer(reorder_buffer_t **resultsBuffer)
{
	reorder_buffer_t *reorder = (* resultsBuffer);

	if(reorder->h_reorderBuffer != NULL){
		CUDA_ERROR(cudaFreeHost(reorder->h_reorderBuffer));
		reorder->h_reorderBuffer = NULL;
	}

	if(reorder->h_initPosPerBucket != NULL){
		CUDA_ERROR(cudaFreeHost(reorder->h_initPosPerBucket));
		reorder->h_initPosPerBucket = NULL;
	}

	if(reorder->h_initWarpPerBucket != NULL){
		CUDA_ERROR(cudaFreeHost(reorder->h_initWarpPerBucket));
		reorder->h_initWarpPerBucket = NULL;
	}

	if(reorder != NULL){
        free(reorder);
        reorder = NULL;
    }

	(* resultsBuffer) = reorder;
	return(SUCCESS);
}

MYERS_INLINE myersError_t freeCandidates(candidates_buffer_t **candidatesBuffer)
{
	candidates_buffer_t *candidates = (* candidatesBuffer);

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

	(* candidatesBuffer) = candidates;
	return(SUCCESS);
}

MYERS_INLINE myersError_t freeDevicesInfo(device_info_t ***deviceList)
{
	device_info_t **devices = (* deviceList);
	uint32_t idSupportedDevice, numSupportedDevices;

	numSupportedDevices = devices[0]->numSupportedDevices;
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		if(devices[idSupportedDevice] != NULL){
			free(devices[idSupportedDevice]);
			devices[idSupportedDevice] = NULL;
		}
	}

    MYERS_ERROR(freeDevicesListHost(&devices));

    (* deviceList) = devices;
    return(SUCCESS);
}

MYERS_INLINE void endMyers(void ***myersBuffer)
{
	buffer_t **mBuff = (buffer_t **) (* myersBuffer);
	uint32_t numSupportedDevices, idSupportedDevice, numBuffers, idBuffer;
	device_info_t **devices = NULL;

	numBuffers = mBuff[0]->numBuffers;

	/*recollect all the supported devices */
	numSupportedDevices = mBuff[0]->device->numSupportedDevices;
	devices = (device_info_t **) malloc(numSupportedDevices * sizeof(device_info_t *));
		if (devices == NULL) MYERS_ERROR(E_ALLOCATE_MEM);

	devices[0] = mBuff[0]->device;
	idSupportedDevice = 0;
	for(idBuffer = 1; idBuffer < numBuffers; ++idBuffer){
		 if(devices[idSupportedDevice] != mBuff[idBuffer]->device){
			 idSupportedDevice++;
			 devices[idSupportedDevice] = mBuff[idBuffer]->device;
		 }
	}

	/* Synchronize all the Devices to the Host */
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
	    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
		CUDA_ERROR(cudaDeviceSynchronize());
	}

	/* Free all the references */
	MYERS_ERROR(freeReference(&mBuff[0]->reference, devices));

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
	    CUDA_ERROR(cudaSetDevice(mBuff[idBuffer]->device->idDevice));
		MYERS_ERROR(freeQueries(&mBuff[idBuffer]->queries));
		MYERS_ERROR(freeCandidates(&mBuff[idBuffer]->candidates));
		MYERS_ERROR(freeReorderBuffer(&mBuff[idBuffer]->reorderBuffer));
		MYERS_ERROR(freeResults(&mBuff[idBuffer]->results));

		CUDA_ERROR(cudaStreamDestroy(mBuff[idBuffer]->idStream));
	}

	/* reset all the device environments */
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
	    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
		CUDA_ERROR(cudaDeviceReset());
	}

	MYERS_ERROR(freeDevicesInfo(&devices));

	if(mBuff != NULL){
    	free(mBuff);
    	mBuff = NULL;
    }

	(* myersBuffer) = (void **) mBuff;
}
