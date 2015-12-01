
#include "../include/gpu_buffers.h"

/************************************************************
Functions to get information from the system
************************************************************/

GPU_INLINE uint32_t gpu_buffer_get_id_device_(const void* const gpuBuffer)
{
	const gpu_buffer_t* const mBuff = (gpu_buffer_t *) gpuBuffer;
	const uint32_t idSupDevice = mBuff->idSupportedDevice;
	return(mBuff->device[idSupDevice]->idDevice);
}

GPU_INLINE uint32_t gpu_buffer_get_id_supported_device_(const void* const gpuBuffer)
{
	const gpu_buffer_t* const mBuff = (gpu_buffer_t *) gpuBuffer;
	return(mBuff->idSupportedDevice);
}

GPU_INLINE gpu_error_t gpu_free_buffer(gpu_buffer_t *mBuff)
{
    if(mBuff->h_rawData != NULL){
        CUDA_ERROR(cudaFreeHost(mBuff->h_rawData));
        mBuff->h_rawData = NULL;
    }

    if(mBuff->d_rawData != NULL){
        CUDA_ERROR(cudaFree(mBuff->d_rawData));
        mBuff->d_rawData = NULL;
    }

    return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_get_min_memory_size_per_buffer(size_t *bytesPerBuffer)
{
	const uint32_t averarageNumPEQEntries = 1;
	const uint32_t candidatesPerQuery     = 1;

	const size_t bytesPerBPMBuffer    = GPU_BPM_MIN_ELEMENTS        * gpu_bpm_size_per_candidate(averarageNumPEQEntries,candidatesPerQuery);
	const size_t bytesPerSearchBuffer = GPU_FMI_SEARCH_MIN_ELEMENTS * gpu_fmi_decode_input_size();
	const size_t bytesPerDecodeBuffer = GPU_FMI_DECODE_MIN_ELEMENTS * gpu_fmi_search_input_size();

	(* bytesPerBuffer) = GPU_MAX(bytesPerBPMBuffer,GPU_MAX(bytesPerSearchBuffer,bytesPerDecodeBuffer));
	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_schedule_buffers(gpu_buffer_t ***gpuBuffer, const uint32_t numBuffers, gpu_device_info_t **device,
										   gpu_reference_buffer_t *reference, gpu_index_buffer_t *index, float maxMbPerBuffer)
{
	uint32_t idSupportedDevice, numBuffersPerDevice, idLocalBuffer;
	const size_t maxBytesPerBuffer = GPU_CONVERT_MB_TO_B(maxMbPerBuffer);
	const uint32_t numSupportedDevices = device[0]->numSupportedDevices;
	int32_t remainderBuffers = numBuffers, idGlobalBuffer = 0;

	gpu_buffer_t **buffer = (gpu_buffer_t **) malloc(numBuffers * sizeof(gpu_buffer_t *));
	if (buffer == NULL) GPU_ERROR(E_ALLOCATE_MEM);

	/* Assigning buffers for each GPU (to adapt the workload) */
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		size_t bytesPerDevice, bytesPerBuffer, minimumMemorySize;
		const uint32_t idDevice = device[idSupportedDevice]->idDevice;
		const size_t freeDeviceMemory = gpu_get_device_free_memory(idDevice);

		numBuffersPerDevice = GPU_ROUND(numBuffers * device[idSupportedDevice]->relativePerformance);
		if(idSupportedDevice == numSupportedDevices-1) numBuffersPerDevice = remainderBuffers;

		/* Resize the buffers to the available memory (Best effort system) */
		if(maxBytesPerBuffer != 0) bytesPerDevice = GPU_MIN(numBuffersPerDevice * maxBytesPerBuffer, freeDeviceMemory);
			else bytesPerDevice = freeDeviceMemory;

		bytesPerBuffer = bytesPerDevice / numBuffersPerDevice;

    	GPU_ERROR(gpu_get_min_memory_per_module(&minimumMemorySize, reference, index, numBuffers, GPU_NONE_MODULES));
		if(freeDeviceMemory < minimumMemorySize) GPU_ERROR(E_INSUFFICIENT_MEM_GPU);

		for(idLocalBuffer = 0; idLocalBuffer < numBuffersPerDevice; ++idLocalBuffer){
			buffer[idGlobalBuffer] = (gpu_buffer_t *) malloc(sizeof(gpu_buffer_t));
		    CUDA_ERROR(cudaSetDevice(idDevice));
		    GPU_ERROR(gpu_configure_buffer(buffer[idGlobalBuffer], idGlobalBuffer, idSupportedDevice, bytesPerBuffer, numBuffers,
										   device, reference, index));
			idGlobalBuffer++;
		}
		remainderBuffers -= numBuffersPerDevice;
	  }

	(* gpuBuffer) = buffer;
	return (SUCCESS);
}

GPU_INLINE void gpu_init_buffers_(gpu_buffers_dto_t *buff, gpu_index_dto_t *rawIndex, gpu_reference_dto_t* rawRef, gpu_info_dto_t *sys, const bool verbose)
{
	/* Buffer info */
	const float 				maxMbPerBuffer  		= buff->maxMbPerBuffer;
	const uint32_t 				numBuffers 				= buff->numBuffers;
	const gpu_module_t 			activeModules 			= buff->activeModules;
	/* System info */
	const gpu_dev_arch_t 		selectedArchitectures 	= sys->selectedArchitectures;
	const gpu_data_location_t	userAllocOption 		= sys->userAllocOption;
	const uint32_t				numSupportedDevices 	= gpu_get_num_supported_devices_(selectedArchitectures);
	/* Internal buffers info */
	gpu_buffer_t				**buffer 				= NULL;
	gpu_reference_buffer_t		*reference 				= NULL;
	gpu_index_buffer_t			*index	 				= NULL;
	gpu_device_info_t			**devices 				= NULL;
	/* Index info */
	const void 	 		    	*fmiRaw 				= rawIndex->fmi;
	const gpu_index_coding_t	indexCoding				= rawIndex->indexCoding;
	const uint64_t 				bwtSize 				= rawIndex->bwtSize;
	/* Reference info */
	const char 		 			*referenceRaw 			= rawRef->reference;
	const gpu_ref_coding_t 		refCoding 				= rawRef->refCoding;
	const uint64_t 		 		refSize 				= rawRef->refSize;

	GPU_ERROR(gpu_fast_driver_awake());
	GPU_ERROR(gpu_init_reference(&reference, referenceRaw, refSize, refCoding, numSupportedDevices, activeModules));
	GPU_ERROR(gpu_init_index(&index, fmiRaw, bwtSize, indexCoding, numSupportedDevices, activeModules));

	GPU_ERROR(gpu_configure_modules(&devices, selectedArchitectures, userAllocOption, numBuffers, reference, index));
	GPU_ERROR(gpu_set_devices_local_memory(devices, cudaFuncCachePreferL1));

	GPU_ERROR(gpu_transfer_reference_CPU_to_GPUs(reference, devices));
	GPU_ERROR(gpu_transfer_index_CPU_to_GPUs(index, devices));

	GPU_ERROR(gpu_schedule_buffers(&buffer, numBuffers, devices, reference, index, maxMbPerBuffer));

	GPU_ERROR(gpu_free_unused_reference_host(reference, devices));
	GPU_ERROR(gpu_free_unused_index_host(index, devices));

	buff->buffer = (void **) buffer;
}

GPU_INLINE void gpu_destroy_buffers_(gpu_buffers_dto_t *buff)
{
	gpu_buffer_t **mBuff = (gpu_buffer_t **) buff->buffer;
	uint32_t idBuffer;
	gpu_device_info_t **devices = mBuff[0]->device;

	const uint32_t numBuffers = mBuff[0]->numBuffers;

	GPU_ERROR(gpu_synchronize_all_devices(devices));

	/* Free all the references */
	GPU_ERROR(gpu_free_reference(&mBuff[0]->reference, devices));
	GPU_ERROR(gpu_free_index(&mBuff[0]->index, devices));

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		const uint32_t idSupDevice = mBuff[idBuffer]->idSupportedDevice;
	    CUDA_ERROR(cudaSetDevice(devices[idSupDevice]->idDevice));
		GPU_ERROR(gpu_free_buffer(mBuff[idBuffer]));
		CUDA_ERROR(cudaStreamDestroy(mBuff[idBuffer]->idStream));
		free(mBuff[idBuffer]);
	}

	GPU_ERROR(gpu_reset_all_devices(devices));
	GPU_ERROR(gpu_free_devices_info(devices));

	if(mBuff != NULL){
    	free(mBuff);
    	mBuff = NULL;
    	buff->buffer = NULL;
    }
}

/************************************************************
Primitives to schedule and manage the buffers
************************************************************/

GPU_INLINE void gpu_alloc_buffer_(void *gpuBuffer)
{
	gpu_buffer_t *mBuff = (gpu_buffer_t *) gpuBuffer;
	const uint32_t idSupDevice = mBuff->idSupportedDevice;
	mBuff->h_rawData = NULL;
	mBuff->d_rawData = NULL;

	//Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

	//ALLOCATE HOST AND DEVICE BUFFER
	CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));
}

GPU_INLINE gpu_error_t gpu_configure_buffer(gpu_buffer_t *mBuff, const uint32_t idBuffer, const uint32_t idSupportedDevice, const size_t bytesPerBuffer,
											const uint32_t numBuffers, gpu_device_info_t **device, gpu_reference_buffer_t *reference, gpu_index_buffer_t *index)
{
	/* Buffer information */
	mBuff->idBuffer 		 = idBuffer;
	mBuff->numBuffers 		 = numBuffers;
	mBuff->idSupportedDevice = idSupportedDevice;
	mBuff->device			 = device;
	mBuff->sizeBuffer		 = bytesPerBuffer;
	mBuff->typeBuffer		 = GPU_NONE_MODULES;

	/* Module structures */
	mBuff->index			 = index;
	mBuff->reference		 = reference;

	/* Chunk of RAW memory for the buffer */
	mBuff->h_rawData		 = NULL;
	mBuff->d_rawData 		 = NULL;

	/* Set in which Device we create and initialize the structures */
    CUDA_ERROR(cudaSetDevice(mBuff->device[idSupportedDevice]->idDevice));

	/* Create the CUDA stream per each buffer */
	CUDA_ERROR(cudaStreamCreate(&mBuff->idStream));

	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_get_min_memory_per_module(size_t *minimumMemorySize, gpu_reference_buffer_t *reference, gpu_index_buffer_t* index,
													 const uint32_t numBuffers, const gpu_module_t activeModules)
{
	size_t memorySize;
	size_t bytesPerReference, bytesPerIndex, bytesPerBuffer, bytesPerAllBuffers;

	gpu_get_min_memory_size_per_buffer(&bytesPerBuffer);

	bytesPerAllBuffers  = numBuffers * bytesPerBuffer;
	bytesPerReference 	= reference->numEntries * GPU_REFERENCE_BYTES_PER_ENTRY;
	bytesPerIndex     	= index->numEntries * sizeof(gpu_fmi_entry_t);

	memorySize = bytesPerAllBuffers;
	if(activeModules &  GPU_BPM) memorySize += bytesPerReference;
	if(activeModules & (GPU_FMI_EXACT_SEARCH |  GPU_FMI_DECODE_POS)) memorySize += bytesPerIndex;

	(* minimumMemorySize) = memorySize;
	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_configure_modules(gpu_device_info_t ***devices, const gpu_dev_arch_t selectedArchitectures,
											 const gpu_data_location_t userAllocOption, const uint32_t numBuffers,
											 gpu_reference_buffer_t *reference, gpu_index_buffer_t *index)
{
	const uint32_t numDevices = gpu_get_num_devices();
	const uint32_t numSupportedDevices = gpu_get_num_supported_devices_(selectedArchitectures);
	uint32_t idDevice, idSupportedDevice = 0;

	gpu_device_info_t **dev = (gpu_device_info_t **) malloc(numSupportedDevices * sizeof(gpu_device_info_t *));

	for(idDevice = 0, idSupportedDevice = 0; idDevice < numDevices; ++idDevice){
		size_t recomendedMemorySize, requiredMemorySize;
		bool localReference = true, localIndex = true, dataFitsMemoryDevice = true;
		const bool deviceArchSupported = gpu_get_device_architecture(idDevice) & selectedArchitectures;
		if(deviceArchSupported){
			GPU_ERROR(gpu_module_memory_manager(idDevice, idSupportedDevice, numBuffers,
												userAllocOption, &dataFitsMemoryDevice, &localReference, &localIndex,
												&recomendedMemorySize, &requiredMemorySize, reference, index));
			if(dataFitsMemoryDevice){
				GPU_ERROR(gpu_init_device(&dev[idDevice], idDevice, idSupportedDevice, selectedArchitectures));
				idSupportedDevice++;
			}
		}
		GPU_ERROR(gpu_screen_status_device(idDevice, deviceArchSupported,
										   dataFitsMemoryDevice, localReference, localIndex,
										   recomendedMemorySize, requiredMemorySize));
	}

	GPU_ERROR(gpu_characterize_devices(dev, idSupportedDevice));

	(* devices) = dev;
	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_module_memory_manager(const uint32_t idDevice, const uint32_t idSupDevice, const uint32_t numBuffers,
												 const gpu_data_location_t userAllocOption, bool *dataFits,
												 bool *lReference, bool *lIndex, size_t *recMemSize, size_t *reqMemSize,
												 gpu_reference_buffer_t *reference, gpu_index_buffer_t *index)
{
	const gpu_module_t activeModules = reference->activeModules | index->activeModules;
	bool localReference = true, localIndex = true, dataFitsMemoryDevice = true;
	size_t memoryFree = gpu_get_device_free_memory(idDevice);
    memory_alloc_t referenceModule = GPU_NONE_MAPPED, indexModule = GPU_NONE_MAPPED;
	size_t minimumMemorySize;

    if (userAllocOption == GPU_REMOTE_DATA){
    	const gpu_module_t localActiveModules = GPU_NONE_MODULES;
    	localReference 		 = false;
    	localIndex 	   		 = false;
    	GPU_ERROR(gpu_get_min_memory_per_module(&minimumMemorySize, reference, index, numBuffers, localActiveModules));
		if (memoryFree < minimumMemorySize) dataFitsMemoryDevice = false;
    }else if (userAllocOption == GPU_LOCAL_DATA){
    	const gpu_module_t localActiveModules = GPU_ALL_MODULES & activeModules;
    	GPU_ERROR(gpu_get_min_memory_per_module(&minimumMemorySize, reference, index, numBuffers, localActiveModules));
		if (memoryFree < minimumMemorySize) dataFitsMemoryDevice = false;
	}else{ // GPU_LOCAL_OR_REMOTE_DATA by default
		if(dataFitsMemoryDevice){ 		/* Put all the structures in the device (LOCAL) */
	    	const gpu_module_t localActiveModules = GPU_ALL_MODULES & activeModules;
	    	GPU_ERROR(gpu_get_min_memory_per_module(&minimumMemorySize, reference, index, numBuffers, localActiveModules));
	    	if(memoryFree < minimumMemorySize) localReference = false;
		}
		if(!localReference){ 			/* Put the reference in the host (REMOTE) */
	    	const gpu_module_t localActiveModules = (GPU_FMI_DECODE_POS | GPU_FMI_EXACT_SEARCH) & activeModules;
	    	GPU_ERROR(gpu_get_min_memory_per_module(&minimumMemorySize, reference, index, numBuffers, localActiveModules));
	    	if(memoryFree < minimumMemorySize) localIndex = false;
		}
		if(!localIndex){				/* Put the reference and the index in the host (REMOTE) */
	    	const gpu_module_t localActiveModules = GPU_NONE_MODULES;
	    	GPU_ERROR(gpu_get_min_memory_per_module(&minimumMemorySize, reference, index, numBuffers, localActiveModules));
	    	if(memoryFree < minimumMemorySize) dataFitsMemoryDevice = false;
		}
	}

    /* Configure the index & reference mapping side (HOST / DEVICE / NONE) */
    if(dataFitsMemoryDevice && (activeModules & GPU_BPM)){
    	if(localReference) referenceModule = GPU_DEVICE_MAPPED;
    	else referenceModule = GPU_HOST_MAPPED;
    }
    if(dataFitsMemoryDevice && (activeModules & (GPU_FMI_EXACT_SEARCH | GPU_FMI_DECODE_POS))){
       	if(localIndex) indexModule = GPU_DEVICE_MAPPED;
       	else indexModule = GPU_HOST_MAPPED;
    }

    if(activeModules & GPU_BPM)
    	reference->memorySpace[idSupDevice] = referenceModule;
    if(activeModules & (GPU_FMI_EXACT_SEARCH | GPU_FMI_DECODE_POS))
    	index->memorySpace[idSupDevice] 	= indexModule;

	/* Return the memory sizes and flags */
	GPU_ERROR(gpu_get_min_memory_per_module(reqMemSize, reference, index, numBuffers, GPU_NONE_MODULES));
	GPU_ERROR(gpu_get_min_memory_per_module(recMemSize, reference, index, numBuffers, activeModules));

	(* lReference)  = localReference;
	(* lIndex)		= localIndex;
	(* dataFits) 	= dataFitsMemoryDevice;

	return(SUCCESS);
}
