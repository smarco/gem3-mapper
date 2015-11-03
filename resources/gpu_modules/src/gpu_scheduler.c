/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Host scheduler for BPM on GPU
 */

#include "../include/gpu_scheduler.h"

/************************************************************
Functions to init all the Myers resources
************************************************************/

GPU_INLINE uint32_t gpu_buffer_get_id_device_(void *gpuBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) gpuBuffer;
	return(mBuff->device->idDevice);
}

GPU_INLINE uint32_t gpu_buffer_get_id_supported_device_(void *gpuBuffer){
	gpu_buffer_t *mBuff = (gpu_buffer_t *) gpuBuffer;
	return(mBuff->device->idSupportedDevice);
}

GPU_INLINE size_t gpu_get_device_free_memory(uint32_t idDevice)
{
	size_t free, total;
    CUDA_ERROR(cudaSetDevice(idDevice));
	CUDA_ERROR(cudaMemGetInfo(&free, &total));
	return (free);
}

GPU_INLINE gpu_dev_arch_t gpu_get_device_architecture(uint32_t idDevice)
{
	struct cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, idDevice);

	if (devProp.major <= 1) return(GPU_ARCH_TESLA);								/* CC 1.X		*/
	if (devProp.major == 2 && devProp.minor == 0) return(GPU_ARCH_FERMI_1G); 	/* CC 2.0		*/
	if (devProp.major == 2 && devProp.minor >  0) return(GPU_ARCH_FERMI_2G); 	/* CC 2.1		*/
	if (devProp.major == 3 && devProp.minor <  5) return(GPU_ARCH_KEPLER_1G); 	/* CC 3.0, 3.2	*/
	if (devProp.major == 3 && devProp.minor >= 5) return(GPU_ARCH_KEPLER_2G); 	/* CC 3.5		*/
	if (devProp.major == 5 && devProp.minor == 0) return(GPU_ARCH_MAXWELL_1G);  /* CC 5.0		*/
	if (devProp.major == 5 && devProp.minor >= 5) return(GPU_ARCH_MAXWELL_2G);  /* CC 5.2		*/
	if (devProp.major == 6 && devProp.minor == 0) return(GPU_ARCH_PASCAL_1G);   /* CC 6.0		*/
	if (devProp.major == 6 && devProp.minor >= 5) return(GPU_ARCH_PASCAL_2G);   /* CC 6.2		*/

	return(GPU_ARCH_NEWGEN);													/* CC X.X		*/
}

GPU_INLINE uint32_t gpu_get_SM_cuda_cores(gpu_dev_arch_t architecture)
{
	switch (architecture) {
		case GPU_ARCH_TESLA:		return(8);
		case GPU_ARCH_FERMI_1G:		return(32);
		case GPU_ARCH_FERMI_2G:		return(48);
		case GPU_ARCH_KEPLER_1G:	return(192);
		case GPU_ARCH_KEPLER_2G:	return(192);
		case GPU_ARCH_MAXWELL_1G:	return(128);
		case GPU_ARCH_MAXWELL_2G:	return(128);
		case GPU_ARCH_PASCAL_1G:	return(128);
		case GPU_ARCH_PASCAL_2G:	return(128);
		default:					return(128);
	}
}

GPU_INLINE uint32_t get_device_cuda_cores(uint32_t idDevice)
{
	uint32_t coresPerSM;
	gpu_dev_arch_t architecture;

	struct cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, idDevice);

	architecture = gpu_get_device_architecture(idDevice);
	coresPerSM = gpu_get_SM_cuda_cores(architecture);

	return(devProp.multiProcessorCount * coresPerSM);
}

GPU_INLINE uint32_t gpu_get_num_devices(){
	uint32_t numDevices;
	CUDA_ERROR(cudaGetDeviceCount(&numDevices));
	return(numDevices);
}

GPU_INLINE gpu_error_t gpu_get_num_supported_devices_(gpu_dev_arch_t selectedArchitectures)
{
	uint32_t idDevice, numDevices, numSupportedDevices = 0;
	gpu_dev_arch_t deviceArch;

	CUDA_ERROR(cudaGetDeviceCount(&numDevices));

	for(idDevice = 0; idDevice < numDevices; ++idDevice){
		deviceArch = gpu_get_device_architecture(idDevice);
	    if(deviceArch & selectedArchitectures) numSupportedDevices++;
	}

	return(numSupportedDevices);
}

GPU_INLINE gpu_error_t gpu_get_min_memory_size_per_buffer(size_t *bytesPerBuffer)
{
	const uint32_t averarageNumPEQEntries = 1;
	const uint32_t candidatesPerQuery     = 1;

	bytesPerBPMBuffer    = GPU_BPM_MIN_ELEMENTS        * gpu_bpm_size_per_candidate(averarageNumPEQEntries,candidatesPerQuery);
	bytesPerSearchBuffer = GPU_FMI_SEARCH_MIN_ELEMENTS * gpu_fmi_decode_input_size();
	bytesPerDecodeBuffer = GPU_FMI_DECODE_MIN_ELEMENTS * gpu_fmi_search_input_size();

	(* bytesPerBuffer) = GPU_MAX(bytesPerBPMBuffer,GPU_MAX(bytesPerBPMBuffer,bytesPerDecodeBuffer));
	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_get_min_memory_size_per_device(size_t *minimumMemorySize, gpu_reference_buffer_t *reference,
													  	  gpu_index_buffer_t* index,  const uint32_t numBuffers,
													  	  const gpu_module_t activeModules)
{
	size_t memorySize;
	size_t bytesPerReference, bytesPerIndex, bytesPerBuffer;

	gpu_get_min_memory_size_per_buffer(&bytesPerBuffer);

	bytesPerBuffers   = numBuffers * bytesPerBuffer;
	bytesPerReference = reference->numEntries * GPU_REFERENCE_BYTES_PER_ENTRY;
	bytesPerIndex     = index->numEntries * sizeof(gpu_fmi_entry_t);

	memorySize = bytesPerBuffers;
	if(activeModules &  GPU_BPM) memorySize += bytesPerReference;
	if(activeModules & (GPU_FMI_EXACT_SEARCH |  GPU_FMI_DECODE_POS)) memorySize += bytesPerIndex;

	(* minimumMemorySize) = memorySize;
	return (SUCCESS);
}


/************************************************************
Functions to init all the Myers resources
************************************************************/

GPU_INLINE gpu_error_t gpu_set_device_local_memory(gpu_device_info_t **devices, enum cudaFuncCache cacheConfig)
{
	uint32_t idSupportedDevice, numSupportedDevices = devices[0]->numSupportedDevices;

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
	    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
		CUDA_ERROR(cudaDeviceSetCacheConfig(cacheConfig));
	}

	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_fast_driver_awake()
{
	//Dummy call to the NVIDIA API to awake earlier the driver.
	void* prt = NULL;
	cudaMalloc(&prt, 0);

	return(SUCCESS);
}


/************************************************************
Functions to init all the Myers resources
************************************************************/

GPU_INLINE gpu_error_t gpu_alloc_buffer_(gpu_buffer_t *mBuff)
{
	mBuff->h_rawData = NULL;
	mBuff->d_rawData = NULL;

	//Select the device of the Multi-GPU platform
    CUDA_ERROR(cudaSetDevice(mBuff->device->idDevice));

	//Create the CUDA stream per each buffer
	CUDA_ERROR(cudaStreamCreate(&buffer->idStream));

	//ALLOCATE HOST AND DEVICE BUFFER
	CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
	CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));

	return(SUCCESS);
}



/************************************************************
Functions to init all the Myers resources
************************************************************/

//GPU_INLINE gpu_error_t gpu_init_device(gpu_device_info_t *device, uint32_t idDevice, uint32_t idSupportedDevice)
//{
//	if (localReference) reference->memorySpace[idSupportedDevice] = GPU_DEVICE_MAPPED;
//		else reference->memorySpace[idSupportedDevice] = GPU_HOST_MAPPED;
//
//	if (localIndex) index->memorySpace[idSupportedDevice] = GPU_DEVICE_MAPPED;
//		else index->memorySpace[idSupportedDevice] = GPU_HOST_MAPPED;
//
//    if ( userReferenceAllocOption == GPU_REMOTE_DATA){
//    	localReference = false;
//    }
//
//	if (userReferenceAllocOption  == GPU_LOCAL_OR_REMOTE_DATA){
//		if (memoryFree < minimumMemorySize) localReference = false;
//	}
//
//	if (userReferenceAllocOption == GPU_LOCAL_DATA){
//		if(memoryFree < minimumMemorySize) dataFitsMemoryDevice = false;
//	}
//
//	return(SUCCESS);
//}

GPU_INLINE gpu_error_t gpu_init_device(gpu_device_info_t **device, uint32_t idDevice, uint32_t idSupportedDevice)
{
	gpu_device_info_t *dev = NULL;
	struct cudaDeviceProp devProp;

	CUDA_ERROR(cudaGetDeviceProperties(&devProp, idDevice));

	dev = (device_info_t *) malloc(sizeof(gpu_device_info_t));
	if(dev == NULL) return(E_ALLOCATE_MEM);

	dev->numDevices 		 = gpu_get_num_devices();
	dev->numSupportedDevices = gpu_get_num_supported_devices_(selectedArchitectures);
	dev->idSupportedDevice 	 = idSupportedDevice;
	dev->idDevice 			 = idDevice;
	dev->architecture 		 = gpu_get_device_architecture(idDevice);
	dev->cudaCores 			 = gpu_get_device_cuda_cores(idDevice);
	dev->coreClockRate		 = devProp.clockRate;
	dev->memoryBusWidth		 = devProp.memoryBusWidth;
	dev->memoryClockRate	 = devProp.memoryClockRate;

	dev->absolutePerformance = dev->cudaCores * dev->coreClockRate;
	dev->absoluteBandwidth   = 2.0 * dev->memoryClockRate * (dev->memoryBusWidth / 8) / 1000;

	(* device) = dev;
	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_screen_status_device(const char * deviceName, const uint32_t idDevice,
												const uint32_t memoryFree, const uint32_t recomendedMemorySize,
												const bool deviceArchSupported, const bool dataFitsMemoryDevice,
												const bool localIndex, const bool localReference)
{
	fprintf(stderr, "GPU Device %d: %s", idDevice, deviceName);
	if((deviceArchSupported && dataFitsMemoryDevice)){
		if ((localReference && localIndex)  == false){
			fprintf(stderr, "\t \t Selected Device [RUNNING] \n");
			fprintf(stderr, "\t \t WARNING PERF. SLOWDOWNS (Mem Recommended: %u MBytes - Mem Available: %u MBytes) \n",
					recomendedMemorySize, memoryFree);
		}else{
			fprintf(stderr, "\t \t Selected Device [RUNNING] \n");
		}
	}
	if(!deviceArchSupported)
		fprintf(stderr, "\t \t UNSUPPORTED Architecture [Compute Capability < 2.0: UNSUITABLE] \n");
	if(localReference && (memoryFree < minimumMemorySize))
		fprintf(stderr, "\t \t INSUFFICIENT DEVICE MEMORY (Mem Required: %u MBytes - Mem Available: %u MBytes) \n",
				minimumMemorySize, memoryFree);

	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_characterize_devices(gpu_device_info_t **dev)
{
	const unint32_t numSupportedDevices = dev->numSupportedDevices[0];
	uint32_t idSupportedDevice, allSystemPerformance = 0, allSystemBandwidth = 0;

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		allSystemPerformance += dev[idSupportedDevice]->absolutePerformance;
		allSystemBandwidth += dev[idSupportedDevice]->absoluteBandwidth;
	}

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		dev[idSupportedDevice]->relativeBandwidth  = (float) dev[idSupportedDevice]->absoluteBandwidth / (float) allSystemBandwidth;
		dev[idSupportedDevice]->allSystemBandwidth = allSystemBandwidth;

		dev[idSupportedDevice]->relativePerformance  = (float) dev[idSupportedDevice]->absolutePerformance / (float) allSystemPerformance;
		dev[idSupportedDevice]->allSystemPerformance = allSystemPerformance;
	}
	return(SUCCESS);
}


GPU_INLINE gpu_error_t gpu_configure_supported_devices(gpu_device_info_t ***devices, const gpu_dev_arch_t selectedArchitectures,
													  const gpu_data_location_t userReferenceAllocOption, gpu_index_buffer_t *index,
													  gpu_reference_buffer_t *reference)
{
	const uint32_t numDevices = gpu_get_num_devices();
	const uint32_t numSupportedDevices = gpu_get_num_supported_devices_(selectedArchitectures);
	uint32_t idDevice, idSupportedDevice;

	gpu_device_info_t **dev = (device_info_t **) malloc(numSupportedDevices * sizeof(gpu_device_info_t *));
	if (dev == NULL) return(E_ALLOCATE_MEM);

	// Sanity check
	dev[0]->numSupportedDevices = 0;

	for(idDevice = 0, idSupportedDevice = 0; idDevice < numDevices; ++idDevice){
		bool localReference = true, localIndex = true, dataFitsMemoryDevice = true;
		const bool deviceArchSupported = get_device_architecture(idDevice) & GPU_ARCH_SUPPORTED;
		if(deviceArchSupported){
			size_t memoryFree = gpu_get_device_free_memory(idDevice);

			if(dataFitsMemoryDevice){
				GPU_ERROR(gpu_init_device(dev, idDevice, idSupportedDevice));
				idSupportedDevice++;
			}
		}
		GPU_ERROR(gpu_screen_status_device(memoryFree, recomendedMemorySize, deviceArchSupported, dataFitsMemoryDevice, localIndex, localReference));
	}

	GPU_ERROR(gpu_characterize_devices(dev));

	(* devices) = dev;
	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_configure_device_buffers(gpu_buffer_t ***gpuBuffer, const uint32_t numBuffers,
													gpu_reference_buffer_t *reference, gpu_index_buffer_t *index,
												    gpu_device_info_t **device, uint32_t maxMbPerBuffer, const gpu_module_t activeModules,
												    const bool verbose)
{
	uint32_t numSupportedDevices, idSupportedDevice, idDevice, idGlobalBuffer, numBuffersPerDevice, idLocalBuffer;
	uint32_t freeDeviceMemory, maxMbPerDevice;
	size_t bytesPerBuffer;
	float bytesPerCandidate;
	int32_t remainderBuffers;

	gpu_buffer_t **buffer = (gpu_buffer_t **) malloc(numBuffers * sizeof(gpu_buffer_t *));
	if (buffer == NULL) GPU_ERROR(E_ALLOCATE_MEM);

	numSupportedDevices = device[0]->numSupportedDevices;
	remainderBuffers = numBuffers;
	idGlobalBuffer = 0;

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		idDevice = device[idSupportedDevice]->idDevice;
		freeDeviceMemory = gpu_get_device_free_memory(idDevice);

		numBuffersPerDevice = GPU_ROUND(numBuffers * device[idSupportedDevice]->relativePerformance);
		if(idSupportedDevice == numSupportedDevices-1) numBuffersPerDevice = remainderBuffers;

		if(maxMbPerBuffer != 0) maxMbPerDevice = GPU_MIN(numBuffersPerDevice * maxMbPerBuffer, freeDeviceMemory);
			else maxMbPerDevice = freeDeviceMemory;

		bytesPerBuffer = (GPU_CONVERT_MB_TO_B(maxMbPerDevice) * 0.95) / numBuffersPerDevice;

		maxCandidates = (uint32_t) ((bytesPerBuffer - bytesPadding) / bytesPerCandidate);
		maxQueries    = maxCandidates / candidatesPerQuery;
		maxPEQEntries = maxQueries * averageNumPEQEntries;

		//if((maxCandidates < 1) || (maxQueries < 1) || (GPU_CONVERT_B_TO_MB(bytesPerBuffer * numBuffersPerDevice) > freeDeviceMemory))
		//	GPU_ERROR(E_INSUFFICIENT_MEM_GPU);

		for(idLocalBuffer = 0; idLocalBuffer < numBuffersPerDevice; ++idLocalBuffer){
			buffer[idGlobalBuffer] = (buffer_t *) malloc(sizeof(buffer_t));
		    CUDA_ERROR(cudaSetDevice(idDevice));
			GPU_ERROR(gpu_configure_buffer(buffer[idGlobalBuffer], idGlobalBuffer, device[idSupportedDevice], numBuffers, maxCandidates,
						  maxQueries, maxPEQEntries, bucketPaddingCandidates, reference));
			idGlobalBuffer++;
		}
		remainderBuffers -= numBuffersPerDevice;
	  }

	(* gpuBuffer) = buffer;
	return (SUCCESS);
}

GPU_INLINE void gpu_init_(void ***gpuBuffer, uint32_t numBuffers, uint32_t maxMbPerBuffer,
						  const char *referenceRaw, gpu_ref_coding_t refCoding, const uint64_t refSize,
						  const char *indexRaw, gpu_index_coding_t indexCoding, const uint64_t bwtSize,
						  gpu_dev_arch_t selectedArchitectures, gpu_data_location_t userReferenceAllocOption,
						  const bool verbose)
{
	gpu_buffer_t			**buffer 	= NULL;
	gpu_reference_buffer_t	*reference 	= NULL;
	gpu_index_buffer_t		*index	 	= NULL;
	gpu_device_info_t		**devices 	= NULL;
	size_t					minimumMemorySize = 0;

	GPU_ERROR(gpu_fast_driver_awake());
	GPU_ERROR(gpu_init_reference(&reference, referenceRaw, refSize, refCoding));
	GPU_ERROR(gpu_init_index(&index, indexRaw, bwtSize, indexCoding));

	GPU_ERROR(gpu_configure_supported_devices(&devices, selectedArchitectures, userReferenceAllocOption, verbose)); //TODO review
	GPU_ERROR(gpu_set_device_local_memory(devices, cudaFuncCachePreferL1));

	GPU_ERROR(gpu_transfer_reference_CPU_to_GPUs(reference, devices));
	GPU_ERROR(gpu_transfer_index_CPU_to_GPUs(index, devices));

	GPU_ERROR(gpu_configure_device_buffers(&buffer, numBuffers, reference, devices, maxMbPerBuffer, verbose)); //TODO review

	GPU_ERROR(gpu_free_unused_reference_host(reference, devices));
	GPU_ERROR(gpu_free_unused_index_host(index, devices));
	GPU_ERROR(gpu_free_devices_list_host(&devices));

	(* gpuBuffer) = (void **) buffer;
}

/************************************************************
Functions to free all the buffer resources (HOST & DEVICE)
************************************************************/

GPU_INLINE void gpu_reset_all_devices(gpu_buffer_t **mBuff)
{
	uint32_t numSupportedDevices = mBuff[0]->device->numSupportedDevices;

	/* reset all the device environments */
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
		CUDA_ERROR(cudaDeviceReset());
	}
}

GPU_INLINE void gpu_synchronize_all_devices(gpu_buffer_t **mBuff)
{
	const uint32_t numSupportedDevices = mBuff[0]->device->numSupportedDevices;
	uint32_t 	   idSupportedDevice;

	/* Synchronize all the Devices to the Host */
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
		CUDA_ERROR(cudaDeviceSynchronize());
	}
}

GPU_INLINE void gpu_build_supported_devices(gpu_buffer_t **mBuff, gpu_device_info_t ***dev)
{
	const uint32_t 	  numBuffers 		  = mBuff[0]->numBuffers;
	const uint32_t 	  numSupportedDevices = mBuff[0]->device->numSupportedDevices;
	gpu_device_info_t devices;
	uint32_t		  idBuffer, idSupportedDevice;

	/*recollect all the supported devices */
	devices = (device_info_t **) malloc(numSupportedDevices * sizeof(gpu_device_info_t *));
		if (devices == NULL) GPU_ERROR(E_ALLOCATE_MEM);

	devices[0] = mBuff[0]->device;
	idSupportedDevice = 0;
	for(idBuffer = 1; idBuffer < numBuffers; ++idBuffer){
		 if(devices[idSupportedDevice] != mBuff[idBuffer]->device){
			 idSupportedDevice++;
			 devices[idSupportedDevice] = mBuff[idBuffer]->device;
		 }
	}

	(* dev) = devices;
}

GPU_INLINE gpu_error_t gpu_free_buffer(gpu_buffer_t *mBuff)
{
    if(mBuff->h_rawData != NULL){
        CUDA_ERROR(cudaFreeHost(mBuff->h_rawData));
        mBuff->h_rawData = NULL;
    }

    if(queries->d_qinfo != NULL){
        CUDA_ERROR(cudaFree(mBuff->d_rawData));
        mBuff->h_rawData = NULL;
    }

    return(SUCCESS);
}

//ok
GPU_INLINE gpu_error_t gpu_free_devices_list_host(gpu_device_info_t ***deviceList)
{
	gpu_device_info_t **device = (* deviceList);

    if(device != NULL){
        free(device);
        device = NULL;
    }

    (* deviceList) = device;
    return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_free_devices_info(gpu_device_info_t ***deviceList)
{
	gpu_device_info_t **devices = (* deviceList);
	uint32_t idSupportedDevice, numSupportedDevices;

	numSupportedDevices = devices[0]->numSupportedDevices;
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		if(devices[idSupportedDevice] != NULL){
			free(devices[idSupportedDevice]);
			devices[idSupportedDevice] = NULL;
		}
	}

    GPU_ERROR(gpu_free_devices_list_host(&devices));

    (* deviceList) = devices;
    return(SUCCESS);
}

GPU_INLINE void gpu_destroy_(void ***gpuBuffer)
{
	gpu_buffer_t **mBuff = (gpu_buffer_t **) (* gpuBuffer);
	uint32_t idBuffer;
	gpu_device_info_t **devices = NULL;

	const uint32_t numBuffers = mBuff[0]->numBuffers;

	GPU_ERROR(gpu_build_supported_devices(mBuff, &devices));
	GPU_ERROR(gpu_synchronize_all_devices(mBuff));

	/* Free all the references */
	GPU_ERROR(gpu_free_reference(&mBuff[0]->reference, devices));
	GPU_ERROR(gpu_free_index(&mBuff[0]->index, devices));

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
	    CUDA_ERROR(cudaSetDevice(mBuff[idBuffer]->device->idDevice));
		GPU_ERROR(gpu_free_buffer(mBuff[idBuffer]));
		CUDA_ERROR(cudaStreamDestroy(mBuff[idBuffer]->idStream));
	}

	GPU_ERROR(gpu_reset_all_devices(mBuff));
	GPU_ERROR(gpu_free_devices_info(&devices));

	if(mBuff != NULL){
    	free(mBuff);
    	mBuff = NULL;
    }

	(* gpuBuffer) = (void **) mBuff;
}
