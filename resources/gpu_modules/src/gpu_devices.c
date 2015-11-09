
#include "../include/gpu_devices.h"

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

GPU_INLINE uint32_t gpu_get_device_cuda_cores(uint32_t idDevice)
{
	uint32_t coresPerSM;
	gpu_dev_arch_t architecture;

	struct cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, idDevice);

	architecture = gpu_get_device_architecture(idDevice);
	coresPerSM = gpu_get_SM_cuda_cores(architecture);

	return(devProp.multiProcessorCount * coresPerSM);
}

GPU_INLINE size_t gpu_get_device_free_memory(uint32_t idDevice)
{
	size_t free, total;
    CUDA_ERROR(cudaSetDevice(idDevice));
	CUDA_ERROR(cudaMemGetInfo(&free, &total));
	return (free * 0.95);
}

GPU_INLINE uint32_t gpu_get_num_devices()
{
	int32_t numDevices;
	CUDA_ERROR(cudaGetDeviceCount(&numDevices));
	return(numDevices);
}

GPU_INLINE uint32_t gpu_get_num_supported_devices_(gpu_dev_arch_t selectedArchitectures)
{
	uint32_t idDevice, numSupportedDevices = 0;
	int32_t numDevices;
	gpu_dev_arch_t deviceArch;


	CUDA_ERROR(cudaGetDeviceCount(&numDevices));

	for(idDevice = 0; idDevice < numDevices; ++idDevice){
		deviceArch = gpu_get_device_architecture(idDevice);
	    if(deviceArch & selectedArchitectures) numSupportedDevices++;
	}

	return(numSupportedDevices);
}

GPU_INLINE gpu_error_t gpu_characterize_devices(gpu_device_info_t **dev, uint32_t numSupportedDevices)
{
	uint32_t idSupportedDevice, allSystemPerformance = 0, allSystemBandwidth = 0;

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		allSystemPerformance += dev[idSupportedDevice]->absolutePerformance;
		allSystemBandwidth += dev[idSupportedDevice]->absoluteBandwidth;
	}

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		dev[idSupportedDevice]->numSupportedDevices = numSupportedDevices;
		dev[idSupportedDevice]->relativeBandwidth  = (float) dev[idSupportedDevice]->absoluteBandwidth / (float) allSystemBandwidth;
		dev[idSupportedDevice]->allSystemBandwidth = allSystemBandwidth;
		dev[idSupportedDevice]->relativePerformance  = (float) dev[idSupportedDevice]->absolutePerformance / (float) allSystemPerformance;
		dev[idSupportedDevice]->allSystemPerformance = allSystemPerformance;
	}
	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_reset_all_devices(gpu_device_info_t **devices)
{
	const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
	uint32_t idSupportedDevice;

	/* reset all the device environments */
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
		CUDA_ERROR(cudaDeviceReset());
	}
	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_synchronize_all_devices(gpu_device_info_t **devices)
{
	const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
	uint32_t 	   idSupportedDevice;

	/* Synchronize all the Devices to the Host */
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
		CUDA_ERROR(cudaDeviceSynchronize());
	}
	return(SUCCESS);
}

/************************************************************
Primitives to manage device driver options
************************************************************/

GPU_INLINE gpu_error_t gpu_set_devices_local_memory(gpu_device_info_t **devices, enum cudaFuncCache cacheConfig)
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
Primitives to schedule and manage the devices
************************************************************/

GPU_INLINE gpu_error_t gpu_init_device(gpu_device_info_t **device, uint32_t idDevice, uint32_t idSupportedDevice,
									   const gpu_dev_arch_t selectedArchitectures)
{
	gpu_device_info_t *dev = NULL;
	struct cudaDeviceProp devProp;

	CUDA_ERROR(cudaGetDeviceProperties(&devProp, idDevice));

	dev = (gpu_device_info_t *) malloc(sizeof(gpu_device_info_t));
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


GPU_INLINE gpu_error_t gpu_screen_status_device(const uint32_t idDevice, const bool deviceArchSupported,
												const bool dataFitsMemoryDevice, const bool localReference, const bool localIndex,
												const size_t recomendedMemorySize, const size_t requiredMemorySize)
{
	struct cudaDeviceProp devProp;
	const size_t memoryFree = gpu_get_device_free_memory(idDevice);
	CUDA_ERROR(cudaGetDeviceProperties(&devProp, idDevice));
	fprintf(stderr, "GPU Device %d: %s", idDevice, devProp.name);

	if(!deviceArchSupported){
		fprintf(stderr, "\t \t UNSUPPORTED Architecture [Compute Capability < 2.0: UNSUITABLE] \n");
		return(SUCCESS);
	}
	if(!dataFitsMemoryDevice){
		fprintf(stderr, "\t \t INSUFFICIENT DEVICE MEMORY (Mem Required: %lu MBytes - Mem Available: %lu MBytes) \n",
				GPU_CONVERT_B_TO_MB(requiredMemorySize), GPU_CONVERT_B_TO_MB(memoryFree));
		return(SUCCESS);
	}
	if((deviceArchSupported && dataFitsMemoryDevice)){
		fprintf(stderr, "\t \t Selected Device [RUNNING] \n");
		if (!localReference || !localIndex){
			fprintf(stderr, "\t \t WARNING PERF. SLOWDOWNS (Mem Recommended: %lu MBytes - Mem Available: %lu MBytes) \n",
					GPU_CONVERT_B_TO_MB(recomendedMemorySize), GPU_CONVERT_B_TO_MB(memoryFree));
		}
		return(SUCCESS);
	}
	return(SUCCESS);
}


/************************************************************
Functions to free all the buffer resources (HOST & DEVICE)
************************************************************/

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

GPU_INLINE gpu_error_t gpu_free_devices_info(gpu_device_info_t **devices)
{
	uint32_t idSupportedDevice, numSupportedDevices;

	numSupportedDevices = devices[0]->numSupportedDevices;
	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		if(devices[idSupportedDevice] != NULL){
			free(devices[idSupportedDevice]);
			devices[idSupportedDevice] = NULL;
		}
	}

    GPU_ERROR(gpu_free_devices_list_host(&devices));

    return(SUCCESS);
}
