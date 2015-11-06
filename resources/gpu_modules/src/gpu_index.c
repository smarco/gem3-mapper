#include "../include/gpu_index.h"

/************************************************************
Functions to initialize the index data on the DEVICE
************************************************************/

GPU_INLINE gpu_error_t gpu_load_index_PROFILE(const char *fn, gpu_index_buffer_t *index)
{
	FILE *fp = NULL;
	size_t result;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (E_OPENING_FILE);

    result = fread(&index->numEntries, sizeof(uint64_t), 1, fp);
	if (result != 1) return (E_READING_FILE);
    result = fread(&index->bwtSize, sizeof(uint64_t), 1, fp);
	if (result != 1) return (E_READING_FILE);

	CUDA_ERROR(cudaHostAlloc((void**) &index->h_fmi, index->numEntries * sizeof(gpu_fmi_entry_t), cudaHostAllocMapped));

	result = fread(index->h_fmi, sizeof(gpu_fmi_entry_t), index->numEntries, fp);
	if (result != index->numEntries) return (E_READING_FILE);

	fclose(fp);
	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_transfer_index_CPU_to_GPUs(gpu_index_buffer_t *index, gpu_device_info_t **devices)
{
	uint32_t deviceFreeMemory, idSupportedDevice;
	uint32_t numSupportedDevices = devices[0]->numSupportedDevices;

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		if(index[idSupportedDevice]->memorySpace == GPU_DEVICE_MAPPED){
			const size_t cpySize = index->numEntries * sizeof(gpu_fmi_entry_t);
			deviceFreeMemory = gpu_get_device_free_memory(devices[idSupportedDevice]->idDevice);
			if ((GPU_CONVERT_B_TO_MB(cpySize)) > deviceFreeMemory) return(E_INSUFFICIENT_MEM_GPU);
	    	CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
			//Synchronous allocate & transfer the FM-index to the GPU
			CUDA_ERROR(cudaMalloc((void**) &index->d_fmi[idSupportedDevice], cpySize));
			CUDA_ERROR(cudaMemcpy(index->d_fmi[idSupportedDevice], index->h_fmi, cpySize, cudaMemcpyHostToDevice));
		}else{
			index->d_fmi[idSupportedDevice] = index->h_fmi;
		}
	}

	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_init_index(gpu_index_buffer_t **index, const char *indexRaw,
									  const uint64_t bwtSize, const gpu_index_coding_t indexCoding,
									  const uint32_t numSupportedDevices)
{
	gpu_index_buffer_t *fmi = (gpu_index_buffer_t *) malloc(sizeof(gpu_index_buffer_t));

	fmi->d_fmi  	 = NULL;
	fmi->h_fmi 		 = NULL;
	fmi->memorySpace = NULL;
	fmi->bwtSize 	 = bwtSize;
	fmi->numEntries  = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;

	fmi->d_fmi = (gpu_fmi_entry_t **) malloc(numSupportedDevices * sizeof(gpu_fmi_entry_t *));
	if (fmi->d_fmi == NULL) GPU_ERROR(E_ALLOCATE_MEM);

	fmi->memorySpace = (gpu_data_location_t *) malloc(numSupportedDevices * sizeof(gpu_data_location_t));
	if (fmi->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);

	switch(indexCoding){
    	case GPU_INDEX_ASCII:
    		//GPU_ERROR(gpu_transform_index_ASCII(indexRaw, fmi));
    		GPU_ERROR(E_NOT_IMPLEMENTED);
    		break;
    	case GPU_INDEX_GEM_FULL:
    		//GPU_ERROR(gpu_transform_index_GEMFR(indexRaw, fmi));
    		GPU_ERROR(E_NOT_IMPLEMENTED);
    		break;
    	case GPU_INDEX_GEM_ONLY_FORWARD:
    		//GPU_ERROR(gpu_transform_index_GEMF(indexRaw, fmi));
    		GPU_ERROR(E_NOT_IMPLEMENTED);
    		break;
    	case GPU_INDEX_MFASTA_FILE:
    		//GPU_ERROR(gpu_load_index_MFASTA(indexRaw, fmi));
    		GPU_ERROR(E_NOT_IMPLEMENTED);
    		break;
    	case GPU_INDEX_PROFILE_FILE:
    		GPU_ERROR(gpu_load_index_PROFILE(indexRaw, fmi));
    		break;
    	default:
    		GPU_ERROR(E_INDEX_CODING);
    	  break;
	}

	(* index) = fmi;
	return (SUCCESS);
}



/************************************************************
 Functions to release the index data from the DEVICE & HOST
************************************************************/

GPU_INLINE gpu_error_t gpu_free_index_host(gpu_index_buffer_t *index)
{
    if(index->h_fmi != NULL){
        CUDA_ERROR(cudaFreeHost(index->h_fmi));
        index->h_fmi = NULL;
    }

    return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_free_unused_index_host(gpu_index_buffer_t *index, gpu_device_info_t **devices)
{
	uint32_t idSupportedDevice, numSupportedDevices;
	bool indexInHostSideUsed = false;

	numSupportedDevices = devices[0]->numSupportedDevices;

	//Free all the unused references in the host side
    for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
			if(index->memorySpace[idSupportedDevice] == GPU_HOST_MAPPED) indexInHostSideUsed = true;
    }

    if(!referenceInHostSideUsed){
    	GPU_ERROR(gpu_free_index_host(index));
    }

    return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_free_index_device(gpu_index_buffer_t *index, gpu_device_info_t **devices)
{
	const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
	uint32_t idSupportedDevice;

	//Free all the references in the devices
    for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
	    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
        if(index->d_fmi[idSupportedDevice] != NULL){
			if(index->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED)
            	CUDA_ERROR(cudaFree(index->d_fmi[idSupportedDevice]));
			index->d_fmi[idSupportedDevice] = NULL;
        }
    }

    //Free the index list
    if(index->d_fmi != NULL){
        free(index->d_fmi);
        index->d_fmi = NULL;
    }

    return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_free_index(gpu_index_buffer_t **index, gpu_device_info_t **devices)
{
	gpu_index_buffer_t *fmi = (* index);

    GPU_ERROR(gpu_free_index_host(fmi));
    GPU_ERROR(gpu_free_index_device(fmi, devices));

    if(fmi != NULL){
        free(fmi);
        fmi = NULL;
    }

    (* index) = fmi;
    return(SUCCESS);
}
