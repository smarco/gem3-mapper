
#include "../include/gpu_reference.h"

/************************************************************
String basic functions
************************************************************/

GPU_INLINE uint64_t gpu_char_to_bin_ASCII(unsigned char base)
{
	switch(base)
	{
    	case 'A':
    	case 'a':
    	    return(GPU_ENC_DNA_CHAR_A);
    	case 'C':
    	case 'c':
    	    return(GPU_ENC_DNA_CHAR_C << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH));
    	case 'G':
    	case 'g':
    	    return(GPU_ENC_DNA_CHAR_G << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH));
    	case 'T':
    	case 't':
    	    return(GPU_ENC_DNA_CHAR_T << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH));
    	default :
    	    return(GPU_ENC_DNA_CHAR_N << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH));
	}
}

GPU_INLINE char gpu_complement_base(const char character)
{
	char referenceChar = character;
	referenceChar = (character == GPU_ENC_DNA_CHAR_A) ? GPU_ENC_DNA_CHAR_T : referenceChar;
	referenceChar = (character == GPU_ENC_DNA_CHAR_C) ? GPU_ENC_DNA_CHAR_G : referenceChar;
	referenceChar = (character == GPU_ENC_DNA_CHAR_G) ? GPU_ENC_DNA_CHAR_C : referenceChar;
	referenceChar = (character == GPU_ENC_DNA_CHAR_T) ? GPU_ENC_DNA_CHAR_A : referenceChar;
	return(referenceChar);
}


/************************************************************
Transform reference functions
************************************************************/

GPU_INLINE gpu_error_t gpu_transform_reference_ASCII(const char *referenceASCII, gpu_reference_buffer_t *reference)
{
	uint64_t indexBase, bitmap;
	uint64_t idEntry, i, referencePosition;
	unsigned char referenceChar;

	CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));

	for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
		bitmap = 0;
		for(i = 0; i < GPU_REFERENCE_CHARS_PER_ENTRY; i++){
			referencePosition = idEntry * GPU_REFERENCE_CHARS_PER_ENTRY + i;
			if (referencePosition < reference->size) referenceChar = referenceASCII[referencePosition];
				else referenceChar = 'N'; //filling reference padding
			indexBase = gpu_char_to_bin_ASCII(referenceChar);
			bitmap = (bitmap >> GPU_REFERENCE_CHAR_LENGTH) | indexBase;
		}
		reference->h_reference[referencePosition / GPU_REFERENCE_CHARS_PER_ENTRY] = bitmap;
	}
	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_transform_reference_GEM_FR(const char *referenceGEM, gpu_reference_buffer_t *reference)
{
	uint64_t indexBase, bitmap;
	uint64_t idEntry, i, referencePosition;
	unsigned char referenceChar;

	CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));

	for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
		bitmap = 0;
		for(i = 0; i < GPU_REFERENCE_CHARS_PER_ENTRY; ++i){
			referencePosition = idEntry * GPU_REFERENCE_CHARS_PER_ENTRY + i;
			if (referencePosition < reference->size) referenceChar = referenceGEM[referencePosition];
				else referenceChar = 'N'; //filling reference padding
			indexBase = ((uint64_t) referenceChar) << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH);
			bitmap = (bitmap >> GPU_REFERENCE_CHAR_LENGTH) | indexBase;
		}
		reference->h_reference[referencePosition / GPU_REFERENCE_CHARS_PER_ENTRY] = bitmap;
	}
	return(SUCCESS);
}


GPU_INLINE gpu_error_t gpu_transform_reference_GEM_F(const char *referenceGEM, gpu_reference_buffer_t *reference)
{
	uint64_t indexBase, bitmap;
	uint64_t idEntry, i, referencePosition;
	char referenceChar;

	// Recompute size of the full reference (forward + reverse-complement)
	const uint64_t forward_ref_size = reference->size;
	const uint64_t total_ref_size = 2 * forward_ref_size;

	reference->size = total_ref_size;
	reference->numEntries = GPU_DIV_CEIL(total_ref_size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

	// Allocate CUDA-HostMem
	CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));

	// Copy reference
	for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
		bitmap = 0;
		for(i = 0; i < GPU_REFERENCE_CHARS_PER_ENTRY; ++i){
			referencePosition = idEntry * GPU_REFERENCE_CHARS_PER_ENTRY + i;
			if (referencePosition < forward_ref_size) {
				referenceChar = referenceGEM[referencePosition];
			} else if (referencePosition < reference->size) {
				referenceChar = gpu_complement_base(referenceGEM[2*forward_ref_size-referencePosition-2]);
			} else {
				referenceChar = 'N'; //filling reference padding
			}
			indexBase = ((uint64_t) referenceChar) << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH);
			bitmap = (bitmap >> GPU_REFERENCE_CHAR_LENGTH) | indexBase;
		}
		reference->h_reference[referencePosition / GPU_REFERENCE_CHARS_PER_ENTRY] = bitmap;
	}

	// Return
	return(SUCCESS);
}


/************************************************************
Input & Output reference functions
************************************************************/

GPU_INLINE gpu_error_t gpu_load_reference_MFASTA(const char *fn, gpu_reference_buffer_t *reference)
{
	FILE *fp = NULL;
	char lineFile[GPU_FILE_SIZE_LINES], *tmp_reference;
	uint64_t sizeFile = 0, position = 0;
	int32_t charsRead = 0;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (E_OPENING_FILE);

	fseek(fp, 0L, SEEK_END);
	sizeFile = ftell(fp);
	rewind(fp);

	tmp_reference = (char*) malloc(sizeFile * sizeof(char));
	if (tmp_reference == NULL) return (E_ALLOCATE_MEM);

	if ((fgets(lineFile, GPU_FILE_SIZE_LINES, fp) == NULL) || (lineFile[0] != '>'))
		return (E_READING_FILE);

	while((!feof(fp)) && (fgets(lineFile, GPU_FILE_SIZE_LINES, fp) != NULL)){
		if (lineFile[0] != '>'){
			charsRead = strlen(lineFile);
			if(charsRead) charsRead--;
			memcpy((tmp_reference + position), lineFile, charsRead);
			position +=  charsRead;
		}
	}

	reference->size = position;
	reference->numEntries = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
	GPU_ERROR(gpu_transform_reference_ASCII(tmp_reference, reference));

	fclose(fp);
	free(tmp_reference);
	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_load_reference_PROFILE(const char *fn, gpu_reference_buffer_t *reference)
{
	FILE *fp = NULL;
	size_t result;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (E_OPENING_FILE);

    result = fread(&reference->numEntries, sizeof(uint64_t), 1, fp);
	if (result != 1) return (E_READING_FILE);
    result = fread(&reference->size, sizeof(uint64_t), 1, fp);
	if (result != 1) return (E_READING_FILE);

	CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));

	result = fread(reference->h_reference, sizeof(uint64_t), reference->numEntries, fp);
	if (result != reference->numEntries) return (E_READING_FILE);

	fclose(fp);
	return (SUCCESS);
}



/************************************************************
Functions to get the GPU FMI buffers
************************************************************/

GPU_INLINE gpu_error_t gpu_init_reference(gpu_reference_buffer_t **reference, const char *referenceRaw,
										  const uint64_t refSize, const gpu_ref_coding_t refCoding,
										  const uint32_t numSupportedDevices, gpu_module_t activeModules)
{
	gpu_reference_buffer_t *ref = (gpu_reference_buffer_t *) malloc(sizeof(gpu_reference_buffer_t));
	uint32_t idSupDevice;

	ref->d_reference   = NULL;
	ref->h_reference   = NULL;
	ref->memorySpace   = NULL;
	ref->size 		   = 0;
	ref->numEntries    = GPU_DIV_CEIL(ref->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
	ref->activeModules = activeModules;

	ref->d_reference = (uint64_t **) malloc(numSupportedDevices * sizeof(uint64_t *));
	if (ref->d_reference == NULL) GPU_ERROR(E_ALLOCATE_MEM);
	ref->memorySpace = (memory_alloc_t *) malloc(numSupportedDevices * sizeof(memory_alloc_t));
	if (ref->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);

	for(idSupDevice = 0; idSupDevice < numSupportedDevices; ++idSupDevice){
		ref->d_reference[idSupDevice] = NULL;
		ref->memorySpace[idSupDevice] = GPU_NONE_MAPPED;
	}

	if(activeModules & GPU_BPM){
		ref->size = refSize;
		ref->numEntries = GPU_DIV_CEIL(ref->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

		switch(refCoding){
			case GPU_REF_ASCII:
				GPU_ERROR(gpu_transform_reference_ASCII(referenceRaw, ref));
				break;
			case GPU_REF_GEM_FULL:
				GPU_ERROR(gpu_transform_reference_GEM_FR(referenceRaw, ref));
				break;
			case GPU_REF_GEM_ONLY_FORWARD:
				GPU_ERROR(gpu_transform_reference_GEM_F(referenceRaw, ref));
				break;
			case GPU_REF_MFASTA_FILE:
				GPU_ERROR(gpu_load_reference_MFASTA(referenceRaw, ref));
				break;
			case GPU_REF_PROFILE_FILE:
				GPU_ERROR(gpu_load_reference_PROFILE(referenceRaw, ref));
				break;
			default:
				GPU_ERROR(E_REFERENCE_CODING);
			  break;
		}
	}
	(* reference) = ref;
	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_transfer_reference_CPU_to_GPUs(gpu_reference_buffer_t *reference, gpu_device_info_t **devices)
{
	uint32_t deviceFreeMemory, idSupportedDevice;
	uint32_t numSupportedDevices = devices[0]->numSupportedDevices;

	for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
		if(reference->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
			const size_t cpySize = reference->numEntries * sizeof(uint64_t);
			deviceFreeMemory = gpu_get_device_free_memory(devices[idSupportedDevice]->idDevice);
			if ((GPU_CONVERT_B_TO_MB(cpySize)) > deviceFreeMemory) return(E_INSUFFICIENT_MEM_GPU);
	    	CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
			//Synchronous allocate & transfer Binary Reference to GPU
			CUDA_ERROR(cudaMalloc((void**) &reference->d_reference[idSupportedDevice], cpySize));
			CUDA_ERROR(cudaMemcpy(reference->d_reference[idSupportedDevice], reference->h_reference, cpySize, cudaMemcpyHostToDevice));
		}else{
			reference->d_reference[idSupportedDevice] = reference->h_reference;
		}
	}

	return (SUCCESS);
}


/************************************************************
Free reference functions
************************************************************/

GPU_INLINE gpu_error_t gpu_free_reference_host(gpu_reference_buffer_t *reference)
{
    if(reference->h_reference != NULL){
        CUDA_ERROR(cudaFreeHost(reference->h_reference));
        reference->h_reference = NULL;
    }

    return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_free_unused_reference_host(gpu_reference_buffer_t *reference, gpu_device_info_t **devices)
{
	const gpu_module_t activeModules = reference->activeModules;
	uint32_t idSupportedDevice, numSupportedDevices;
	bool referenceInHostSideUsed = false;

	if(activeModules & GPU_BPM){
		numSupportedDevices = devices[0]->numSupportedDevices;
		//Free all the unused references in the host side
		for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
				if(reference->memorySpace[idSupportedDevice] == GPU_HOST_MAPPED) referenceInHostSideUsed = true;
		}

		if(!referenceInHostSideUsed){
			GPU_ERROR(gpu_free_reference_host(reference));
		}
	}

    return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_free_reference_device(gpu_reference_buffer_t *reference, gpu_device_info_t **devices)
{
	const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
	uint32_t idSupportedDevice;

	//Free all the references in the devices
    for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
	    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
        if(reference->d_reference[idSupportedDevice] != NULL){
			if(reference->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED)
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

GPU_INLINE gpu_error_t gpu_free_reference(gpu_reference_buffer_t **reference, gpu_device_info_t **devices)
{
	gpu_reference_buffer_t *ref = (* reference);
	const gpu_module_t activeModules = ref->activeModules;

	if(activeModules & GPU_BPM){
		GPU_ERROR(gpu_free_reference_host(ref));
    	GPU_ERROR(gpu_free_reference_device(ref, devices));
	}

    if(ref != NULL){
    	free(ref->memorySpace);
        free(ref);
        ref = NULL;
    }

    (* reference) = ref;
    return(SUCCESS);
}
