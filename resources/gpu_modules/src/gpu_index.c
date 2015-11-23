#include "../include/gpu_index.h"

/************************************************************
Functions to initialize the index data on the DEVICE
************************************************************/

GPU_INLINE uint32_t gpu_char_to_bin(const char base)
{
	uint32_t indexBase = GPU_ENC_DNA_CHAR_X;
	indexBase = ((base =='A') || (base =='a')) ? GPU_ENC_DNA_CHAR_A : indexBase;
	indexBase = ((base =='C') || (base =='c')) ? GPU_ENC_DNA_CHAR_C : indexBase;
	indexBase = ((base =='G') || (base =='g')) ? GPU_ENC_DNA_CHAR_G : indexBase;
	indexBase = ((base =='T') || (base =='t')) ? GPU_ENC_DNA_CHAR_T : indexBase;
	return(indexBase);
}

GPU_INLINE char gpu_bin_to_char(const uint32_t indexBase)
{
	const char LUT[5] = {'A', 'C', 'G', 'T', 'N'};
	const char base = (indexBase < 5) ? LUT[indexBase] : 'N';
	return (base);
}

GPU_INLINE void gpu_encode_entry_BWT_to_PEQ(gpu_index_bitmap_entry_t *bwt_entry, char base, const uint32_t size)
{
	const uint32_t binBase = gpu_char_to_bin(base);
	uint32_t idBase;
	for(idBase = 0; idBase < GPU_FMI_BWT_CHAR_LENGTH; ++idBase){
		bwt_entry->bitmaps[idBase] |= ((binBase >> idBase) & 0x1) << (GPU_UINT32_LENGTH - size - 1);
	}
}

GPU_INLINE void gpu_bining_bases(gpu_index_counter_entry_t *counterEntry, const char base)
{
	const uint32_t indexBase = gpu_char_to_bin(base);
	if(indexBase < 4) counterEntry->counters[indexBase]++;
}

GPU_INLINE gpu_error_t gpu_index_build_PEQ(gpu_index_buffer_t *fmi, const char *h_ascii_BWT, gpu_index_bitmap_entry_t *h_bitmap_BWT)
{
	uint64_t idEntry, i, bwtPosition;
	unsigned char bwtChar;

	//Padded to FMI_ENTRY_SIZE (128 bases)
	const uint32_t bwtNumEntries = fmi->numEntries * (GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH);

	for(idEntry = 0; idEntry < bwtNumEntries; ++idEntry){
		for(i = 0; i < GPU_UINT32_LENGTH; ++i){
			bwtPosition = (idEntry * GPU_UINT32_LENGTH) + i;
			if (bwtPosition < fmi->bwtSize) bwtChar = h_ascii_BWT[bwtPosition];
				else bwtChar = 'N'; //filling BWT padding
			gpu_encode_entry_BWT_to_PEQ(&h_bitmap_BWT[idEntry], bwtChar, i);
		}
	}
	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_index_build_counters(gpu_index_buffer_t *fmi, gpu_index_counter_entry_t *h_counters_FMI, const char *h_ascii_BWT)
{
	uint64_t idEntry, i, bwtPosition;
	uint32_t idBase, idCounter;
	char bwtChar;

	gpu_index_counter_entry_t  localCounters;
	const uint32_t  		   countersNumEntries = fmi->numEntries * (GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH);

	// Initialize first local BWT entry
	for(idEntry = 0; idEntry < countersNumEntries; ++idEntry){
		for(idBase = 0; idBase < GPU_FMI_NUM_COUNTERS; ++idBase){
			h_counters_FMI[idEntry].counters[idBase] = 0;
		}
	}

	// Accumulate values locally for each BWT entry
	for(idEntry = 0; idEntry < countersNumEntries - 1; ++idEntry){
		for(i = 0; i < GPU_UINT32_LENGTH; ++i){
			bwtPosition = (idEntry * GPU_UINT32_LENGTH) + i;
			if (bwtPosition < fmi->bwtSize) bwtChar = h_ascii_BWT[bwtPosition];
				else bwtChar = 'N'; //filling BWT padding (will be N)
			gpu_bining_bases(&h_counters_FMI[idEntry + 1], bwtChar);
		}
	}

	// Accumulate values globally for all the BWT
	for(idEntry = 1; idEntry < countersNumEntries; ++idEntry){
		for(idBase = 0; idBase < GPU_FMI_NUM_COUNTERS; ++idBase){
			h_counters_FMI[idEntry].counters[idBase] += h_counters_FMI[idEntry - 1].counters[idBase];
		}
	}

	// Prepare the local counters with the previous accumulative letters
	localCounters.counters[0] = 0;
	for(idBase = 1; idBase < GPU_FMI_NUM_COUNTERS; ++idBase){
		localCounters.counters[idBase] = localCounters.counters[idBase - 1] + h_counters_FMI[countersNumEntries - 4].counters[idBase - 1];
	}

	// Accumulate the previous alphabet letters to the global counters
	for(idEntry = 0; idEntry < countersNumEntries; ++idEntry){
		for(idBase = 1; idBase < GPU_FMI_NUM_COUNTERS; ++idBase){
			h_counters_FMI[idEntry].counters[idBase] += localCounters.counters[idBase];
		}
	}

	return (SUCCESS);
}

GPU_INLINE void gpu_index_set_layout_counter(gpu_index_counter_entry_t *h_counters_FMI, gpu_fmi_entry_t *h_fmi, const uint64_t idEntry)
{
	uint32_t idCounter;
	for(idCounter = 0; idCounter < GPU_FMI_COUNTERS_PER_ENTRY; ++idCounter){
		const uint32_t FMI_SET_ALT_COUNTER = idEntry % GPU_FMI_COUNTERS_PER_ENTRY;
		h_fmi->counters[idCounter] = h_counters_FMI->counters[(FMI_SET_ALT_COUNTER * GPU_FMI_COUNTERS_PER_ENTRY) + idCounter];
	}
}

GPU_INLINE void gpu_index_set_layout_bitmap(gpu_index_bitmap_entry_t *h_bitmap_BWT, gpu_fmi_entry_t *h_fmi)
{
	const uint32_t LUT[12] = {3,7,11,0,1,2,4,5,6,8,9,10}; // 1 1 (1) 0 2 2 (2) 0 3 3 (3) (0)
	uint32_t idPacket, idFMIBucket, padding = 0;
	const uint32_t FMI_NUM_BITMAPS        = GPU_FMI_ENTRY_SIZE * GPU_FMI_BWT_CHAR_LENGTH / GPU_UINT32_LENGTH; // 12 words             (4 FMI entries x 3 bits)
	const uint32_t NUM_PACKET_BMP_ENTRIES = GPU_FMI_ENTRY_SIZE  / GPU_UINT32_LENGTH;           				  // 4 BMP entries        (128 bases / 32 bits)
	const uint32_t NUM_PACKET_FMI_ENTRIES = FMI_NUM_BITMAPS / NUM_PACKET_BMP_ENTRIES;         				  // 3 FMI bitmap entries (12 bitmaps / 4 packets)

	for(idFMIBucket = 0; idFMIBucket < NUM_PACKET_BMP_ENTRIES; ++idFMIBucket)                 // Iterate over BMP entries (4)
		h_bitmap_BWT[idFMIBucket].bitmaps[NUM_PACKET_FMI_ENTRIES - 1] = ~ h_bitmap_BWT[idFMIBucket].bitmaps[NUM_PACKET_FMI_ENTRIES - 1];

		for(idPacket = 0; idPacket < NUM_PACKET_BMP_ENTRIES; ++idPacket){          		      // Iterate over BMP entries (4)
			for(idFMIBucket = 0; idFMIBucket < NUM_PACKET_FMI_ENTRIES; ++idFMIBucket){        // Iterate over FMI bitmaps (3)
				h_fmi->bitmaps[LUT[idPacket * NUM_PACKET_FMI_ENTRIES + idFMIBucket]] =  h_bitmap_BWT[idPacket].bitmaps[idFMIBucket];
		}
	}
}


GPU_INLINE gpu_error_t gpu_index_build_FMI(gpu_index_buffer_t *fmi, gpu_index_bitmap_entry_t *h_bitmap_BWT, gpu_index_counter_entry_t *h_counters_FMI)
{
	const uint32_t BITMAPS_PER_FMI = GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH; // 4 FMI bitmap entries
	uint64_t idEntry;

	for(idEntry = 0; idEntry < fmi->numEntries; ++idEntry){
		gpu_index_set_layout_counter(&h_counters_FMI[idEntry * BITMAPS_PER_FMI], &fmi->h_fmi[idEntry], idEntry);
		gpu_index_set_layout_bitmap(&h_bitmap_BWT[idEntry * BITMAPS_PER_FMI], &fmi->h_fmi[idEntry]);
	}

	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_transform_index_ASCII(const char *h_BWT, gpu_index_buffer_t *fmi)
{
	gpu_index_bitmap_entry_t  *h_bitmaps_BWT  = NULL;
	gpu_index_counter_entry_t *h_counters_FMI = NULL;

	const uint32_t bwtNumEntries 	  = fmi->numEntries * (GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH);
	const uint32_t countersNumEntries = fmi->numEntries * (GPU_FMI_ENTRY_SIZE / GPU_UINT32_LENGTH);

	h_bitmaps_BWT = (gpu_index_bitmap_entry_t *) malloc(bwtNumEntries * sizeof(gpu_index_bitmap_entry_t));
	if (h_bitmaps_BWT == NULL) return (E_ALLOCATE_MEM);
	h_counters_FMI = (gpu_index_counter_entry_t *) malloc(countersNumEntries * sizeof(gpu_index_counter_entry_t));
	if (h_counters_FMI == NULL) return (E_ALLOCATE_MEM);
	CUDA_ERROR(cudaHostAlloc((void**) &fmi->h_fmi, fmi->numEntries * sizeof(gpu_fmi_entry_t), cudaHostAllocMapped));
	if (fmi->h_fmi == NULL) return (E_ALLOCATE_MEM);

	GPU_ERROR(gpu_index_build_PEQ(fmi, h_BWT, h_bitmaps_BWT));
	GPU_ERROR(gpu_index_build_counters(fmi, h_counters_FMI, h_BWT));
	GPU_ERROR(gpu_index_build_FMI(fmi, h_bitmaps_BWT, h_counters_FMI));

	free(h_bitmaps_BWT);
	free(h_counters_FMI);
	return(SUCCESS);
}

GPU_INLINE gpu_error_t gpu_transform_index_GEM_FULL(gpu_gem_fmi_dto_t* const gpu_gem_fmi_dto, gpu_index_buffer_t *fmi)
{
  // BWT Parameters
  const uint64_t BWT_MINOR_BLOCKS_PER_MAYOR_BLOCK = (1<<10); /* 1024 */
  const uint64_t BWT_MINOR_BLOCK_LENGTH = 64;
  const uint64_t* c = gpu_gem_fmi_dto->c;
  const uint64_t* C = gpu_gem_fmi_dto->C;
  const uint64_t* mayor_counters = gpu_gem_fmi_dto->mayor_counters;
  const uint32_t* bwt_mem = (uint32_t*) gpu_gem_fmi_dto->bwt_mem;
  const uint64_t bwt_length = gpu_gem_fmi_dto->bwt_length;
  // Allocate fmi memory
  CUDA_ERROR(cudaHostAlloc((void**) &fmi->h_fmi, fmi->numEntries * sizeof(gpu_fmi_entry_t), cudaHostAllocMapped));
  if (fmi->h_fmi == NULL) return (E_ALLOCATE_MEM);
  gpu_fmi_entry_t* h_fmi = fmi->h_fmi; // Host FMI
  // Iterate thought the BWT memory layout
  uint16_t* minor_counters;
  uint32_t x0, x1, x2;
  uint32_t y0, y1, y2;
  uint32_t z0, z1, z2;
  uint32_t w0, w1, w2;
  uint64_t minor_block = 0, h_fmi_entry = 0;
  uint64_t bwt_pos, i;
  for (bwt_pos=0;bwt_pos<bwt_length;) {
    // Get next block (Block0)
    minor_counters = (uint16_t*)bwt_mem; bwt_mem += 4; // Minor Counters
    x0 = *(bwt_mem); ++bwt_mem;
    y0 = *(bwt_mem); ++bwt_mem;
    x1 = *(bwt_mem); ++bwt_mem;
    y1 = *(bwt_mem); ++bwt_mem;
    x2 = *(bwt_mem); ++bwt_mem;
    y2 = *(bwt_mem); ++bwt_mem;
    bwt_mem += 2; // Sampling bitmap
    bwt_pos+=BWT_MINOR_BLOCK_LENGTH; // Next Block
    // Check end-of-bwt
    if (bwt_pos<bwt_length) {
      // Get next block (Block1)
      bwt_mem += 4; // Skip Minor Counters
      z0 = *(bwt_mem); ++bwt_mem;
      w0 = *(bwt_mem); ++bwt_mem;
      z1 = *(bwt_mem); ++bwt_mem;
      w1 = *(bwt_mem); ++bwt_mem;
      z2 = *(bwt_mem); ++bwt_mem;
      w2 = *(bwt_mem); ++bwt_mem;
      bwt_mem += 2; // Sampling bitmap
      bwt_pos+=BWT_MINOR_BLOCK_LENGTH; // Next Block
    } else {
      // Padding zeros
      z0 = 0; z1 = 0; z2 = 0;
      w0 = 0; w1 = 0; w2 = 0;
    }
    // Write Counters
    if (h_fmi_entry % 2 == 0) {
      h_fmi->counters[0] = minor_counters[0] + mayor_counters[0]; // 'A'
      h_fmi->counters[1] = minor_counters[1] + mayor_counters[1]; // 'C'
    } else {
      h_fmi->counters[0] = minor_counters[2] + mayor_counters[2]; // 'G'
      h_fmi->counters[1] = minor_counters[3] + mayor_counters[3]; // 'T'
    }
    // Write Bitmap
    h_fmi->bitmaps[0] = __builtin_bswap32(y0);
    h_fmi->bitmaps[1] = __builtin_bswap32(y1);
    h_fmi->bitmaps[2] = __builtin_bswap32(~y2);
    h_fmi->bitmaps[3] = __builtin_bswap32(x0);
    h_fmi->bitmaps[4] = __builtin_bswap32(z0);
    h_fmi->bitmaps[5] = __builtin_bswap32(z1);
    h_fmi->bitmaps[6] = __builtin_bswap32(~z2);
    h_fmi->bitmaps[7] = __builtin_bswap32(x1);
    h_fmi->bitmaps[8] = __builtin_bswap32(w0);
    h_fmi->bitmaps[9] = __builtin_bswap32(w1);
    h_fmi->bitmaps[10] = __builtin_bswap32(~w2);
    h_fmi->bitmaps[11] = __builtin_bswap32(~x2);
    // Next Entry
    ++h_fmi; ++h_fmi_entry;
    minor_block += 2;
    if (minor_block >= BWT_MINOR_BLOCKS_PER_MAYOR_BLOCK) {
      mayor_counters += 8;
    }
  }
  // Ghost entry for alternate counters
  if (h_fmi_entry % 2 == 1) {
    // Write Counters
    h_fmi->counters[0] = c[2] + C[2]; // 'G'
    h_fmi->counters[1] = c[3] + C[3]; // 'T'
    // Write Bitmap
    h_fmi->bitmaps[0] = 0; h_fmi->bitmaps[1] = 0;
    h_fmi->bitmaps[2] = 0; h_fmi->bitmaps[3] = 0;
    h_fmi->bitmaps[4] = 0; h_fmi->bitmaps[5] = 0;
    h_fmi->bitmaps[6] = 0; h_fmi->bitmaps[7] = 0;
    h_fmi->bitmaps[8] = 0; h_fmi->bitmaps[9] = 0;
    h_fmi->bitmaps[10] = 0; h_fmi->bitmaps[11] = 0;
  }
  // Return SUCCESS
  return (SUCCESS);
}
GPU_INLINE gpu_error_t gpu_load_BWT_MFASTA(const char *fn, gpu_index_buffer_t *fmi, char **h_BWT)
{
	FILE *fp = NULL;
	char *h_ascii_BWT = NULL;
	char lineFile[GPU_FILE_SIZE_LINES];
	uint64_t sizeFile = 0, position = 0;
	int32_t charsRead = 0;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (E_OPENING_FILE);

	fseek(fp, 0L, SEEK_END);
	sizeFile = ftell(fp);
	rewind(fp);

	h_ascii_BWT = (char*) malloc(sizeFile * sizeof(char));
	if (h_ascii_BWT == NULL) return (E_ALLOCATE_MEM);

	while((!feof(fp)) && (fgets(lineFile, GPU_FILE_SIZE_LINES, fp) != NULL)){
		if (lineFile[0] != '>'){
			charsRead = strlen(lineFile);
			if(charsRead) charsRead--;
			memcpy((h_ascii_BWT + position), lineFile, charsRead);
			position +=  charsRead;
		}
	}

	fmi->bwtSize    = position;
	fmi->numEntries = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
	(* h_BWT) 		= h_ascii_BWT;

	fclose(fp);
	return (SUCCESS);
}

GPU_INLINE gpu_error_t gpu_save_index_PROFILE(const char *fn, gpu_index_buffer_t *index)
{
	const uint32_t sizeFileName = 512;
    char fileName[sizeFileName];
    FILE *fp = NULL;

    sprintf(fileName, "%s.%lu.%u.fmi", fn, index->bwtSize, GPU_FMI_ENTRY_SIZE);
    fp = fopen(fileName, "wb");
    if (fp == NULL) return (E_WRITING_FILE);

    fwrite(&index->numEntries, sizeof(uint64_t), 1, fp);
    fwrite(&index->bwtSize, sizeof(uint64_t), 1, fp);
    fwrite(index->h_fmi, sizeof(gpu_fmi_entry_t), index->numEntries, fp);

    fclose(fp);
    return (SUCCESS);
}

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
		if(index->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
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

GPU_INLINE gpu_error_t gpu_init_index(gpu_index_buffer_t **index, const void *indexRaw,
									  const uint64_t bwtSize, const gpu_index_coding_t indexCoding,
									  const uint32_t numSupportedDevices, gpu_module_t activeModules)
{
	gpu_index_buffer_t *fmi = (gpu_index_buffer_t *) malloc(sizeof(gpu_index_buffer_t));
	uint32_t idSupDevice;

	fmi->d_fmi  		 = NULL;
	fmi->h_fmi 		 	 = NULL;
	fmi->memorySpace 	 = NULL;
	fmi->bwtSize 	 	 = 0;
	fmi->numEntries  	 = 0;
	fmi->activeModules	 = activeModules;

	fmi->d_fmi = (gpu_fmi_entry_t **) malloc(numSupportedDevices * sizeof(gpu_fmi_entry_t *));
	if (fmi->d_fmi == NULL) GPU_ERROR(E_ALLOCATE_MEM);
	fmi->memorySpace = (memory_alloc_t *) malloc(numSupportedDevices * sizeof(memory_alloc_t));
	if (fmi->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);

	for(idSupDevice = 0; idSupDevice < numSupportedDevices; ++idSupDevice){
		fmi->d_fmi[idSupDevice] 	  = NULL;
		fmi->memorySpace[idSupDevice] = GPU_NONE_MAPPED;
	}

	if(activeModules & (GPU_FMI_DECODE_POS | GPU_FMI_EXACT_SEARCH)){
		fmi->bwtSize 	 	 = bwtSize;
		fmi->numEntries  	 = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
		char *h_BWT = NULL;

		switch(indexCoding){
			case GPU_INDEX_ASCII:
				GPU_ERROR(gpu_transform_index_ASCII((const char*)indexRaw, fmi));
				break;
			case GPU_INDEX_GEM_FULL:
				GPU_ERROR(gpu_transform_index_GEM_FULL((gpu_gem_fmi_dto_t*)indexRaw, fmi));
				break;
			case GPU_INDEX_MFASTA_FILE:
				GPU_ERROR(gpu_load_BWT_MFASTA((const char*)indexRaw, fmi, &h_BWT));
				GPU_ERROR(gpu_transform_index_ASCII(h_BWT, fmi));
				free(h_BWT);
				break;
			case GPU_INDEX_PROFILE_FILE:
				GPU_ERROR(gpu_load_index_PROFILE((const char*)indexRaw, fmi));
				break;
			default:
				GPU_ERROR(E_INDEX_CODING);
			break;
		}
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
	const gpu_module_t activeModules = index->activeModules;
	uint32_t idSupportedDevice, numSupportedDevices;
	bool indexInHostSideUsed = false;

	if(activeModules & (GPU_FMI_DECODE_POS | GPU_FMI_EXACT_SEARCH)){
		numSupportedDevices = devices[0]->numSupportedDevices;
		//Free all the unused references in the host side
		for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
				if(index->memorySpace[idSupportedDevice] == GPU_HOST_MAPPED) indexInHostSideUsed = true;
		}

		if(!indexInHostSideUsed){
			GPU_ERROR(gpu_free_index_host(index));
		}
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
	const gpu_module_t activeModules = fmi->activeModules;

	if(activeModules & (GPU_FMI_EXACT_SEARCH | GPU_FMI_DECODE_POS)){
		GPU_ERROR(gpu_free_index_host(fmi));
		GPU_ERROR(gpu_free_index_device(fmi, devices));
	}

    if(fmi != NULL){
    	free(fmi->memorySpace);
        free(fmi);
        fmi = NULL;
    }

    (* index) = fmi;
    return(SUCCESS);
}
