
#include "../include/gpu_reference.h"

/************************************************************
String basic functions
************************************************************/

uint64_t gpu_char_to_bin_ASCII(const unsigned char base)
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
      return(GPU_ENC_DNA_CHAR_A << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH));
  }
}

char gpu_complement_base(const char character)
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

gpu_error_t gpu_transform_reference_ASCII(const char* const referenceASCII, gpu_reference_buffer_t* const reference)
{
  uint64_t indexBase, bitmap;
  uint64_t idEntry, i, referencePosition;
  unsigned char referenceChar;

  reference->numEntries = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));

  for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
    bitmap = 0;
    for(i = 0; i < GPU_REFERENCE_CHARS_PER_ENTRY; ++i){
      referencePosition = idEntry * GPU_REFERENCE_CHARS_PER_ENTRY + i;
      if (referencePosition < reference->size) referenceChar = referenceASCII[referencePosition];
        else referenceChar = 'A'; //filling reference padding
      indexBase = gpu_char_to_bin_ASCII(referenceChar);
      bitmap = (bitmap >> GPU_REFERENCE_CHAR_LENGTH) | indexBase;
    }
    reference->h_reference[referencePosition / GPU_REFERENCE_CHARS_PER_ENTRY] = bitmap;
  }
  return(SUCCESS);
}

gpu_error_t gpu_transform_reference_GEM(const gpu_gem_ref_dto_t* const gem_reference, gpu_reference_buffer_t* const reference)
{
  uint64_t bitmap, base;
  uint64_t idEntry, i, referencePosition;

  // Compute size of the forward reference
  const char* const h_gem_reference = gem_reference->reference;
  const uint64_t total_ref_size     = gem_reference->ref_length;
  const uint64_t baseMask           = GPU_UINT64_ONES << GPU_REFERENCE_CHAR_LENGTH;

  reference->size = total_ref_size;
  reference->numEntries = GPU_DIV_CEIL(total_ref_size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

  CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));

  for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
    bitmap = 0;
    for(i = 0; i < GPU_REFERENCE_CHARS_PER_ENTRY; ++i){
      referencePosition = idEntry * GPU_REFERENCE_CHARS_PER_ENTRY + i;
      if (referencePosition < reference->size) base = (uint64_t) h_gem_reference[referencePosition];
        else base = GPU_ENC_DNA_CHAR_A; //filling reference padding
      base = (base & baseMask) ? GPU_ENC_DNA_CHAR_A : base;
      base = base << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH);
      bitmap = (bitmap >> GPU_REFERENCE_CHAR_LENGTH) | base;
    }
    reference->h_reference[referencePosition / GPU_REFERENCE_CHARS_PER_ENTRY] = bitmap;
  }
  return(SUCCESS);
}


gpu_error_t gpu_transform_reference_GEM_FULL(const gpu_gem_ref_dto_t* const gem_reference, gpu_reference_buffer_t* const reference)
{
  uint64_t base, bitmap;
  uint64_t idEntry, i, referencePosition;

  // Recompute size of the full reference (forward + reverse-complement)
  const char* const h_gem_reference = gem_reference->reference;
  const uint64_t forward_ref_size   = gem_reference->ref_length;
  const uint64_t total_ref_size     = 2 * forward_ref_size;
  const uint64_t baseMask           = GPU_UINT64_ONES << GPU_REFERENCE_CHAR_LENGTH;

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
        base = (uint64_t) h_gem_reference[referencePosition];
      } else if (referencePosition < reference->size) {
        base = (uint64_t) gpu_complement_base(h_gem_reference[2*forward_ref_size-referencePosition-2]);
      } else {
        base = GPU_ENC_DNA_CHAR_A; //filling reference padding
      }
      base = (base & baseMask) ? GPU_ENC_DNA_CHAR_A : base;
      base = base << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH);
      bitmap = (bitmap >> GPU_REFERENCE_CHAR_LENGTH) | base;
    }
    reference->h_reference[referencePosition / GPU_REFERENCE_CHARS_PER_ENTRY] = bitmap;
  }

  // Return
  return(SUCCESS);
}


/************************************************************
Input & Output reference functions
************************************************************/

gpu_error_t gpu_read_reference(FILE* fp, gpu_reference_buffer_t* const reference)
{
  size_t result;

  result = fread(&reference->numEntries, sizeof(uint64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&reference->size, sizeof(uint64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));

  result = fread(reference->h_reference, sizeof(uint64_t), reference->numEntries, fp);
  if (result != reference->numEntries) return (E_READING_FILE);

  return (SUCCESS);
}

gpu_error_t gpu_write_reference(FILE* fp, const gpu_reference_buffer_t* const reference)
{
  size_t result;

  result = fwrite(&reference->numEntries, sizeof(uint64_t), 1, fp);
  if (result != 1) return (E_WRITING_FILE);
  result = fwrite(&reference->size, sizeof(uint64_t), 1, fp);
  if (result != 1) return (E_WRITING_FILE);

  result = fwrite(reference->h_reference, sizeof(uint64_t), reference->numEntries, fp);
  if (result != reference->numEntries) return (E_WRITING_FILE);

  return (SUCCESS);
}


/************************************************************
Initialize reference functions
************************************************************/

gpu_error_t gpu_init_reference_dto(gpu_reference_buffer_t* const ref)
{
  //Initialize the reference structure
  ref->d_reference   = NULL;
  ref->h_reference   = NULL;
  ref->memorySpace   = NULL;
  ref->size          = 0;
  ref->numEntries    = 0;
  ref->activeModules = GPU_NONE_MODULES;
  return (SUCCESS);
}

gpu_error_t gpu_transform_reference(const char* const referenceRaw, gpu_reference_buffer_t* const ref, const gpu_ref_coding_t refCoding)
{
  if((ref->activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  switch(refCoding){
    case GPU_REF_ASCII:
      GPU_ERROR(gpu_transform_reference_ASCII(referenceRaw, ref));
      break;
    case GPU_REF_GEM_ONLY_FORWARD:
      GPU_ERROR(gpu_transform_reference_GEM((gpu_gem_ref_dto_t*)referenceRaw, ref));
      break;
    case GPU_REF_GEM_FULL:
      GPU_ERROR(gpu_transform_reference_GEM_FULL((gpu_gem_ref_dto_t*)referenceRaw, ref));
      break;
    case GPU_REF_GEM_FILE:
      GPU_ERROR(gpu_load_reference_GEM_FULL(referenceRaw, ref));
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

  return(SUCCESS);
}

gpu_error_t gpu_init_reference(gpu_reference_buffer_t **reference, const char* const referenceRaw,
                               const uint64_t refSize, const gpu_ref_coding_t refCoding,
                               const uint32_t numSupportedDevices, const gpu_module_t activeModules)
{
  gpu_reference_buffer_t* const ref = (gpu_reference_buffer_t *) malloc(sizeof(gpu_reference_buffer_t));
  uint32_t idSupDevice;

  ref->d_reference   = NULL;
  ref->h_reference   = NULL;
  ref->memorySpace   = NULL;
  ref->size          = refSize;
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

  GPU_ERROR(gpu_transform_reference(referenceRaw, ref, refCoding));

  (* reference) = ref;
  return (SUCCESS);
}

gpu_error_t gpu_transfer_reference_CPU_to_GPUs(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices)
{
  uint32_t deviceFreeMemory, idSupportedDevice;
  uint32_t numSupportedDevices = devices[0]->numSupportedDevices;

  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(reference->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
      const size_t cpySize = reference->numEntries * sizeof(uint64_t);
      deviceFreeMemory = gpu_get_device_free_memory(devices[idSupportedDevice]->idDevice);
      if ((GPU_CONVERT__B_TO_MB(cpySize)) > deviceFreeMemory) return(E_INSUFFICIENT_MEM_GPU);
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

gpu_error_t gpu_free_reference_host(gpu_reference_buffer_t* const reference)
{
    if(reference->h_reference != NULL){
      CUDA_ERROR(cudaFreeHost(reference->h_reference));
      reference->h_reference = NULL;
    }

    return(SUCCESS);
}

gpu_error_t gpu_free_unused_reference_host(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices)
{
  const gpu_module_t activeModules = reference->activeModules;
  uint32_t idSupportedDevice, numSupportedDevices;
  bool referenceInHostSideUsed = false;

  if(activeModules & GPU_REFERENCE){
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

gpu_error_t gpu_free_reference_device(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices)
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

gpu_error_t gpu_free_reference(gpu_reference_buffer_t **reference, gpu_device_info_t** const devices)
{
  gpu_reference_buffer_t* ref = (* reference);
  const gpu_module_t activeModules = ref->activeModules;

  if(activeModules & GPU_REFERENCE){
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
