/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_REFERENCE_C_
#define GPU_REFERENCE_C_

#include "../include/gpu_reference.h"
#include "../include/gpu_io.h"

/************************************************************
Get information functions
************************************************************/

gpu_error_t gpu_reference_get_size(gpu_reference_buffer_t* const reference, size_t* bytesPerReference)
{
  (* bytesPerReference) = reference->numEntries * GPU_REFERENCE_BYTES_PER_ENTRY;
  return (SUCCESS);
}


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

gpu_error_t gpu_reference_transform_ASCII(const char* const referenceASCII, gpu_reference_buffer_t* const reference)
{
  uint64_t indexBase, bitmap;
  uint64_t idEntry, i, referencePosition;
  unsigned char referenceChar;

  reference->numEntries = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

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

gpu_error_t gpu_reference_transform_GEM(const gpu_gem_ref_dto_t* const gem_reference, gpu_reference_buffer_t* const reference)
{
  uint64_t bitmap, base;
  uint64_t idEntry, i, referencePosition;

  // Compute size of the forward reference
  const char* const h_gem_reference = gem_reference->reference;
  const uint64_t total_ref_size     = gem_reference->ref_length;
  const uint64_t baseMask           = GPU_UINT64_ONES << GPU_REFERENCE_CHAR_LENGTH;

  reference->size = total_ref_size;
  reference->numEntries = GPU_DIV_CEIL(total_ref_size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

  for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
    bitmap = 0;
    for(i = 0; i < GPU_REFERENCE_CHARS_PER_ENTRY; ++i){
      referencePosition = idEntry * GPU_REFERENCE_CHARS_PER_ENTRY + i;
      if (referencePosition < reference->size) {
        base = (uint64_t) h_gem_reference[referencePosition];
      } else {
        base = GPU_ENC_DNA_CHAR_A; //filling reference padding
      }
      base = (base & baseMask) ? GPU_ENC_DNA_CHAR_A : base;
      base = base << (GPU_UINT64_LENGTH - GPU_REFERENCE_CHAR_LENGTH);
      bitmap = (bitmap >> GPU_REFERENCE_CHAR_LENGTH) | base;
    }
    reference->h_reference[referencePosition / GPU_REFERENCE_CHARS_PER_ENTRY] = bitmap;
  }
  return(SUCCESS);
}


gpu_error_t gpu_reference_transform_GEM_FULL(const gpu_gem_ref_dto_t* const gem_reference, gpu_reference_buffer_t* const reference)
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

gpu_error_t gpu_reference_read_specs(FILE* fp, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules)
{
  size_t result;

  if(activeModules & GPU_REFERENCE)
    return(E_MODULE_NOT_FOUND);

  result = fread(&reference->numEntries, sizeof(uint64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&reference->size, sizeof(uint64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);

  return (SUCCESS);
}

gpu_error_t gpu_reference_read(FILE* fp, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules)
{
  size_t result;

  if(activeModules & GPU_REFERENCE)
    return(E_MODULE_NOT_FOUND);

  result = fread(&reference->numEntries, sizeof(uint64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(&reference->size, sizeof(uint64_t), 1, fp);
  if (result != 1) return (E_READING_FILE);
  result = fread(reference->h_reference, sizeof(uint64_t), reference->numEntries, fp);
  if (result != reference->numEntries) return (E_READING_FILE);

  return (SUCCESS);
}

gpu_error_t gpu_reference_write(FILE* fp, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules)
{
  size_t result;

  if(activeModules & GPU_REFERENCE)
    return(E_MODULE_NOT_FOUND);

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

gpu_error_t gpu_reference_init_dto(gpu_reference_buffer_t* const ref)
{
  //Initialize the reference structure
  ref->d_reference    = NULL;
  ref->h_reference    = NULL;
  ref->hostAllocStats = GPU_PAGE_UNLOCKED;
  ref->memorySpace    = NULL;
  ref->size           = 0;
  ref->numEntries     = 0;
  ref->activeModules  = GPU_NONE_MODULES;

  return (SUCCESS);
}

gpu_error_t gpu_reference_set_specs(gpu_reference_buffer_t* const ref, const char* const referenceRaw,
                                    const gpu_ref_coding_t refCoding, const gpu_module_t activeModules)
{
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  switch(refCoding){
    case GPU_REF_ASCII:
      /* Not require special I/O initialization */
      break;
    case GPU_REF_GEM_FULL:
      /* Not require special I/O initialization */
      break;
    case GPU_REF_GEM_ONLY_FORWARD:
      /* Not require special I/O initialization */
      break;
    case GPU_REF_GEM_FILE:
      GPU_ERROR(gpu_io_load_reference_specs_GEM_FULL(referenceRaw, ref));
      break;
    case GPU_REF_MFASTA_FILE:
      GPU_ERROR(gpu_io_load_reference_specs_MFASTA(referenceRaw, ref));
      break;
    case GPU_REF_PROFILE_FILE:
      GPU_ERROR(gpu_io_load_reference_specs_PROFILE(referenceRaw, ref));
      break;
    default:
      GPU_ERROR(E_REFERENCE_CODING);
      break;
    }

  return(SUCCESS);
}

gpu_error_t gpu_reference_transform(gpu_reference_buffer_t* const ref, const char* const referenceRaw,
                                    const gpu_ref_coding_t refCoding, const gpu_module_t activeModules)
{
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  switch(refCoding){
    case GPU_REF_ASCII:
      GPU_ERROR(gpu_reference_transform_ASCII(referenceRaw, ref));
      break;
    case GPU_REF_GEM_ONLY_FORWARD:
      GPU_ERROR(gpu_reference_transform_GEM((gpu_gem_ref_dto_t*)referenceRaw, ref));
      break;
    case GPU_REF_GEM_FULL:
      GPU_ERROR(gpu_reference_transform_GEM_FULL((gpu_gem_ref_dto_t*)referenceRaw, ref));
      break;
    case GPU_REF_GEM_FILE:
      GPU_ERROR(gpu_io_load_reference_GEM_FULL(referenceRaw, ref));
      break;
    case GPU_REF_MFASTA_FILE:
      GPU_ERROR(gpu_io_load_reference_MFASTA(referenceRaw, ref));
      break;
    case GPU_REF_PROFILE_FILE:
      GPU_ERROR(gpu_io_load_reference_PROFILE(referenceRaw, ref));
      break;
    default:
      GPU_ERROR(E_REFERENCE_CODING);
      break;
  }

  return(SUCCESS);
}

gpu_error_t gpu_reference_init(gpu_reference_buffer_t **reference, const gpu_reference_dto_t* const referenceRaw,
                               const uint32_t numSupportedDevices, const gpu_module_t activeModules)
{
  gpu_reference_buffer_t* const ref = (gpu_reference_buffer_t *) malloc(sizeof(gpu_reference_buffer_t));
  uint32_t idSupDevice;

  GPU_ERROR(gpu_reference_init_dto(ref));

  ref->activeModules  = activeModules;
  ref->size           = referenceRaw->refSize;
  ref->numEntries     = GPU_DIV_CEIL(ref->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

  ref->d_reference = (uint64_t **) malloc(numSupportedDevices * sizeof(uint64_t *));
  if (ref->d_reference == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  ref->memorySpace = (memory_alloc_t *) malloc(numSupportedDevices * sizeof(memory_alloc_t));
  if (ref->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);

  for(idSupDevice = 0; idSupDevice < numSupportedDevices; ++idSupDevice){
    ref->d_reference[idSupDevice] = NULL;
    ref->memorySpace[idSupDevice] = GPU_NONE_MAPPED;
  }

  GPU_ERROR(gpu_reference_set_specs(ref, referenceRaw->reference, referenceRaw->refCoding, activeModules));

  (* reference) = ref;
  return (SUCCESS);
}

gpu_error_t gpu_reference_allocate(gpu_reference_buffer_t *reference, const gpu_module_t activeModules)
{
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  reference->numEntries = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  if(reference->hostAllocStats & GPU_PAGE_LOCKED){
    CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));
  }else{
    reference->h_reference = malloc(reference->numEntries * sizeof(uint64_t));
    if (reference->h_reference == NULL) return (E_DATA_NOT_ALLOCATED);
  }

  return(SUCCESS);
}

gpu_error_t gpu_reference_load(gpu_reference_buffer_t *reference, const gpu_reference_dto_t* const referenceRaw,
                               const gpu_module_t activeModules)
{
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  GPU_ERROR(gpu_reference_allocate(reference, activeModules));
  GPU_ERROR(gpu_reference_transform(reference, referenceRaw->reference, referenceRaw->refCoding, activeModules));

  return(SUCCESS);
}

gpu_error_t gpu_reference_transfer_CPU_to_GPUs(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices,
                                               const gpu_module_t activeModules)
{
  uint32_t deviceFreeMemory, idSupportedDevice;
  uint32_t numSupportedDevices = devices[0]->numSupportedDevices;

  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(reference->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
      const size_t cpySize = reference->numEntries * sizeof(uint64_t);
      deviceFreeMemory = gpu_device_get_free_memory(devices[idSupportedDevice]->idDevice);
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

gpu_error_t gpu_reference_free_host(gpu_reference_buffer_t* const reference)
{
    if(reference->h_reference != NULL){
      if(reference->hostAllocStats == GPU_PAGE_LOCKED) CUDA_ERROR(cudaFreeHost(reference->h_reference));
      else free(reference->h_reference);
      reference->h_reference = NULL;
    }

    return(SUCCESS);
}

gpu_error_t gpu_reference_free_unused_host(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices,
                                           const gpu_module_t activeModules)
{
  uint32_t idSupportedDevice, numSupportedDevices;
  bool referenceInHostSideUsed = false;

  if(activeModules & GPU_REFERENCE){
    numSupportedDevices = devices[0]->numSupportedDevices;
    //Free all the unused references in the host side
    for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
      if(reference->memorySpace[idSupportedDevice] == GPU_HOST_MAPPED) referenceInHostSideUsed = true;
    }

    if(!referenceInHostSideUsed){
      GPU_ERROR(gpu_reference_free_host(reference));
    }
  }

  return(SUCCESS);
}

gpu_error_t gpu_reference_free_device(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices)
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

gpu_error_t gpu_reference_free(gpu_reference_buffer_t **reference, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  gpu_reference_buffer_t* ref = (* reference);

  if(activeModules & GPU_REFERENCE){
    GPU_ERROR(gpu_reference_free_host(ref));
    GPU_ERROR(gpu_reference_free_device(ref, devices));
  }

  if(ref != NULL){
    free(ref->memorySpace);
    free(ref);
    ref = NULL;
  }

  (* reference) = ref;
  return(SUCCESS);
}

#endif /* GPU_REFERENCE_C_ */

