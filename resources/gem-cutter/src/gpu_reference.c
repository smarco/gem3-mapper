/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
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

size_t gpu_reference_get_size(gpu_reference_buffer_t* const reference, size_t* const referenceSize, const gpu_module_t activeModules)
{
  size_t bytesPerReference = 0;
  const size_t bytesPerReferencePlain  = reference->numEntriesPlain  * GPU_REFERENCE_PLAIN__ENTRY_SIZE;
  const size_t bytesPerReferenceMasked = reference->numEntriesMasked * GPU_REFERENCE_MASKED__ENTRY_SIZE;
  if(activeModules & GPU_REFERENCE_PLAIN)  bytesPerReference += bytesPerReferencePlain;
  if(activeModules & GPU_REFERENCE_MASKED) bytesPerReference += bytesPerReferenceMasked;
  (* referenceSize) = bytesPerReference;
  return(SUCCESS);
}


/************************************************************
String basic functions
************************************************************/

uint64_t gpu_char_to_bin_ASCII(const unsigned char base)
{
  // Remove the higher bits of the base representation
  const uint64_t shiftReverseBases = GPU_REFERENCE_PLAIN__ENTRY_LENGTH - GPU_REFERENCE_PLAIN__CHAR_LENGTH;
  switch(base)
  {
    case 'A':
    case 'a':
      return(GPU_ENC_DNA_CHAR_A << shiftReverseBases);
    case 'C':
    case 'c':
      return(GPU_ENC_DNA_CHAR_C << shiftReverseBases);
    case 'G':
    case 'g':
      return(GPU_ENC_DNA_CHAR_G << shiftReverseBases);
    case 'T':
    case 't':
      return(GPU_ENC_DNA_CHAR_T << shiftReverseBases);
    default :
      return(GPU_ENC_DNA_CHAR_A << shiftReverseBases);
  }
}

uint64_t gpu_char_is_not_base_ASCII(const unsigned char base)
{
  switch(base)
  {
    case 'A':
    case 'a':
    case 'C':
    case 'c':
    case 'G':
    case 'g':
    case 'T':
    case 't':
      return(GPU_UINT64_ZEROS);
    default :
      return(GPU_UINT64_MASK_ONE_LOW);
  }
}

uint64_t gpu_char_is_not_base(const unsigned char base)
{
  if((base == GPU_ENC_DNA_CHAR_A) || (base == GPU_ENC_DNA_CHAR_C) ||
     (base == GPU_ENC_DNA_CHAR_G) || (base == GPU_ENC_DNA_CHAR_T))
    return(GPU_UINT64_ZEROS);
  else
    return(GPU_UINT64_MASK_ONE_LOW);
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

gpu_error_t gpu_reference_transform_ASCII(const char* const referenceASCII, gpu_reference_buffer_t* const reference,
										  const gpu_module_t activeModules)
{
  uint64_t indexBase, bitmap;
  uint64_t idEntry, i, referencePosition;
  unsigned char referenceChar;

  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  reference->numEntriesPlain = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  for(idEntry = 0; idEntry < reference->numEntriesPlain; ++idEntry){
    bitmap = 0;
    for(i = 0; i < GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY; ++i){
      referencePosition = idEntry * GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY + i;
      if (referencePosition < reference->size){
        referenceChar = referenceASCII[referencePosition];
      }else{
        referenceChar = 'A'; //filling reference padding
      }
      indexBase = gpu_char_to_bin_ASCII(referenceChar);
      bitmap = (bitmap >> GPU_REFERENCE_PLAIN__CHAR_LENGTH) | indexBase;
    }
    reference->h_reference_plain[referencePosition / GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY] = bitmap;
  }
  return(SUCCESS);
}

gpu_error_t gpu_reference_transform_GEM(const gpu_gem_ref_dto_t* const gem_reference, gpu_reference_buffer_t* const reference,
										const gpu_module_t activeModules)
{
  uint64_t bitmap, base;
  uint64_t idEntry, i, referencePosition;

  // Compute size of the forward reference
  const char* const h_gem_reference = gem_reference->reference;
  const uint64_t total_ref_size     = gem_reference->ref_length;
  const uint64_t baseMask           = GPU_UINT64_ONES << GPU_REFERENCE_CHAR_LENGTH;

  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  reference->size = total_ref_size;
  reference->numEntriesPlain = GPU_DIV_CEIL(total_ref_size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;

  for(idEntry = 0; idEntry < reference->numEntriesPlain; ++idEntry){
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
    reference->h_reference_plain[referencePosition / GPU_REFERENCE_CHARS_PER_ENTRY] = bitmap;
  }
  return(SUCCESS);
}


gpu_error_t gpu_reference_transform_masked_GEM_FULL(const char* const h_gem_reference, uint64_t* const h_reference,
												   const uint64_t refForwardSize, const uint64_t refCompleteSize, const uint64_t numEntries)
{
  // Process the reference compacting the alphabet (4 bases in 2 bits)
  for(uint64_t idEntry = 0; idEntry < numEntries; ++idEntry){
    uint64_t bitmap = 0, referencePosition;
    for(uint64_t idBase = 0; idBase < GPU_REFERENCE_MASKED__CHARS_PER_ENTRY; ++idBase){
      // By default fills reference padding with non-base
      uint64_t base = GPU_UINT64_ZEROS;
      referencePosition = idEntry * GPU_REFERENCE_MASKED__CHARS_PER_ENTRY + idBase;
      if (referencePosition < refForwardSize)
          base = gpu_char_is_not_base(h_gem_reference[referencePosition]);
      else if (referencePosition < refCompleteSize)
          base = (uint64_t) gpu_char_is_not_base(h_gem_reference[(2 * refForwardSize) - referencePosition - 2]);
      // Packing the current base in a reference entry representation
      base = base << (GPU_REFERENCE_MASKED__ENTRY_LENGTH - GPU_REFERENCE_MASKED__CHAR_LENGTH);
      bitmap = (bitmap >> GPU_REFERENCE_MASKED__CHAR_LENGTH) | base;
    }
    // Saving entry in the packed reference representation
    h_reference[referencePosition / GPU_REFERENCE_MASKED__CHARS_PER_ENTRY] = bitmap;
  }
  // Return
  return(SUCCESS);
}

gpu_error_t gpu_reference_transform_plain_GEM_FULL(const char* const h_gem_reference, uint64_t* const h_reference,
		                                           const uint64_t refForwardSize, const uint64_t refCompleteSize, const uint64_t numEntries)
{
  // Recompute size of the full reference (forward + reverse-complement)
  const uint64_t baseMask = GPU_UINT64_ONES << GPU_REFERENCE_PLAIN__CHAR_LENGTH;
  // Process the reference compacting the alphabet (4 bases in 2 bits)
  for(uint64_t idEntry = 0; idEntry < numEntries; ++idEntry){
    uint64_t bitmap = 0, referencePosition;
    for(uint64_t idBase = 0; idBase < GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY; ++idBase){
      // By default fills reference padding
      uint64_t base = GPU_ENC_DNA_CHAR_A;
      referencePosition = idEntry * GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY + idBase;
      if (referencePosition < refForwardSize)
        base = (uint64_t) h_gem_reference[referencePosition];
      else if (referencePosition < refCompleteSize)
        base = (uint64_t) gpu_complement_base(h_gem_reference[(2 * refForwardSize) - referencePosition - 2]);
      // Setting non-bases (different to A,C,G,T) to As
      base = (base & baseMask) ? GPU_ENC_DNA_CHAR_A : base;
      // Packing the current base in a reference entry representation
      base = base << (GPU_REFERENCE_PLAIN__ENTRY_LENGTH - GPU_REFERENCE_PLAIN__CHAR_LENGTH);
      bitmap = (bitmap >> GPU_REFERENCE_PLAIN__CHAR_LENGTH) | base;
    }
    // Saving entry in the packed reference representation
    h_reference[referencePosition / GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY] = bitmap;
  }
  // Return
  return(SUCCESS);
}

gpu_error_t gpu_reference_transform_GEM_FULL(const gpu_gem_ref_dto_t* const gem_reference, gpu_reference_buffer_t* const reference,
											 const gpu_module_t activeModules)
{
  // Recompute size of the full reference (forward + reverse-complement)
  const char* const h_gem_reference = gem_reference->reference;
  const uint64_t forward_ref_size   = gem_reference->ref_length;
  const uint64_t total_ref_size     = 2 * forward_ref_size;
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);
  //Transform the reference content in a compact representation
  reference->size = total_ref_size;
  reference->numEntriesPlain  = GPU_DIV_CEIL(total_ref_size, GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  reference->numEntriesMasked = GPU_DIV_CEIL(total_ref_size, GPU_REFERENCE_MASKED__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  //Transform the reference content in a compact representation
  GPU_ERROR(gpu_reference_transform_plain_GEM_FULL(h_gem_reference, reference->h_reference_plain, forward_ref_size, total_ref_size, reference->numEntriesPlain));
  GPU_ERROR(gpu_reference_transform_masked_GEM_FULL(h_gem_reference, reference->h_reference_masked, forward_ref_size, total_ref_size, reference->numEntriesMasked));
  // Return
  return(SUCCESS);
}


/************************************************************
Input & Output reference functions
************************************************************/

gpu_error_t gpu_reference_read_specs(int fp, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules)
{
  size_t result, bytesRequest;
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);
  // Read the reference specifications
  bytesRequest = sizeof(uint64_t);
  result = read(fp, (void* )&reference->numEntriesPlain, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  bytesRequest = sizeof(uint64_t);
  result = read(fp, (void* )&reference->numEntriesMasked, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  bytesRequest = sizeof(uint64_t);
  result = read(fp, (void* )&reference->size, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_reference_write_specs(int fp, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules)
{
  size_t result, bytesRequest;
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);
  // Write the reference specifications
  bytesRequest = sizeof(uint64_t);
  result = write(fp, (void* )&reference->numEntriesPlain, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  bytesRequest = sizeof(uint64_t);
  result = write(fp, (void* )&reference->numEntriesMasked, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  bytesRequest = sizeof(uint64_t);
  result = write(fp, (void* )&reference->size, bytesRequest);
  if (result != bytesRequest) return (E_READING_FILE);
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_reference_read(int fp, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules)
{
  // Request size definition
  size_t bytesRequest = 0;
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);
  // Read the reference specifications
  GPU_ERROR(gpu_reference_read_specs(fp, reference, activeModules));
  // Read the plain reference in a block iterative way
  bytesRequest = GPU_REFERENCE_PLAIN__ENTRY_SIZE * reference->numEntriesPlain;
  GPU_ERROR(gpu_io_read_buffered(fp, (void* )reference->h_reference_plain,  bytesRequest));
  // Read the masked reference in a block iterative way
  bytesRequest = GPU_REFERENCE_MASKED__ENTRY_SIZE * reference->numEntriesMasked;
  GPU_ERROR(gpu_io_read_buffered(fp, (void* )reference->h_reference_masked, bytesRequest));
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_reference_write(int fp, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules)
{
  // Request size definition
  size_t bytesRequest = 0;
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);
  // Write the reference specifications
  GPU_ERROR(gpu_reference_write_specs(fp, reference, activeModules));
  // Write the plain reference in a block iterative way
  bytesRequest = GPU_REFERENCE_PLAIN__ENTRY_SIZE * reference->numEntriesPlain;
  GPU_ERROR(gpu_io_write_buffered(fp, (void* )reference->h_reference_plain,  bytesRequest));
  // Write the masked reference in a block iterative way
  bytesRequest = GPU_REFERENCE_MASKED__ENTRY_SIZE * reference->numEntriesMasked;
  GPU_ERROR(gpu_io_write_buffered(fp, (void* )reference->h_reference_masked, bytesRequest));
  // Succeed
  return (SUCCESS);
}


/************************************************************
Initialize reference functions
************************************************************/

gpu_error_t gpu_reference_init_dto(gpu_reference_buffer_t* const ref)
{
  //Initialize the reference structure
  ref->size                = 0;
  //Initialize basic reference
  ref->numEntriesPlain     = 0;
  ref->h_reference_plain   = NULL;
  ref->d_reference_plain   = NULL;
  //Initialize masked reference
  ref->numEntriesMasked    = 0;
  ref->h_reference_masked  = NULL;
  ref->d_reference_masked  = NULL;
  //Initialize reference allocation policies
  ref->hostAllocStats 	   = GPU_PAGE_UNLOCKED;
  ref->memorySpace         = NULL;
  ref->activeModules       = GPU_NONE_MODULES;
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_reference_set_specs(gpu_reference_buffer_t* const ref, const char* const referenceRaw,
                                    const gpu_ref_coding_t refCoding, const gpu_module_t activeModules)
{
  // Module sanity checker
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
      GPU_ERROR(gpu_io_load_reference_specs_GEM_FULL(referenceRaw, ref, activeModules));
      break;
    case GPU_REF_MFASTA_FILE:
      GPU_ERROR(gpu_io_load_reference_specs_MFASTA(referenceRaw, ref, activeModules));
      break;
    case GPU_REF_PROFILE_FILE:
      GPU_ERROR(gpu_io_load_reference_specs_PROFILE(referenceRaw, ref, activeModules));
      break;
    default:
      GPU_ERROR(E_REFERENCE_CODING);
      break;
    }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_reference_transform(gpu_reference_buffer_t* const ref, const char* const referenceRaw,
                                    const gpu_ref_coding_t refCoding, const gpu_module_t activeModules)
{
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);

  switch(refCoding){
    case GPU_REF_ASCII:
      GPU_ERROR(gpu_reference_transform_ASCII(referenceRaw, ref, activeModules));
      break;
    case GPU_REF_GEM_ONLY_FORWARD:
      GPU_ERROR(gpu_reference_transform_GEM((gpu_gem_ref_dto_t*)referenceRaw, ref, activeModules));
      break;
    case GPU_REF_GEM_FULL:
      GPU_ERROR(gpu_reference_transform_GEM_FULL((gpu_gem_ref_dto_t*)referenceRaw, ref, activeModules));
      break;
    case GPU_REF_GEM_FILE:
      GPU_ERROR(gpu_io_load_reference_GEM_FULL(referenceRaw, ref, activeModules));
      break;
    case GPU_REF_MFASTA_FILE:
      GPU_ERROR(gpu_io_load_reference_MFASTA(referenceRaw, ref, activeModules));
      break;
    case GPU_REF_PROFILE_FILE:
      GPU_ERROR(gpu_io_load_reference_PROFILE(referenceRaw, ref, activeModules));
      break;
    default:
      GPU_ERROR(E_REFERENCE_CODING);
      break;
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_reference_init(gpu_reference_buffer_t **reference, const gpu_reference_dto_t* const referenceRaw,
                               const uint32_t numSupportedDevices, const gpu_module_t activeModules)
{
  gpu_reference_buffer_t* const ref = (gpu_reference_buffer_t *) malloc(sizeof(gpu_reference_buffer_t));
  uint32_t idSupDevice;
  // Creating reference structure
  GPU_ERROR(gpu_reference_init_dto(ref));
  // Initializing reference attributes
  ref->activeModules    = activeModules & GPU_REFERENCE;
  ref->size             = referenceRaw->refSize;
  ref->numEntriesPlain  = GPU_DIV_CEIL(ref->size, GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY)  + GPU_REFERENCE_END_PADDING;
  ref->numEntriesMasked = GPU_DIV_CEIL(ref->size, GPU_REFERENCE_MASKED__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  // Allocating the description for the GPU reference collections
  ref->d_reference_plain = (uint64_t **) malloc(numSupportedDevices * sizeof(uint64_t *));
  if (ref->d_reference_plain == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  ref->d_reference_masked = (uint64_t **) malloc(numSupportedDevices * sizeof(uint64_t *));
  if (ref->d_reference_masked == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  ref->memorySpace = (memory_alloc_t *) malloc(numSupportedDevices * sizeof(memory_alloc_t));
  if (ref->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  // Devices initializations
  for(idSupDevice = 0; idSupDevice < numSupportedDevices; ++idSupDevice){
    ref->d_reference_plain[idSupDevice]  = NULL;
    ref->d_reference_masked[idSupDevice] = NULL;
    ref->memorySpace[idSupDevice] = GPU_NONE_MAPPED;
  }
  // Setting file headers and descriptions
  GPU_ERROR(gpu_reference_set_specs(ref, referenceRaw->reference, referenceRaw->refCoding, activeModules));
  (* reference) = ref;
  // Succeed
  return (SUCCESS);
}

gpu_error_t gpu_reference_allocate(gpu_reference_buffer_t *reference, const gpu_module_t activeModules)
{
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);
  // Setting size for host allocations dedicated to device transference
  const uint64_t numEntriesPlain  = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_PLAIN__CHARS_PER_ENTRY)  + GPU_REFERENCE_END_PADDING;
  const uint64_t numEntriesMasked = GPU_DIV_CEIL(reference->size, GPU_REFERENCE_MASKED__CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  const size_t   cpySizeRefPlain  = numEntriesPlain  * GPU_REFERENCE_PLAIN__ENTRY_SIZE;
  const size_t   cpySizeRefMasked = numEntriesMasked * GPU_REFERENCE_MASKED__ENTRY_SIZE;
  // Setting reference sizes
  reference->numEntriesPlain  = numEntriesPlain;
  reference->numEntriesMasked = numEntriesMasked;
  // Pinned or non-pinned allocation to optimize zero transfer requests
  if(reference->hostAllocStats & GPU_PAGE_LOCKED){
   CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference_plain, cpySizeRefPlain, cudaHostAllocMapped));
   CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference_masked, cpySizeRefMasked, cudaHostAllocMapped));
  }else{
    reference->h_reference_plain = malloc(cpySizeRefPlain);
    if (reference->h_reference_plain == NULL) return (E_DATA_NOT_ALLOCATED);
    reference->h_reference_masked = malloc(cpySizeRefMasked);
    if (reference->h_reference_masked == NULL) return (E_DATA_NOT_ALLOCATED);
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_reference_load(gpu_reference_buffer_t *reference, const gpu_reference_dto_t* const referenceRaw,
                               const gpu_module_t activeModules)
{
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);
  // Managing the device references
  GPU_ERROR(gpu_reference_allocate(reference, activeModules));
  GPU_ERROR(gpu_reference_transform(reference, referenceRaw->reference, referenceRaw->refCoding, activeModules));
  // Succeed
  return(SUCCESS);
}

//
gpu_error_t gpu_reference_transfer_CPU_to_GPUs(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices,
                                               const gpu_module_t activeModules)
{
  uint32_t deviceFreeMemory = 0, idSupportedDevice = 0;
  uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
  // Module sanity checker
  if((activeModules & GPU_REFERENCE) == 0)
    return(E_MODULE_NOT_FOUND);
  // Transfer data for each GPU on the system
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    if(reference->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
      size_t cpySize = 0;
      gpu_reference_get_size(reference, &cpySize, activeModules);
      deviceFreeMemory = gpu_device_get_free_memory(devices[idSupportedDevice]->idDevice);
      if ((GPU_CONVERT__B_TO_MB(cpySize)) > deviceFreeMemory) return(E_INSUFFICIENT_MEM_GPU);
      CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
      //Synchronous allocate & transfer plain reference to GPU
      if(activeModules & GPU_REFERENCE_PLAIN){
        gpu_reference_get_size(reference, &cpySize, GPU_REFERENCE_PLAIN);
    	CUDA_ERROR(cudaMalloc((void**) &reference->d_reference_plain[idSupportedDevice], cpySize));
        CUDA_ERROR(cudaMemcpy(reference->d_reference_plain[idSupportedDevice], reference->h_reference_plain, cpySize, cudaMemcpyHostToDevice));
      }
      //Synchronous allocate & transfer masked reference to GPU
      if(activeModules & GPU_REFERENCE_MASKED){
        gpu_reference_get_size(reference, &cpySize, GPU_REFERENCE_MASKED);
    	CUDA_ERROR(cudaMalloc((void**) &reference->d_reference_masked[idSupportedDevice], cpySize));
        CUDA_ERROR(cudaMemcpy(reference->d_reference_masked[idSupportedDevice], reference->h_reference_masked, cpySize, cudaMemcpyHostToDevice));
      }
    }else{
      if(activeModules & GPU_REFERENCE_PLAIN)
    	  reference->d_reference_plain[idSupportedDevice] = reference->h_reference_plain;
      if(activeModules & GPU_REFERENCE_MASKED)
    	  reference->d_reference_masked[idSupportedDevice] = reference->h_reference_masked;
    }
  }
  // Succeed
  return (SUCCESS);
}


/************************************************************
Free reference functions
************************************************************/

gpu_error_t gpu_reference_free_host(gpu_reference_buffer_t* const reference)
{
  // Deallocate plain reference from host
  if(reference->h_reference_plain != NULL){
    if(reference->hostAllocStats == GPU_PAGE_LOCKED) CUDA_ERROR(cudaFreeHost(reference->h_reference_plain));
    else free(reference->h_reference_plain);
    reference->h_reference_plain = NULL;
  }
  // Deallocate masked reference from host
  if(reference->h_reference_masked != NULL){
    if(reference->hostAllocStats == GPU_PAGE_LOCKED) CUDA_ERROR(cudaFreeHost(reference->h_reference_masked));
    else free(reference->h_reference_masked);
    reference->h_reference_masked = NULL;
  }
  // Succeed
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
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_reference_free_device(gpu_reference_buffer_t* const reference, gpu_device_info_t** const devices)
{
  const uint32_t numSupportedDevices = devices[0]->numSupportedDevices;
  uint32_t idSupportedDevice;
  // Free all the references in the devices
  for(idSupportedDevice = 0; idSupportedDevice < numSupportedDevices; ++idSupportedDevice){
    // Selecting GPU device
    CUDA_ERROR(cudaSetDevice(devices[idSupportedDevice]->idDevice));
    if(reference->memorySpace[idSupportedDevice] == GPU_DEVICE_MAPPED){
      // Free the device plain reference
      if(reference->d_reference_plain[idSupportedDevice] != NULL)
        CUDA_ERROR(cudaFree(reference->d_reference_plain[idSupportedDevice]));
      // Free the device masked reference
      if(reference->d_reference_masked[idSupportedDevice] != NULL)
        CUDA_ERROR(cudaFree(reference->d_reference_masked[idSupportedDevice]));
      // Resetting pointer values
      reference->d_reference_plain[idSupportedDevice]  = NULL;
      reference->d_reference_masked[idSupportedDevice] = NULL;
    }
  }
  // Free the device plain reference list
  if(reference->d_reference_plain != NULL){
    free(reference->d_reference_plain);
    reference->d_reference_plain = NULL;
  }
  // Free the device masked reference list
  if(reference->d_reference_masked != NULL){
    free(reference->d_reference_masked);
    reference->d_reference_masked = NULL;
  }
  // Succeed
  return(SUCCESS);
}

gpu_error_t gpu_reference_free(gpu_reference_buffer_t **reference, gpu_device_info_t** const devices, const gpu_module_t activeModules)
{
  gpu_reference_buffer_t* ref = (* reference);
  // Free device and host references
  if(activeModules & GPU_REFERENCE){
    GPU_ERROR(gpu_reference_free_host(ref));
    GPU_ERROR(gpu_reference_free_device(ref, devices));
  }
  // Free memory space specifications
  if(ref != NULL){
    free(ref->memorySpace);
    free(ref);
    ref = NULL;
  }
  // Return reference changes
  (* reference) = ref;
  // Succeed
  return(SUCCESS);
}

#endif /* GPU_REFERENCE_C_ */

