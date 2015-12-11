#include "../include/gpu_io.h"


/************************************************************
Primitives for input/output
************************************************************/

// Input-file
void gpu_load_indexed_structures_()
{
  FILE *fp = NULL;
  fp = fopen(fn, "rb");
  if (fp == NULL) return (E_OPENING_FILE);
  //TODO
  //Check active modules
  //Conditional load
  //Pinned allocation
  gpu_load_reference_PROFILE(fp, gpu_reference_buffer_t *reference);
}


//gpu_error_t gpu_load_reference(const char *fn, gpu_reference_buffer_t *reference)
//{
//  FILE *fp = NULL;
//  size_t result;
//
//  fp = fopen(fn, "rb");
//  if (fp == NULL) return (E_OPENING_FILE);
//
//  result = fread(&reference->numEntries, sizeof(uint64_t), 1, fp);
//  if (result != 1) return (E_READING_FILE);
//  result = fread(&reference->size, sizeof(uint64_t), 1, fp);
//  if (result != 1) return (E_READING_FILE);
//
//  CUDA_ERROR(cudaHostAlloc((void**) &reference->h_reference, reference->numEntries * sizeof(uint64_t), cudaHostAllocMapped));
//
//  result = fread(reference->h_reference, sizeof(uint64_t), reference->numEntries, fp);
//  if (result != reference->numEntries) return (E_READING_FILE);
//
//  fclose(fp);
//  return (SUCCESS);
//}

void gpu_load_indexed_structures(const char* const fileName, gpu_index_buffer_t *rawIndex, gpu_reference_buffer_t* rawRef)
{
}

void gpu_store_indexed_structures_(const char* const fileName, gpu_gem_fmi_dto_t *gemIndex, gpu_gem_ref_dto_t* gemRef)
{
  gpu_reference_buffer_t ref;
  gpu_index_buffer_t     index;

  ref->d_reference   = NULL;
  ref->h_reference   = NULL;
  ref->memorySpace   = NULL;
  ref->size          = 0;
  ref->numEntries    = GPU_DIV_CEIL(ref->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
  ref->activeModules = activeModules;

  ref->d_reference = (uint64_t **) malloc(numSupportedDevices * sizeof(uint64_t *));
  if (ref->d_reference == NULL) GPU_ERROR(E_ALLOCATE_MEM);
  ref->memorySpace = (memory_alloc_t *) malloc(numSupportedDevices * sizeof(memory_alloc_t));
  if (ref->memorySpace == NULL) GPU_ERROR(E_ALLOCATE_MEM);

  GPU_ERROR(gpu_init_reference(&reference, referenceRaw, refSize, refCoding, numSupportedDevices, activeModules));
  GPU_ERROR(gpu_init_index(&index, fmiRaw, bwtSize, indexCoding, numSupportedDevices, activeModules));

  if(activeModules & GPU_BPM){
    ref->size = refSize;
    ref->numEntries = GPU_DIV_CEIL(ref->size, GPU_REFERENCE_CHARS_PER_ENTRY) + GPU_REFERENCE_END_PADDING;
    GPU_ERROR(gpu_transform_reference_GEM_FR(gemRef, ref));
  }

  if(activeModules & (GPU_FMI_DECODE_POS | GPU_FMI_EXACT_SEARCH)){
    fmi->bwtSize       = bwtSize;
    fmi->numEntries    = GPU_DIV_CEIL(fmi->bwtSize, GPU_FMI_ENTRY_SIZE) + 1;
    char *h_BWT = NULL;
    GPU_ERROR(gpu_transform_index_GEM_FULL(gemIndex, fmi));
  }
}
