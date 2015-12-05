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
