/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Common headers and data structures for BPM on GPU library
 */

#ifndef GPU_ERRORS_H_
#define GPU_ERRORS_H_

#include "gpu_commons.h"

typedef enum
{
  SUCCESS,
  E_OPENING_FILE,
  E_READING_FILE,
  E_INSUFFICIENT_MEM_GPU,
  E_ALLOCATE_MEM,
  E_INCOMPATIBLE_GPU,
  E_NO_SUPPORTED_GPUS,
  E_REFERENCE_CODING,
  E_INDEX_CODING,
  E_WRITING_FILE,
  E_MODULE_NOT_FOUND,
  E_INSUFFICIENT_MEM_PER_BUFFER,
  E_NOT_IMPLEMENTED
} gpu_error_t;

#define CUDA_ERROR(error)   (cudaError(error, __FILE__, __LINE__ ))
#define GPU_ERROR(error)    (gpuError(error, __FILE__, __LINE__ ))


/************************************************************
Functions to handle errors
************************************************************/

void cudaError(cudaError_t err, const char *file,  int line);
const char* gpuGetErrorString(gpu_error_t error);
void gpuError(gpu_error_t err, const char *file,  int line);

#endif /* GPU_ERRORS_H_ */
