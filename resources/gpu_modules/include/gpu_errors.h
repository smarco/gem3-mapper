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
	E_NOT_IMPLEMENTED
} gpu_error_t;

#define	CUDA_ERROR(error)		(cudaError(error, __FILE__, __LINE__ ))
#define	GPU_ERROR(error)		(gpuError(error, __FILE__, __LINE__ ))

void gpuError(gpu_error_t error, const char *file, int32_t line);
char *gpuGetErrorString(gpu_error_t error);

#endif /* GPU_ERRORS_H_ */
