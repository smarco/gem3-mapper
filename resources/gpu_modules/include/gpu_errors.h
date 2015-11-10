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
	E_NOT_IMPLEMENTED
} gpu_error_t;

#define	CUDA_ERROR(error)		(cudaError(error, __FILE__, __LINE__ ))
#define	GPU_ERROR(error)		(gpuError(error, __FILE__, __LINE__ ))


/************************************************************
Functions to handle errors
************************************************************/

GPU_INLINE void cudaError(cudaError_t err, const char *file,  int line )
{
   	if (err != cudaSuccess) {
      		fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err),  file, line );
       		exit(EXIT_FAILURE);
   	}
}

GPU_INLINE const char* gpuGetErrorString(gpu_error_t error)
{
    switch(error) {
        case E_OPENING_FILE:  			return "GEM GPU - Error: opening file";
        case E_READING_FILE:  			return "GEM GPU - Error: reading file";
        case E_INSUFFICIENT_MEM_GPU:	return "GEM GPU - Error: there aren't enough GPU memory space";
        case E_ALLOCATE_MEM: 			return "GEM GPU - Error: allocating data";
        case E_INCOMPATIBLE_GPU:		return "GEM GPU - Error: incompatible GPU (old CC version)";
        case E_REFERENCE_CODING:		return "GEM GPU - Error: reference coding not supported";
        default: 						return "GEM GPU - Unknown error";
    }
}

GPU_INLINE void gpuError(gpu_error_t err, const char *file,  int line)
{
   	if (err != SUCCESS) {
      		fprintf(stderr, "%s in %s at line %d\n", gpuGetErrorString(err),  file, line );
       		exit(EXIT_FAILURE);
   	}
}

#endif /* GPU_ERRORS_H_ */
