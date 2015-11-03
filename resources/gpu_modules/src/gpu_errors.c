/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Host scheduler for BPM on GPU
 */

#include "../include/gpu_errors.h"

/************************************************************
Functions to handle errors
************************************************************/

GPU_INLINE void cudaError(cudaError_t err, const char *file,  int line ) {
   	if (err != cudaSuccess) {
      		fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err),  file, line );
       		exit(EXIT_FAILURE);
   	}
}

GPU_INLINE char *gpuGetErrorString(gpu_error_t error){
    switch(error) {
        case E_OPENING_FILE:  			return "GEM GPU - Error: opening file"; break;
        case E_READING_FILE:  			return "GEM GPU - Error: reading file"; break;
        case E_INSUFFICIENT_MEM_GPU:	return "GEM GPU - Error: there aren't enough GPU memory space"; break;
        case E_ALLOCATE_MEM: 			return "GEM GPU - Error: allocating data"; break;
        case E_INCOMPATIBLE_GPU:		return "GEM GPU - Error: incompatible GPU (old CC version)"; break;
        case E_REFERENCE_CODING:		return "GEM GPU - Error: reference coding not supported"; break;
        default: 						return "GEM GPU - Unknown error";
    }
}

GPU_INLINE void gpuError(gpu_error_t err, const char *file,  int line ) {
   	if (err != 0) {
      		fprintf(stderr, "%s in %s at line %d\n", gpuGetErrorString(err),  file, line );
       		exit(EXIT_FAILURE);
   	}
}
