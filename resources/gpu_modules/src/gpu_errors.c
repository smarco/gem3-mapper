/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: gpu_errors.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: TODO
 */

#include "../include/gpu_errors.h"

void cudaError(cudaError_t err, const char *file,  int line)
{
  if (err != cudaSuccess) {
    fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err),  file, line);
    exit(EXIT_FAILURE);
  }
}

const char* gpuGetErrorString(gpu_error_t error)
{
  switch(error){
    case E_OPENING_FILE:                return "GEM GPU - Error: opening file";
    case E_READING_FILE:                return "GEM GPU - Error: reading file";
    case E_WRITING_FILE:                return "GEM GPU - Error: writing file";
    case E_INSUFFICIENT_MEM_GPU:        return "GEM GPU - Error: there aren't enough GPU memory space";
    case E_ALLOCATE_MEM:                return "GEM GPU - Error: allocating data";
    case E_INCOMPATIBLE_GPU:            return "GEM GPU - Error: incompatible GPU (old CC version)";
    case E_NO_SUPPORTED_GPUS:           return "GEM GPU - Error: there aren't supported GPUs in the system";
    case E_REFERENCE_CODING:            return "GEM GPU - Error: reference coding not supported";
    case E_INDEX_CODING:                return "GEM GPU - Error: index coding not supported";
    case E_INSUFFICIENT_MEM_PER_BUFFER: return "GEM GPU - Error: GPU buffer size too small";
    case E_MODULE_NOT_FOUND:            return "GEM GPU - Error: module not found";
    case E_NOT_IMPLEMENTED:             return "GEM GPU - Error: functionality not implemented";
    default:                            return "GEM GPU - Unknown error";
  }
}

void gpuError(gpu_error_t err, const char *file,  int line)
{
  if (err != SUCCESS) {
    fprintf(stderr, "%s in %s at line %d\n", gpuGetErrorString(err),  file, line );
    exit(EXIT_FAILURE);
  }
}
