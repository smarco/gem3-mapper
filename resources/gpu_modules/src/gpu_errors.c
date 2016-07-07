/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_ERRORS_C_
#define GPU_ERRORS_C_

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
    case E_NOT_SUPPORTED_GPUS:          return "GEM GPU - Error: there aren't supported GPUs in the system";
    case E_REFERENCE_CODING:            return "GEM GPU - Error: reference coding not supported";
    case E_INDEX_CODING:                return "GEM GPU - Error: index coding not supported";
    case E_INSUFFICIENT_MEM_PER_BUFFER: return "GEM GPU - Error: GPU buffer size too small";
    case E_MODULE_NOT_FOUND:            return "GEM GPU - Error: module not found";
    case E_NOT_SUPPORTED_ALLOC_POLICY:  return "GEM GPU - Error: undefined memory allocation policy";
    case E_NOT_IMPLEMENTED:             return "GEM GPU - Error: functionality not implemented";
    case E_DATA_NOT_ALLOCATED:          return "GEM GPU - Error: structure not allocated to any memory";
    case E_OVERFLOWING_BUFFER:          return "GEM GPU - Error: overflowing elements per buffer";
    case E_FMI_TABLE_INCOMPATIBLE_SIZE: return "GEM GPU - Error: fmi table num levels incompatible";
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

#endif /* GPU_ERRORS_C_ */
