/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
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
  E_NOT_SUPPORTED_GPUS,
  E_REFERENCE_CODING,
  E_INDEX_CODING,
  E_WRITING_FILE,
  E_MODULE_NOT_FOUND,
  E_INSUFFICIENT_MEM_PER_BUFFER,
  E_DATA_NOT_ALLOCATED,
  E_NOT_SUPPORTED_ALLOC_POLICY,
  E_OVERFLOWING_BUFFER,
  E_FMI_TABLE_INCOMPATIBLE_SIZE,
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
