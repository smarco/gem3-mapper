/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include "gpu_commons.h"
#include "gpu_module.h"
/* Include the required objects */
#include "gpu_devices.h"
#include "gpu_reference.h"
#include "gpu_index.h"
/* Include the required modules */
#include "gpu_buffer_modules.h"

#ifndef GPU_BUFFER_H_
#define GPU_BUFFER_H_

typedef struct {
  gpu_module_t            typeBuffer;
  uint32_t                numBuffers;
  uint32_t                idBuffer;
  uint32_t                idStream;
  cudaStream_t            *listStreams;
  uint32_t                idSupportedDevice;
  gpu_device_info_t       **device;
  gpu_reference_buffer_t  *reference;
  gpu_index_buffer_t      *index;
  size_t                  sizeBuffer;
  void                    *h_rawData;
  void                    *d_rawData;
  gpu_buffer_modules_t    data;
} gpu_buffer_t;


/* Primitives to get information from buffers */
gpu_error_t gpu_buffer_get_min_memory_size(size_t *bytesPerBuffer);

/* Primitives to schedule and manage the buffers */
gpu_error_t gpu_buffer_configuration(gpu_buffer_t* const mBuff, const uint32_t idBuffer, const uint32_t idSupportedDevice,
                                     const size_t bytesPerBuffer, const uint32_t numBuffers, gpu_device_info_t** const device,
                                     cudaStream_t* const listStreams, gpu_reference_buffer_t* const reference, gpu_index_buffer_t* const index);
gpu_error_t gpu_buffer_scheduling(gpu_buffer_t ***gpuBuffer, const uint32_t numBuffers, gpu_device_info_t** const device,
                                  gpu_reference_buffer_t *reference, gpu_index_buffer_t *index, float maxMbPerBuffer);

/* Functions to free all the buffer resources (HOST & DEVICE) */
gpu_error_t gpu_buffer_free(gpu_buffer_t *mBuff);


#endif /* GPU_BUFFER_H_ */
