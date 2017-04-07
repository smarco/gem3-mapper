/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_SCHEDULER_CORE_H_
#define GPU_SCHEDULER_CORE_H_

extern "C" {
#include "gpu_commons.h"
}

GPU_INLINE __device__ void gpu_scheduler_scatter_work(const uint32_t globalThreadIdx, const uint32_t* const d_initWarpPerBucket, const uint32_t* const d_initPosPerBucket,
		                                   	   	   	  uint32_t* const idTaskRes, uint32_t* const intraTaskThreadIdxRes, uint32_t* const threadsPerTaskRes)
{
  const uint32_t globalWarpIdx   = globalThreadIdx / GPU_WARP_SIZE;

  uint32_t bucketIdx = 0;
  uint32_t localThreadInTheBucket, idTask, intraTaskThreadIdx, threadsPerTask, tasksPerWarp, localIdTaskInTheBucket;

  //Scan in which bucket is matched this warp
  while((bucketIdx != (GPU_WARP_SIZE + 1)) && (d_initWarpPerBucket[bucketIdx] <= globalWarpIdx)){
    bucketIdx++;
  }
  bucketIdx--;

  threadsPerTask            = bucketIdx + 1;
  tasksPerWarp              = GPU_WARP_SIZE / threadsPerTask;
  localThreadInTheBucket    = globalThreadIdx - (d_initWarpPerBucket[bucketIdx] * GPU_WARP_SIZE);
  localIdTaskInTheBucket    = ((localThreadInTheBucket / GPU_WARP_SIZE) * tasksPerWarp) + ((threadIdx.x % GPU_WARP_SIZE) / threadsPerTask);
  idTask                    = d_initPosPerBucket[bucketIdx] + localIdTaskInTheBucket;
  intraTaskThreadIdx        = (threadIdx.x % GPU_WARP_SIZE) % threadsPerTask;

  (* idTaskRes)             = idTask;
  (* intraTaskThreadIdxRes) = intraTaskThreadIdx;
  (* threadsPerTaskRes)     = threadsPerTask;
}

#endif /* GPU_SCHEDULER_CORE_H_ */

