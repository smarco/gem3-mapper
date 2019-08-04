/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_SCHEDULER_CORE_H_
#define GPU_SCHEDULER_CORE_H_

#define GPU_SCHEDULER_DISABLED_TASK      GPU_UINT32_ONES
#define GPU_SCHEDULER_NONASSIGNED_TASK	(GPU_UINT32_ONES - 1)

extern "C" {
#include "gpu_commons.h"
}

#define GPU_MAX_BUCKETS  GPU_WARP_SIZE

typedef struct {
  uint32_t source;
  uint32_t remapped;
} gpu_task_id_t;

GPU_INLINE __device__ void gpu_scheduler_scatter_work(const uint32_t globalThreadIdx, const uint32_t* const d_initWarpPerBucket,
		                                      const uint32_t* const d_initPosPerBucket, const uint32_t* const d_endPosPerBucket,
                                                      const uint32_t* const d_reorderBuffer, const bool updateScheduling,
		                                      uint32_t* const idTaskRes, uint32_t* const intraTaskThreadIdxRes, uint32_t* const threadsPerTaskRes)
{
  //Warp identification of the current cuda thread
  const uint32_t globalWarpIdx = globalThreadIdx / GPU_WARP_SIZE;
  //Thread initializations
  uint32_t bucketIdx = 0, tasksPerWarp;
  uint32_t idTask = GPU_SCHEDULER_NONASSIGNED_TASK, intraTaskThreadIdx, threadsPerTask;
  uint32_t localThreadInTheBucket, localIdTaskInTheWarp, localIdTaskInTheBucket, startIdTaskPerWarp;
  //Scan in which bucket is matched this warp
  while((bucketIdx != (GPU_MAX_BUCKETS + 1)) && (d_initWarpPerBucket[bucketIdx] <= globalWarpIdx)){
    bucketIdx++;
  }
  bucketIdx--;
  //Rescheduling task ID to thread ID
  threadsPerTask            = bucketIdx + 1;
  tasksPerWarp              = GPU_WARP_SIZE / threadsPerTask;
  localThreadInTheBucket    = globalThreadIdx - (d_initWarpPerBucket[bucketIdx] * GPU_WARP_SIZE);
  // Discards tasks - padded threads in the warp are idle and not assigned to a task
  startIdTaskPerWarp        = (localThreadInTheBucket / GPU_WARP_SIZE) * tasksPerWarp;
  // Setting internal tasks inside a warp
  localIdTaskInTheWarp      = ((threadIdx.x % GPU_WARP_SIZE) / threadsPerTask);
  localIdTaskInTheBucket    = startIdTaskPerWarp + localIdTaskInTheWarp;
  // idTask recalculation
  idTask                    = d_initPosPerBucket[bucketIdx] + localIdTaskInTheBucket;
  // Disabling excess idTasks (tasks between warp buckets)
  idTask		    = (idTask < d_endPosPerBucket[bucketIdx]) ? idTask : GPU_SCHEDULER_DISABLED_TASK;
  intraTaskThreadIdx        = (threadIdx.x % GPU_WARP_SIZE) % threadsPerTask;
  // Update the buffer input/output for the thread re-scheduling
  if(idTask <= d_endPosPerBucket[GPU_MAX_BUCKETS - 1])
    idTask                  = (updateScheduling) ? d_reorderBuffer[idTask] : idTask;
  // Returning new task id and thread work configuration (internal thread ID inside a thread group + Number of threads assigned to a task)
  (* idTaskRes)             = idTask;
  (* intraTaskThreadIdxRes) = intraTaskThreadIdx;
  (* threadsPerTaskRes)     = threadsPerTask;
}

#endif /* GPU_SCHEDULER_CORE_H_ */

