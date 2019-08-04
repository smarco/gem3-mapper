/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */


#ifndef GPU_SCHEDULER_H_
#define GPU_SCHEDULER_H_

/* Defines to define the the task rescheduling order point of view*/
#define GPU_SCHEDULER_NUM_REORDER_BUFFERS                    2
#define GPU_SCHEDULER_TASK_MAPPED                            0
#define GPU_SCHEDULER_THREAD_MAPPED                          1

typedef struct {
  uint32_t            elementsPerBuffer;
  uint32_t            *h_reorderBuffer;
  uint32_t            *d_reorderBuffer;
} gpu_scheduler_t;

typedef struct {
  uint32_t            numBuckets;
  uint32_t            numWarps;
  uint32_t            *h_initPosPerBucket;
  uint32_t            *h_initWarpPerBucket;
  uint32_t			      *h_endPosPerBucket;
  uint32_t            *d_initPosPerBucket;
  uint32_t            *d_initWarpPerBucket;
  uint32_t			      *d_endPosPerBucket;
  //Re-scheduling: thread perspective
  gpu_scheduler_t     threadMapScheduler;
  //Re-scheduling: task perspective
  gpu_scheduler_t     taskMapScheduler;
} gpu_scheduler_buffer_t;



#endif /* GPU_SCHEDULER_H_ */
