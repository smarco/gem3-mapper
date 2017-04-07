/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */


#ifndef GPU_SCHEDULER_H_
#define GPU_SCHEDULER_H_

typedef struct {
  uint32_t            numBuckets;
  uint32_t            elementsPerBuffer;
  uint32_t            numWarps;
  uint32_t            *h_reorderBuffer;
  uint32_t            *d_reorderBuffer;
  uint32_t            *h_initPosPerBucket;
  uint32_t            *h_initWarpPerBucket;
  uint32_t            *d_initPosPerBucket;
  uint32_t            *d_initWarpPerBucket;
} gpu_scheduler_buffer_t;



#endif /* GPU_SCHEDULER_H_ */
