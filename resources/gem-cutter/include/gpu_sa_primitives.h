/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include "gpu_index_modules.h"
#include "gpu_commons.h"
#include "gpu_sa_index.h"

#ifndef GPU_SA_PRIMITIVES_H_
#define GPU_SA_PRIMITIVES_H_

typedef struct {
  uint32_t                  numIntervals;
  gpu_sa_search_inter_t     *h_intervals;
  gpu_sa_search_inter_t     *d_intervals;
} gpu_sa_search_sa_inter_buffer_t;

typedef struct {
  uint32_t                   numDecodings;
  gpu_sa_decode_text_pos_t   *h_textPos;
  gpu_sa_decode_text_pos_t   *d_textPos;
} gpu_sa_decode_text_pos_buffer_t;

#include "gpu_buffer.h"


gpu_sa_entry_t* gpu_sa_buffer_get_index_(const void* const saBuffer);

/* DEVICE Kernels */
gpu_error_t gpu_sa_decode_process_buffer(gpu_buffer_t* const saBuffer);


#endif /* GPU_SA_PRIMITIVES_H_ */
