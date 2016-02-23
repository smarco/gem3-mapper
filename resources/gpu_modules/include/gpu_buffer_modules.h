/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

/* Include all the supported modules */
#include "gpu_fmi_primitives.h"
#include "gpu_sa_primitives.h"
#include "gpu_bpm_primitives.h"

#ifndef GPU_BUFFER_MODULES_H_
#define GPU_BUFFER_MODULES_H_

typedef union{
  gpu_bpm_buffer_t        bpm;
  gpu_fmi_search_buffer_t search;
  gpu_fmi_decode_buffer_t decode;
} gpu_buffer_modules_t;

#endif /* GPU_BUFFER_MODULES_H_ */
