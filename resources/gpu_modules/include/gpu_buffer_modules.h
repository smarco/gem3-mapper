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
#include "gpu_bpm_primitives_filter.h"
#include "gpu_bpm_primitives_align.h"
#include "gpu_kmer_primitives_filter.h"

#ifndef GPU_BUFFER_MODULES_H_
#define GPU_BUFFER_MODULES_H_

typedef union{
  gpu_bpm_align_buffer_t   abpm;
  gpu_bpm_filter_buffer_t  fbpm;
  gpu_kmer_filter_buffer_t fkmer;
  gpu_fmi_asearch_buffer_t asearch;
  gpu_fmi_ssearch_buffer_t ssearch;
  gpu_fmi_decode_buffer_t  decode;
} gpu_buffer_modules_t;

#endif /* GPU_BUFFER_MODULES_H_ */
