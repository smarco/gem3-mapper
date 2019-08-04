/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_PRIMITIVES_C_
#define GPU_FMI_PRIMITIVES_C_

#include "../include/gpu_fmi_primitives.h"
#include "../include/gpu_sa_primitives.h"

/************************************************************
Functions to get the GPU FMI buffers
************************************************************/

gpu_fmi_entry_t* gpu_fmi_buffer_get_index_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->index->fmi.h_fmi);
}

#endif /* GPU_FMI_PRIMITIVES_C_ */

