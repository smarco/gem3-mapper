/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
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


/* Functions to init the buffers (DECODE) */
/*
size_t gpu_fmi_buffer_input_size(gpu_buffer_t* const mBuff)
{
  size_t gpu_fmi_asearch_size_per_query(const uint32_t averageQuerySize, const uint32_t averageRegionsPerQuery)
}
void gpu_fmi_buffer_reallocate_host_layout(gpu_buffer_t* const mBuff)
{

}
void gpu_fmi_buffer_reallocate_device_layout(gpu_buffer_t* const mBuff)
{
}*/

/* Functions to transfer data HOST <-> DEVICE (DECODE) */
/*
gpu_error_t gpu_fmi_buffer_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff)
{

}
gpu_error_t gpu_fmi_buffer_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff)
{

}*/

/* DEVICE Kernels */
/*
gpu_error_t gpu_fmi_buffer_process(gpu_buffer_t* const mBuff)
{

}*/


#endif /* GPU_FMI_PRIMITIVES_C_ */

