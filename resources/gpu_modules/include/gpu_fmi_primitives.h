/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include "gpu_index_modules.h"
#include "gpu_commons.h"
#include "gpu_fmi_structure.h"
#include "gpu_sa_index.h"

//FMI Modules
#include "gpu_fmi_primitives_asearch.h"
#include "gpu_fmi_primitives_ssearch.h"
#include "gpu_fmi_primitives_decode.h"

#ifndef GPU_FMI_PRIMITIVES_H_
#define GPU_FMI_PRIMITIVES_H_

/* Functions to init the buffers (DECODE) */
size_t      gpu_fmi_buffer_input_size(gpu_buffer_t* const mBuff);
void        gpu_fmi_buffer_reallocate_host_layout(gpu_buffer_t* const mBuff);
void        gpu_fmi_buffer_reallocate_device_layout(gpu_buffer_t* const mBuff);

/* Functions to transfer data HOST <-> DEVICE (DECODE) */
gpu_error_t gpu_fmi_buffer_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff);
gpu_error_t gpu_fmi_buffer_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff);

/* DEVICE Kernels */
gpu_error_t gpu_fmi_buffer_process(gpu_buffer_t* const mBuff);

#endif /* GPU_FMI_PRIMITIVES_H_ */

