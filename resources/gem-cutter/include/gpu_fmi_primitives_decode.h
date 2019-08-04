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
#include "gpu_fmi_structure.h"
#include "gpu_sa_index.h"

#ifndef GPU_FMI_PRIMITIVES_DECODE_H_
#define GPU_FMI_PRIMITIVES_DECODE_H_

/********************************
Common constants for Device & Host
*********************************/

/* Defines related to FMI DECODE primitives */
#define GPU_FMI_ENTRIES_PER_DECODE          2                                                           // Necessary FMI entries for unaligned FMI threads
#define GPU_FMI_THREADS_PER_COUNTERS        2                                                           // Necessary threads to load all the counters

#define GPU_FMI_DECODE_THREADS_PER_ENTRY    (GPU_FMI_THREADS_PER_ENTRY * GPU_FMI_ENTRIES_PER_DECODE)    // FMI_THREADS_PER_ENTRY + 1 (rounded to pow of 2)
#define GPU_FMI_DECODE_THREADS_PER_LOAD     (GPU_FMI_THREADS_PER_ENTRY + GPU_FMI_THREADS_PER_COUNTERS)  // Threads to load an entry for the decode primitive
#define GPU_FMI_DECODE_POS_BUFFER_PADDING   10
#define GPU_FMI_DECODE_MIN_ELEMENTS         (2048 / GPU_FMI_DECODE_THREADS_PER_ENTRY)                   // Min elements per buffer (related to the SM)


/*****************************
Internal Objects (Decode)
*****************************/

typedef struct {
  uint32_t                    numDecodings;
  gpu_fmi_decode_init_pos_t   *h_initBWTPos;
  gpu_fmi_decode_init_pos_t   *d_initBWTPos;
} gpu_fmi_decode_init_pos_buffer_t;

typedef struct {
  uint32_t                    numDecodings;
  gpu_fmi_decode_end_pos_t    *h_endBWTPos;
  gpu_fmi_decode_end_pos_t    *d_endBWTPos;
} gpu_fmi_decode_end_pos_buffer_t;

typedef struct {
  uint32_t                    numDecodings;
  gpu_fmi_decode_text_pos_t   *h_textPos;
  gpu_fmi_decode_text_pos_t   *d_textPos;
} gpu_fmi_decode_text_pos_buffer_t;

/*****************************
Internal Objects (General)
*****************************/

typedef struct {
  uint32_t                          numMaxInitPositions;
  uint32_t                          numMaxEndPositions;
  uint32_t                          numMaxTextPositions;
  uint32_t                          samplingRate;
  gpu_fmi_decode_init_pos_buffer_t  initPositions;
  gpu_fmi_decode_end_pos_buffer_t   endPositions;
  gpu_fmi_decode_text_pos_buffer_t  textPositions;
} gpu_fmi_decode_buffer_t;

#include "gpu_buffer.h"

/* Functions to init the buffers (DECODE) */
size_t      gpu_fmi_decode_input_size();
void        gpu_fmi_decode_reallocate_host_buffer_layout(gpu_buffer_t* const mBuff);
void        gpu_fmi_decode_reallocate_device_buffer_layout(gpu_buffer_t* const mBuff);

/* Functions to transfer data HOST <-> DEVICE (DECODE) */
gpu_error_t gpu_fmi_decode_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff);
gpu_error_t gpu_fmi_decode_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff);

/* DEVICE Kernels */
gpu_error_t gpu_fmi_decode_process_buffer(gpu_buffer_t* const mBuff);


#endif /* GPU_FMI_PRIMITIVES_DECODE_H_ */

