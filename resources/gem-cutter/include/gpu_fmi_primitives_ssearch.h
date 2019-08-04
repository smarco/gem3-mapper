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

#ifndef GPU_FMI_PRIMITIVES_SSEARCH_H_
#define GPU_FMI_PRIMITIVES_SSEARCH_H_

/********************************
Common constants for Device & Host
*********************************/

/* Defines related to FMI static BACKWARD-SEARCH primitives */
#define GPU_FMI_ENTRIES_PER_SEED            2
#define GPU_FMI_SEED_ENTRY_LENGTH           128                                                         // 128 bits
#define GPU_FMI_SEED_FIELD_SIZE             8                                                           // 8 bits
#define GPU_FMI_SEED_MAX_CHARS              60                                                          // 60 bases
#define GPU_FMI_SEED_CHAR_LENGTH            2                                                           // 2 bits

#define GPU_FMI_SEED_THREADS_PER_ENTRY      (GPU_FMI_THREADS_PER_ENTRY * GPU_FMI_ENTRIES_PER_SEED)
#define GPU_FMI_SEED_ENTRIES_PER_WARP       (GPU_WARP_SIZE / GPU_FMI_SEED_THREADS_PER_ENTRY)
#define GPU_FMI_SEED_BASES_PER_ENTRY        (GPU_UINT64_LENGTH / GPU_FMI_SEED_CHAR_LENGTH)
#define GPU_FMI_SEARCH_SEEDS_BUFFER_PADDING 10
#define GPU_FMI_SEARCH_MIN_ELEMENTS         (2048 / GPU_FMI_SEED_THREADS_PER_ENTRY)                     // Min elements per buffer (related to the SM -2048th-)


/*****************************
Internal Objects (Static Search)
*****************************/

typedef struct {
  uint32_t                    numSeeds;
  gpu_fmi_search_seed_t       *d_seeds;
  gpu_fmi_search_seed_t       *h_seeds;
} gpu_fmi_search_seeds_buffer_t;

typedef struct {
  uint32_t                    numIntervals;
  gpu_sa_search_inter_t       *h_intervals;
  gpu_sa_search_inter_t       *d_intervals;
} gpu_fmi_search_sa_inter_buffer_t;


/*****************************
Internal Objects (General)
*****************************/

typedef struct {
  uint32_t                          numMaxSeeds;
  uint32_t                          numMaxIntervals;
  gpu_fmi_search_seeds_buffer_t     seeds;
  gpu_fmi_search_sa_inter_buffer_t  saIntervals;
} gpu_fmi_ssearch_buffer_t;

#include "gpu_buffer.h"

/* Functions to init the buffers (E. SEARCH) */
size_t      gpu_fmi_ssearch_input_size();
void        gpu_fmi_ssearch_reallocate_host_buffer_layout(gpu_buffer_t* mBuff);
void        gpu_fmi_ssearch_reallocate_device_buffer_layout(gpu_buffer_t* mBuff);

/* Functions to transfer data HOST <-> DEVICE (E. SEARCH) */
gpu_error_t gpu_fmi_ssearch_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff);
gpu_error_t gpu_fmi_ssearch_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff);

/* DEVICE Kernels */
gpu_error_t gpu_fmi_ssearch_process_buffer(gpu_buffer_t* const mBuff);

/* DEBUG */
uint32_t    gpu_fmi_ssearch_print_buffer(const void* const mBuff);
uint32_t    gpu_fmi_ssearch_print_seed(const gpu_fmi_search_seed_t seed, const uint32_t seedSize);
char        gpu_fmi_ssearch_bin_to_char(const uint32_t base);


#endif /* GPU_FMI_PRIMITIVES_SSEARCH_H_ */

