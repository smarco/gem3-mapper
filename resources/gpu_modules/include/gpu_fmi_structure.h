/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_STRUCTURE_H_
#define GPU_FMI_STRUCTURE_H_

/* Defines related to BWT representation */
#define GPU_FMI_BITS_LOAD_PER_THREAD  128     // accesses of 16 Bytes / thread

#define GPU_FMI_ENTRY_LENGTH        ((GPU_FMI_COUNTERS_PER_ENTRY * GPU_UINT64_LENGTH) \
                                    + (GPU_FMI_ENTRY_SIZE * GPU_FMI_BWT_CHAR_LENGTH))       //   512 bits
#define GPU_FMI_THREADS_PER_ENTRY   (GPU_FMI_ENTRY_LENGTH / GPU_FMI_BITS_LOAD_PER_THREAD)   //   4 threads
#define GPU_FMI_ENTRIES_PER_WARP    (GPU_WARP_SIZE / GPU_FMI_THREADS_PER_ENTRY)
#define GPU_FMI_ENTRIES_PER_BLOCK   (GPU_THREADS_PER_BLOCK / GPU_FMI_THREADS_PER_ENTRY)
#define GPU_FMI_SLICES_PER_THREAD   (GPU_FMI_BITS_LOAD_PER_THREAD / GPU_UINT32_LENGTH)      // slices of 16 Bytes / thread


/*************************************
Specific types for the Devices (GPUs)
**************************************/

// Just to cast and exchange memory between threads for the shared mem
typedef union {
  uint4 s;
  uint32_t v[GPU_FMI_SLICES_PER_THREAD];
} gpu_fmi_exch_bmp_mem_t;

typedef union {
  uint4 v[GPU_FMI_THREADS_PER_ENTRY];
} gpu_fmi_device_entry_t;

typedef struct {
  uint32_t bitmaps[GPU_FMI_BWT_CHAR_LENGTH];
} gpu_index_bitmap_entry_t;

typedef struct {
  uint64_t counters[GPU_FMI_NUM_COUNTERS];
} gpu_index_counter_entry_t;


#endif /* GPU_FMI_STRUCTURE_H_ */
