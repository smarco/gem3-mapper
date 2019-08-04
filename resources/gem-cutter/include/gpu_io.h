/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_IO_H_
#define GPU_IO_H_

/* Unix-Linux portability  */
#ifndef O_BINARY
  #define O_BINARY  0
  #define O_TEXT    0
#endif

#define _LARGEFILE_SOURCE   1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS   64
#define GPU_FILE_BASIC_MODE (O_BINARY | O_CREAT)
#define GPU_FILE_PERMISIONS (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP)
#define off64_t off_t

#include "gpu_commons.h"
/* Include the required objects */
#include "gpu_reference.h"
#include "gpu_index.h"

/* I/O Primitives to read/write in buffered files */
gpu_error_t gpu_io_read_buffered(int fp, void* const buffer, const size_t bytesRequest);
gpu_error_t gpu_io_write_buffered(int fp, void* const buffer, const size_t bytesRequest);

/* Input & Output Multi-FASTA functions (Indexes) */
gpu_error_t gpu_io_load_specs_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const index, char **h_BWT);
/* Input & Output Multi-FASTA functions (Reference) */
gpu_error_t gpu_io_load_reference_specs_MFASTA(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_reference_MFASTA(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);

/* Input & Output BENCHMARK PROFILE functions (Indexes) */
gpu_error_t gpu_io_save_index_PROFILE(const char* const fn, const gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_index_specs_PROFILE(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_index_PROFILE(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
/* Input & Output BENCHMARK PROFILE functions (Reference) */
gpu_error_t gpu_io_load_reference_specs_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_reference_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_save_reference_PROFILE(const char* const fn, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);

/* Input & Output GEM-CUDA functions (Indexes) */
gpu_error_t gpu_io_load_index_specs_GEM_FULL(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_index_GEM_FULL(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_save_index_GEM_FULL(const char* const fn, const gpu_index_buffer_t* const index, const gpu_module_t activeModules);
/* Input & Output GEM-CUDA functions (Reference) */
gpu_error_t gpu_io_load_reference_specs_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_reference_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_save_reference_GEM_FULL(const char* const fn, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);

/* Local definitions */
gpu_error_t gpu_io_save_module_info_GEM_FULL(const int fp, const gpu_module_t fileActiveModules);
gpu_error_t gpu_io_load_module_info_GEM_FULL(const int fp, gpu_module_t* const fileActiveModules);
gpu_error_t gpu_io_save_offsets_info_GEM_FULL(const int fp, const off64_t fileOffsetFMIndex, const off64_t fileOffsetSAIndex, const off64_t fileOffsetRef);
gpu_error_t gpu_io_load_offsets_info_GEM_FULL(const int fp, off64_t* const fileOffsetFMIndex, off64_t* const fileOffsetSAIndex, off64_t* const fileOffsetRef);

#endif /* GPU_IO_H_ */
