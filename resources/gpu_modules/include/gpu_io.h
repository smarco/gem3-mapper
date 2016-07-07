/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_IO_H_
#define GPU_IO_H_

#define _LARGEFILE_SOURCE   1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS   64
#define off64_t off_t

#include "gpu_commons.h"
/* Include the required objects */
#include "gpu_reference.h"
#include "gpu_index.h"


/* Input & Output Multi-FASTA functions (Indexes) */
gpu_error_t gpu_io_load_specs_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const index, char **h_BWT);
/* Input & Output Multi-FASTA functions (Reference) */
gpu_error_t gpu_io_load_reference_specs_MFASTA(const char *fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_reference_MFASTA(const char *fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);

/* Input & Output BENCHMARK PROFILE functions (Indexes) */
gpu_error_t gpu_io_save_index_PROFILE(const char* const fn, const gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_index_specs_PROFILE(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_index_PROFILE(const char* const fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
/* Input & Output BENCHMARK PROFILE functions (Reference) */
gpu_error_t gpu_io_load_reference_specs_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_reference_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_save_reference_PROFILE(const char* const fn, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);

/* Input & Output GEM-CUDA functions (Indexes) */
gpu_error_t gpu_io_load_index_specs_GEM_FULL(const char *fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_index_GEM_FULL(const char *fn, gpu_index_buffer_t* const index, const gpu_module_t activeModules);
gpu_error_t gpu_io_save_index_GEM_FULL(const char* const fn, const gpu_index_buffer_t* const index, const gpu_module_t activeModules);
/* Input & Output GEM-CUDA functions (Reference) */
gpu_error_t gpu_io_load_reference_specs_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_load_reference_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);
gpu_error_t gpu_io_save_reference_GEM_FULL(const char* const fn, const gpu_reference_buffer_t* const reference, const gpu_module_t activeModules);


#endif /* GPU_IO_H_ */
