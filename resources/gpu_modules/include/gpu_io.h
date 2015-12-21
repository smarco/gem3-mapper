/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Common headers and data structures for BPM on GPU library
 */

#ifndef GPU_IO_H_
#define GPU_IO_H_

#define _LARGEFILE_SOURCE   1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS   64

#include "gpu_commons.h"
/* Include the required objects */
#include "gpu_reference.h"
#include "gpu_index.h"

/* Input & Output reference functions  */
gpu_error_t gpu_load_BWT_MFASTA(const char* const fn, gpu_index_buffer_t* const fmi, char **h_BWT);
gpu_error_t gpu_load_index_PROFILE(const char* const fn, gpu_index_buffer_t* const index);
gpu_error_t gpu_save_index_PROFILE(const char* const fn, const gpu_index_buffer_t* const index);

gpu_error_t gpu_load_reference_MFASTA(const char *fn, gpu_reference_buffer_t* const reference);
gpu_error_t gpu_load_reference_PROFILE(const char* const fn, gpu_reference_buffer_t* const reference);
gpu_error_t gpu_save_reference_PROFILE(const char* const fn, const gpu_reference_buffer_t* const reference);

gpu_error_t gpu_load_index_GEM_FULL(const char *fn, gpu_index_buffer_t* const nindex);
gpu_error_t gpu_save_index_GEM_FULL(const char* const fn, const gpu_index_buffer_t* const index);

gpu_error_t gpu_load_reference_GEM_FULL(const char* const fn, gpu_reference_buffer_t* const reference);
gpu_error_t gpu_save_reference_GEM_FULL(const char* const fn, const gpu_reference_buffer_t* const reference);

#endif /* GPU_IO_H_ */
