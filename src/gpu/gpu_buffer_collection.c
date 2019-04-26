/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 * DESCRIPTION:
 *   GPU-adaptor module provides support data-structures and functions
 *   to handle the intermediate buffers that offload computation to the GPU(s)
 */

#include "gpu/gpu_buffer_collection.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA
/*
 * Setup
 */
gpu_buffer_collection_t* gpu_buffer_collection_new(
    char* const gpu_index_name,
    const uint64_t num_buffers,
    const uint64_t buffer_size,
    const uint64_t gpu_devices,
    const bool verbose) {
  PROFILE_START(GP_GPU_BUFFER_COLLECTION_INIT,PROFILE_LEVEL);
  // Allocate Buffer Collection
  gpu_buffer_collection_t* const buffer_collection = mm_alloc(gpu_buffer_collection_t);
  buffer_collection->num_buffers = num_buffers;
  // Initialize GPU Runtime
  gpu_buffers_dto_t* const gpu_buffers_dto = mm_alloc(gpu_buffers_dto_t);
  gpu_buffers_dto->buffer = NULL;
  gpu_buffers_dto->numBuffers = num_buffers;
  gpu_buffers_dto->maxMbPerBuffer = buffer_size;
  gpu_buffers_dto->activeModules =
      GPU_FMI_ADAPT_SEARCH |
      GPU_FMI_DECODE_POS |
      GPU_SA_DECODE_POS |
      GPU_KMER_FILTER |
      GPU_BPM_FILTER |
      GPU_BPM_ALIGN;
  buffer_collection->gpu_buffers_dto = gpu_buffers_dto;
  gpu_index_dto_t gpu_index_dto = {
    .filename             = gpu_index_name,
    //Initialize Table for FM-Index
    .fmi.h_offsetsTable   = NULL,
    .fmi.h_table          = NULL,
    .fmi.numLevelsTable   = 0,
    .fmi.skipLevelsTable  = 0,
    .fmi.numElementsTable = 0,
    //Initialize FM-Index
    .fmi.h_plain          = NULL,
    .fmi.h_fmi            = NULL,
    .fmi.bwtSize          = 0,
    .fmi.indexCoding      = GPU_INDEX_GEM_FILE,
    //Initialize Suffix-Array
    .sa.h_plain           = NULL,
    .sa.h_sa              = NULL,
    .sa.numEntries        = 0,
    .sa.samplingRate      = 0,
    .sa.indexCoding       = GPU_INDEX_GEM_FILE,
  };
  gpu_reference_dto_t gpu_reference_dto = {
    .reference            = gpu_index_name,
    .refCoding            = GPU_REF_GEM_FILE,
    .refSize              = 0,
  };
  gpu_info_dto_t gpu_info_dto = {
    .selectedArchitectures = GPU_ARCH_SUPPORTED,
    .userAllocOption       = GPU_GEM_POLICY,
    .activatedModules      = GPU_NONE_MODULES,
    .allocatedStructures   = GPU_NONE_MODULES,
    .verbose               = verbose,
  };
  gpu_init_buffers_(gpu_buffers_dto,&gpu_index_dto,&gpu_reference_dto,&gpu_info_dto);
  buffer_collection->internal_buffers = gpu_buffers_dto->buffer;
  // Set active modules
  buffer_collection->gpu_region_profile_available = gpu_info_dto.activatedModules & GPU_FMI_ADAPT_SEARCH;
  buffer_collection->gpu_decode_sa_available = gpu_info_dto.activatedModules & GPU_FMI_DECODE_POS;
  buffer_collection->gpu_decode_text_available = gpu_info_dto.activatedModules & GPU_SA_DECODE_POS;
  /*if (buffer_collection->gpu_decode_sa_available !=
      buffer_collection->gpu_decode_text_available) {
    buffer_collection->gpu_decode_sa_available   = false; // Decode goes all or none
    buffer_collection->gpu_decode_text_available = false; // Decode goes all or none
  }*/
  buffer_collection->gpu_kmer_filter_available = gpu_info_dto.activatedModules & GPU_KMER_FILTER;
  buffer_collection->gpu_bpm_distance_available = gpu_info_dto.activatedModules & GPU_BPM_FILTER;
  buffer_collection->gpu_bpm_align_available = gpu_info_dto.activatedModules & GPU_BPM_ALIGN;
  // Return
  PROFILE_STOP(GP_GPU_BUFFER_COLLECTION_INIT,PROFILE_LEVEL);
  return buffer_collection;
}
void gpu_buffer_collection_delete(
    gpu_buffer_collection_t* const gpu_buffer_collection) {
  gpu_destroy_buffers_((gpu_buffers_dto_t*)gpu_buffer_collection->gpu_buffers_dto); // Destroy buffers
  mm_free(gpu_buffer_collection->gpu_buffers_dto); // Free DTO
  mm_free(gpu_buffer_collection); // Free Handler
}
/*
 * Accessors
 */
void* gpu_buffer_collection_get_buffer(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no) {
  return gpu_buffer_collection->internal_buffers[buffer_no];
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * Setup
 */
gpu_buffer_collection_t* gpu_buffer_collection_new(
    char* const gpu_index_name,
    const uint64_t num_buffers,
    const uint64_t buffer_size,
    const uint64_t gpu_devices,
    const bool verbose) {
  GEM_CUDA_NOT_SUPPORTED();
  return NULL;
}
void gpu_buffer_collection_delete(
    gpu_buffer_collection_t* const gpu_buffer_collection) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Accessors
 */
void* gpu_buffer_collection_get_buffer(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
#endif /* HAVE_CUDA */
