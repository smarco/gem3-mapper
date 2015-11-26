/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer_collection.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu_buffer_collection.h"
#include "../resources/gpu_modules/gpu_interface.h"

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA
/*
 * GPU Support
 */
GEM_INLINE bool gpu_supported() {
  return (gpu_get_num_supported_devices_(GPU_ARCH_SUPPORTED) > 0);
}
/*
 * Setup
 */
GEM_INLINE gpu_buffer_collection_t* gpu_buffer_collection_new(
    archive_t* const archive,const uint64_t num_buffers,
    const uint64_t buffer_size,const bool verbose) {
  PROFILE_START(GP_GPU_BUFFER_COLLECTION_INIT);
  // Allocate Buffer Collection
  gpu_buffer_collection_t* const buffer_collection = mm_alloc(gpu_buffer_collection_t);
  buffer_collection->num_buffers = num_buffers;
  // Initialize GPU Runtime
  bwt_t* const bwt = archive->fm_index->bwt;
  gpu_buffers_dto_t* const gpu_buffers_dto = mm_alloc(gpu_buffers_dto_t);
  gpu_buffers_dto->buffer = NULL;
  gpu_buffers_dto->numBuffers = num_buffers;
  gpu_buffers_dto->maxMbPerBuffer = CONVERT_B_TO_MB(buffer_size);
  gpu_buffers_dto->activeModules = GPU_ALL_MODULES;
  buffer_collection->gpu_buffers_dto = gpu_buffers_dto;
  gpu_gem_fmi_dto_t gpu_gem_fmi_dto = {
      .c = bwt->c,
      .C = bwt->C,
      .mayor_counters = bwt->mayor_counters,
      .bwt_mem = bwt->bwt_mem,
      .bwt_length = bwt->length,
  };
  gpu_index_dto_t gpu_index_dto = {
      .fmi = &gpu_gem_fmi_dto,
      .indexCoding = GPU_INDEX_GEM_FULL,
      .bwtSize = bwt->length,
  };
  gpu_reference_dto_t gpu_reference_dto = {
      .reference = (char*)dna_text_get_text(archive->text->enc_text),
      .refCoding = (archive->text->explicit_complement) ? GPU_REF_GEM_FULL : GPU_REF_GEM_ONLY_FORWARD,
      .refSize = (archive->text->explicit_complement) ?
          dna_text_get_length(archive->text->enc_text) : archive->text->forward_text_length,
  };
  gpu_info_dto_t gpu_info_dto = {
      .selectedArchitectures = GPU_ARCH_SUPPORTED,
      .userAllocOption = GPU_LOCAL_OR_REMOTE_DATA,
  };
  gpu_init_buffers_(gpu_buffers_dto,&gpu_index_dto,&gpu_reference_dto,&gpu_info_dto,verbose);
  buffer_collection->internal_buffers = gpu_buffers_dto->buffer;
  // Return
  PROFILE_STOP(GP_GPU_BUFFER_COLLECTION_INIT);
  return buffer_collection;
}
GEM_INLINE void gpu_buffer_collection_delete(gpu_buffer_collection_t* const gpu_buffer_collection) {
  gpu_destroy_buffers_((gpu_buffers_dto_t*)gpu_buffer_collection->gpu_buffers_dto); // Destroy buffers
  mm_free(gpu_buffer_collection->gpu_buffers_dto); // Free DTO
  mm_free(gpu_buffer_collection); // Free Handler
}
/*
 * Accessors
 */
GEM_INLINE void* gpu_buffer_collection_get_buffer(
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_no) {
  return gpu_buffer_collection->internal_buffers[buffer_no];
}
/*
 * CUDA NOT-Supported
 */
#else
/*
 * GPU Support
 */
GEM_INLINE bool gpu_supported() { return false; }
/*
 * Setup
 */
GEM_INLINE gpu_buffer_collection_t* gpu_buffer_collection_new(
    archive_t* const archive,const uint64_t num_buffers,
    const uint64_t buffer_size,const bool verbose) {
  GEM_CUDA_NOT_SUPPORTED();
  return NULL;
}
GEM_INLINE void gpu_buffer_collection_delete(gpu_buffer_collection_t* const gpu_buffer_collection) { GEM_CUDA_NOT_SUPPORTED(); }
/*
 * Accessors
 */
GEM_INLINE void* gpu_buffer_collection_get_buffer(
    const gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_no) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
#endif /* HAVE_CUDA */
