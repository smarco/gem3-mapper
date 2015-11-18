/*
 * PROJECT: GEMMapper
 * FILE: gpu_buffer.c
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
  return (gpu_buffer_align_bpm_get_num_supported_devices_() > 0);
}
/*
 * Setup
 */
GEM_INLINE gpu_buffer_collection_t* gpu_buffer_collection_new(
    archive_t* const archive,const uint64_t num_buffers,
    const uint64_t buffer_size,const bool verbose) {
  PROF_START(GP_BPM_GPU_INIT);
  // Allocate Buffer Collection
  gpu_buffer_collection_t* const buffer_collection = mm_alloc(gpu_buffer_collection_t);
  buffer_collection->num_buffers = num_buffers;
  // Parameters
  const char* const text = (const char* const) dna_text_get_text(archive->text->enc_text);
  const uint64_t text_length = dna_text_get_length(archive->text->enc_text);
  const gpu_ref_coding_t text_encoding =
      (archive->text->explicit_complement) ? GPU_REF_GEM_FULL : GPU_REF_GEM_ONLY_FORWARD;
  bwt_t* const bwt = archive->fm_index->bwt;
  const uint64_t bwt_length = bwt->length;
  gpu_fmi_gem_dto_t gpu_fmi_gem_dto = {
      .mayor_counters = bwt->mayor_counters,
      .bwt_mem = bwt->bwt_mem,
  }
  const gpu_data_location_t gpu_data_location = GPU_LOCAL_DATA;
  // Initialize GPU Runtime
  gpu_init_buffers_(&buffer_collection->internal_buffers,num_buffers,
      CONVERT_B_TO_MB(buffer_size),text,text_encoding,text_encoding,
      &gpu_fmi_gem_dto,GPU_INDEX_GEM_FULL,bwt_length,
      ARCH_SUPPORTED,LOCAL_OR_REMOTE_REFERENCE,verbose);
  // Return
  PROF_STOP(GP_BPM_GPU_INIT);
  return buffer_collection;
}
GEM_INLINE void gpu_buffer_collection_delete(gpu_buffer_collection_t* const gpu_buffer_collection) {
  gpu_destroy_(&gpu_buffer_collection->internal_buffers); // Destroy buffers
  mm_free(gpu_buffer_collection); // Free Handler
}
/*
 * Accessors
 */
GEM_INLINE void* gpu_buffer_collection_get_buffer(
    gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_no) {
  return gpu_buffer_collection->internal_buffers + buffer_no;
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
    gpu_buffer_collection_t* const gpu_buffer_collection,const uint64_t buffer_no) { GEM_CUDA_NOT_SUPPORTED(); return NULL; }
#endif /* HAVE_CUDA */
