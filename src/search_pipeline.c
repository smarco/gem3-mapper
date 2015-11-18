/*
 * PROJECT: GEMMapper
 * FILE: search_pipeline.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_pipeline.h"
#include "archive_search_se.h"
#include "archive_search_pe.h"

/*
 * Setup
 */
GEM_INLINE search_pipeline_t* search_pipeline_new(
    mapper_parameters_t* const mapper_parameters,
    gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset) {
  // Parameters
  const mapper_parameters_cuda_t* const cuda = &mapper_parameters->cuda;
  const bool cpu_emulated = cuda->cpu_emulated;
  // Alloc
  search_pipeline_t* search_pipeline = mm_alloc(search_pipeline_t);
  // Allocate pipeline stages
  uint64_t acc_buffers_offset = buffers_offset;
  search_pipeline->stage_region_profile = search_stage_region_profile_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_fmi_bsearch_buffers,cpu_emulated);
  acc_buffers_offset += cuda->num_fmi_bsearch_buffers;
  search_pipeline->stage_decode_candidates = search_stage_decode_candidates_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_fmi_decode_buffers,cpu_emulated);
  acc_buffers_offset += cuda->num_fmi_decode_buffers;
  search_pipeline->stage_verify_candidates = search_stage_verify_candidates_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_bpm_buffers,cpu_emulated);
  // Allocate archive-search cache
  search_pipeline->archive_search_cache = archive_search_cache_new(mapper_parameters);
  // Allocate mm-search
  search_pipeline->mm_search = mm_search_new(mm_pool_get_slab(mm_pool_32MB));
  // Return
  return search_pipeline;
}
GEM_INLINE void search_pipeline_clear(search_pipeline_t* const search_pipeline) {
  search_stage_region_profile_clear(search_pipeline->stage_region_profile,search_pipeline->archive_search_cache);
  search_stage_decode_candidates_clear(search_pipeline->stage_decode_candidates,search_pipeline->archive_search_cache);
  search_stage_verify_candidates_clear(search_pipeline->stage_verify_candidates,search_pipeline->archive_search_cache);
  mm_search_clear(search_pipeline->mm_search);
}
GEM_INLINE void search_pipeline_delete(search_pipeline_t* const search_pipeline) {
  search_stage_region_profile_delete(search_pipeline->stage_region_profile,search_pipeline->archive_search_cache);
  search_stage_decode_candidates_delete(search_pipeline->stage_decode_candidates,search_pipeline->archive_search_cache);
  search_stage_verify_candidates_delete(search_pipeline->stage_verify_candidates,search_pipeline->archive_search_cache);
  archive_search_cache_delete(search_pipeline->archive_search_cache);
  mm_search_delete(search_pipeline->mm_search);
  mm_free(search_pipeline);
}
/*
 * Archive-Search allocation
 */
GEM_INLINE archive_search_t* search_pipeline_allocate_se(search_pipeline_t* const search_pipeline) {
  // Alloc
  archive_search_t* const archive_search = archive_search_cache_alloc(search_pipeline->archive_search_cache);
  // Init archive search
  archive_search_se_configure(archive_search,search_pipeline->mm_search);
  text_collection_clear(&search_pipeline->mm_search->text_collection); // Clear text-collection
  // Return
  return archive_search;
}
GEM_INLINE void search_pipeline_allocate_pe(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  // Alloca
  *archive_search_end1 = archive_search_cache_alloc(search_pipeline->archive_search_cache);
  *archive_search_end2 = archive_search_cache_alloc(search_pipeline->archive_search_cache);
  // Init archive search
  archive_search_pe_configure(*archive_search_end1,*archive_search_end2,search_pipeline->mm_search);
  text_collection_clear(&search_pipeline->mm_search->text_collection); // Clear text-collection
}
/*
 * Accessors
 */
GEM_INLINE text_collection_t* search_pipeline_get_text_collection(search_pipeline_t* const search_pipeline) {
  return &search_pipeline->mm_search->text_collection;
}

