/*
 * PROJECT: GEMMapper
 * FILE: search_pipeline.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_pipeline/search_pipeline.h"
#include "archive/archive_search_se.h"
#include "archive/archive_search_pe.h"

/*
 * Setup
 */
search_pipeline_t* search_pipeline_new(
    mapper_parameters_t* const mapper_parameters,
    gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const bool paired_end) {
  // Parameters
  const mapper_parameters_cuda_t* const cuda = &mapper_parameters->cuda;
  archive_t* const archive = mapper_parameters->archive;
  fm_index_t* const fm_index = archive->fm_index;
  const bool cpu_emulated = cuda->cpu_emulated;
  // Alloc
  search_pipeline_t* search_pipeline = mm_alloc(search_pipeline_t);
  // Allocate archive-search cache
  search_pipeline->archive_search_cache =
      archive_search_cache_new(archive,&mapper_parameters->base_search_parameters);
  // Allocate Support Data Structures
  search_pipeline->mm_stack = mm_stack_new(mm_pool_get_slab(mm_pool_32MB)); // Memory-Stack allocator
  search_pipeline->mapper_stats = mapper_stats_new(); // Mapping Statistics
  interval_set_init(&search_pipeline->interval_set); // Interval Set
  // Allocate pipeline stages
  uint64_t acc_buffers_offset = buffers_offset;
  search_pipeline->stage_region_profile = search_stage_region_profile_new(
      gpu_buffer_collection,acc_buffers_offset,
      cuda->num_fmi_bsearch_buffers,fm_index,
      cpu_emulated || !gpu_buffer_collection->region_profile);
  acc_buffers_offset += cuda->num_fmi_bsearch_buffers;
  search_pipeline->stage_decode_candidates = search_stage_decode_candidates_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_fmi_decode_buffers,fm_index,
      gpu_buffer_collection->decode_candidates_sa && !cpu_emulated,
      gpu_buffer_collection->decode_candidates_text && !cpu_emulated,
      search_pipeline->mm_stack);
  acc_buffers_offset += cuda->num_fmi_decode_buffers;
  search_pipeline->stage_verify_candidates = search_stage_verify_candidates_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_bpm_buffers,
      paired_end,cpu_emulated || !gpu_buffer_collection->verify_candidates,
      archive->text,search_pipeline->mm_stack);
  // Return
  return search_pipeline;
}
void search_pipeline_clear(search_pipeline_t* const search_pipeline) {
  mm_stack_clear(search_pipeline->mm_stack);
}
void search_pipeline_delete(search_pipeline_t* const search_pipeline) {
  search_stage_region_profile_delete(search_pipeline->stage_region_profile,search_pipeline->archive_search_cache);
  search_stage_decode_candidates_delete(search_pipeline->stage_decode_candidates,search_pipeline->archive_search_cache);
  search_stage_verify_candidates_delete(search_pipeline->stage_verify_candidates,search_pipeline->archive_search_cache);
  archive_search_cache_delete(search_pipeline->archive_search_cache);
  mm_stack_delete(search_pipeline->mm_stack);
  mapper_stats_delete(search_pipeline->mapper_stats);
  interval_set_destroy(&search_pipeline->interval_set);
  mm_free(search_pipeline);
}
/*
 * Archive-Search allocation
 */
void search_pipeline_allocate_se(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search) {
  // Alloc
  archive_search_cache_se_alloc(search_pipeline->archive_search_cache,archive_search);
  // Inject Support Data Structures
  archive_search_inject_mm_stack(*archive_search,search_pipeline->mm_stack);
  archive_search_inject_mapper_stats(*archive_search,search_pipeline->mapper_stats);
  archive_search_inject_interval_set(*archive_search,&search_pipeline->interval_set);
}
void search_pipeline_allocate_pe(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  // Alloc
  archive_search_cache_pe_alloc(
      search_pipeline->archive_search_cache,
      archive_search_end1,archive_search_end2);
  // Inject Support Data Structures
  archive_search_inject_mm_stack(*archive_search_end1,search_pipeline->mm_stack);
  archive_search_inject_mm_stack(*archive_search_end2,search_pipeline->mm_stack);
  archive_search_inject_mapper_stats(*archive_search_end1,search_pipeline->mapper_stats);
  archive_search_inject_mapper_stats(*archive_search_end2,search_pipeline->mapper_stats);
  archive_search_inject_interval_set(*archive_search_end1,&search_pipeline->interval_set);
  archive_search_inject_interval_set(*archive_search_end2,&search_pipeline->interval_set);
}

