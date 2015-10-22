/*
 * PROJECT: GEMMapper
 * FILE: search_group.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_group.h"
#include "archive_search_se.h"
#include "archive_search_pe.h"
#include "align_bpm_gpu.h"

/*
 * Setup
 */
GEM_INLINE search_group_t* search_group_new(
    mapper_parameters_t* const mapper_parameters,bpm_gpu_buffer_t* const bpm_gpu_buffers,
    const uint64_t num_verify_candidates_buffers) {
  // Alloc
  search_group_t* search_group = mm_alloc(search_group_t);
  // Allocate verification-candidates
  search_group->search_group_vc = search_group_verify_candidates_new(
      bpm_gpu_buffers,num_verify_candidates_buffers,mapper_parameters->cuda.cpu_emulated);
  // Allocate archive-search cache
  search_group->archive_search_cache = archive_search_cache_new(mapper_parameters);
  // Allocate mm-search
  search_group->mm_search = mm_search_new(mm_pool_get_slab(mm_pool_32MB));
  // Return
  return search_group;
}
GEM_INLINE void search_group_init(search_group_t* const search_group) {
  // Init buffers
  search_group_verify_candidates_init(search_group->search_group_vc);
}
GEM_INLINE void search_group_clear(search_group_t* const search_group) {
  search_group_verify_candidates_clear(search_group->search_group_vc,search_group->archive_search_cache);
  mm_search_clear(search_group->mm_search);
}
GEM_INLINE void search_group_delete(search_group_t* const search_group) {
  search_group_verify_candidates_delete(search_group->search_group_vc);
  archive_search_cache_delete(search_group->archive_search_cache);
  mm_search_delete(search_group->mm_search);
  mm_free(search_group);
}
/*
 * Archive-Search allocation
 */
GEM_INLINE archive_search_t* search_group_allocate_se(search_group_t* const search_group) {
  // Alloc
  archive_search_t* const archive_search = archive_search_cache_alloc(search_group->archive_search_cache);
  // Init archive search
  archive_search_se_configure(archive_search,search_group->mm_search);
  text_collection_clear(&search_group->mm_search->text_collection); // Clear text-collection
  // Return
  return archive_search;
}
GEM_INLINE void search_group_allocate_pe(
    search_group_t* const search_group,
    archive_search_t** const archive_search_end1,archive_search_t** const archive_search_end2) {
  // Alloca
  *archive_search_end1 = archive_search_cache_alloc(search_group->archive_search_cache);
  *archive_search_end2 = archive_search_cache_alloc(search_group->archive_search_cache);
  // Init archive search
  archive_search_pe_configure(*archive_search_end1,*archive_search_end2,search_group->mm_search);
  text_collection_clear(&search_group->mm_search->text_collection); // Clear text-collection
}
/*
 * Accessors
 */
GEM_INLINE text_collection_t* search_group_get_text_collection(search_group_t* const search_group) {
  return &search_group->mm_search->text_collection;
}

