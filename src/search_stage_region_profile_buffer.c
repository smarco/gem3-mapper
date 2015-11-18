/*
 * PROJECT: GEMMapper
 * FILE: search_stage_region_profile_buffer.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_stage_region_profile_buffer.h"
#include "archive_search_se_stepwise.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Setup
 */
GEM_INLINE search_stage_region_profile_buffer_t* search_stage_region_profile_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,const bool cpu_emulated) {
  // Alloc
  search_stage_region_profile_buffer_t* const region_profile_buffer = mm_alloc(search_stage_region_profile_buffer_t);
  // Init
  region_profile_buffer->gpu_buffer_fmi_search = gpu_buffer_fmi_search_new(gpu_buffer_collection,buffer_no);
  gpu_buffer_fmi_search_device(region_profile_buffer->gpu_buffer_fmi_search,cpu_emulated?DEVICE_CPU:DEVICE_GPU);
  const uint64_t max_queries = gpu_buffer_fmi_search_get_max_queries(region_profile_buffer->gpu_buffer_fmi_search);
  region_profile_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return region_profile_buffer;
}
GEM_INLINE void search_stage_region_profile_buffer_clear(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_cache_t* const archive_search_cache) {
  gpu_buffer_fmi_search_clear(region_profile_buffer->gpu_buffer_fmi_search);
  // Return searches to the cache
  if (archive_search_cache!=NULL) {
    VECTOR_ITERATE(region_profile_buffer->archive_searches,archive_search,n,archive_search_t*) {
      archive_search_cache_free(archive_search_cache,*archive_search);
    }
  }
  // Clear searches vector
  vector_clear(region_profile_buffer->archive_searches);
}
GEM_INLINE void search_stage_region_profile_buffer_delete(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_cache_t* const archive_search_cache) {
  gpu_buffer_fmi_search_delete(region_profile_buffer->gpu_buffer_fmi_search);
  // Return searches to the cache
  VECTOR_ITERATE(region_profile_buffer->archive_searches,archive_search,n,archive_search_t*) {
    archive_search_cache_free(archive_search_cache,*archive_search);
  }
  // Delete searches vector
  vector_delete(region_profile_buffer->archive_searches);
}
/*
 * Occupancy
 */
GEM_INLINE bool search_stage_region_profile_buffer_fits(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
// gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search = region_profile_buffer->gpu_buffer_fmi_search;
//  // Compute dimensions
//  uint64_t total_entries = 0,total_query_chunks = 0,total_candidate_chunks = 0;
//  gpu_buffer_fmi_search_compute_dimensions(gpu_buffer_align_bpm,
//      &archive_search_end1->forward_search_state.pattern,
//      archive_search_get_search_canditates(archive_search_end1),
//      &total_entries,&total_query_chunks,&total_candidate_chunks);
//  if (archive_search_end2!=NULL) {
//    gpu_buffer_fmi_search_compute_dimensions(gpu_buffer_align_bpm,
//        &archive_search_end2->forward_search_state.pattern,
//        archive_search_get_search_canditates(archive_search_end2),
//        &total_entries,&total_query_chunks,&total_candidate_chunks);
//  }
//  // Return if current search fits in buffer
//  return gpu_buffer_fmi_search_fits_in_buffer(gpu_buffer_align_bpm,
//      total_entries,total_query_chunks,total_candidate_chunks);
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  return false;
}
/*
 * Send/Receive
 */
GEM_INLINE void search_stage_region_profile_buffer_send(
    search_stage_region_profile_buffer_t* const region_profile_buffer) {
  gpu_buffer_fmi_search_send(region_profile_buffer->gpu_buffer_fmi_search);
}
GEM_INLINE void search_stage_region_profile_buffer_receive(
    search_stage_region_profile_buffer_t* const region_profile_buffer) {
  gpu_buffer_fmi_search_receive(region_profile_buffer->gpu_buffer_fmi_search);
}
/*
 * Accessors
 */
GEM_INLINE void search_stage_region_profile_buffer_add(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_t* const archive_search) {
  // Add archive-search
  vector_insert(region_profile_buffer->archive_searches,archive_search,archive_search_t*);
  // Copy profile-partitions to the buffer
  gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search = region_profile_buffer->gpu_buffer_fmi_search;
  archive_search_se_stepwise_region_profile_copy(archive_search,gpu_buffer_fmi_search);
}
GEM_INLINE void search_stage_region_profile_buffer_retrieve(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    const uint64_t search_idx,archive_search_t** const archive_search) {
  // Retrieve archive-search
  *archive_search = *vector_get_elm(region_profile_buffer->archive_searches,search_idx,archive_search_t*);
  // Retrieve searched profile-partitions from the buffer
  gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search = region_profile_buffer->gpu_buffer_fmi_search;
  archive_search_se_stepwise_region_profile_retrieve(*archive_search,gpu_buffer_fmi_search);
}


