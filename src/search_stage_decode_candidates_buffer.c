/*
 * PROJECT: GEMMapper
 * FILE: search_stage_decode_candidates_buffer.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_stage_decode_candidates_buffer.h"
#include "archive_search_se_stepwise.h"

/*
 * Setup
 */
GEM_INLINE search_stage_decode_candidates_buffer_t* search_stage_decode_candidates_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,const bool cpu_emulated) {
  // Alloc
  search_stage_decode_candidates_buffer_t* const decode_candidates_buffer =
      mm_alloc(search_stage_decode_candidates_buffer_t);
  // Init
  decode_candidates_buffer->gpu_buffer_fmi_decode = gpu_buffer_fmi_decode_new(gpu_buffer_collection,buffer_no);
  gpu_buffer_fmi_decode_device(decode_candidates_buffer->gpu_buffer_fmi_decode,cpu_emulated?DEVICE_CPU:DEVICE_GPU);
  const uint64_t max_queries = gpu_buffer_fmi_decode_get_max_queries(decode_candidates_buffer->gpu_buffer_fmi_decode);
  decode_candidates_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return decode_candidates_buffer;
}
GEM_INLINE void search_stage_decode_candidates_buffer_clear(
    search_stage_decode_candidates_buffer_t* const decode_candidates_buffer,
    archive_search_cache_t* const archive_search_cache) {
  gpu_buffer_fmi_decode_clear(decode_candidates_buffer->gpu_buffer_fmi_decode);
  // Return searches to the cache
  if (archive_search_cache!=NULL) {
    VECTOR_ITERATE(decode_candidates_buffer->archive_searches,archive_search,n,archive_search_t*) {
      archive_search_cache_free(archive_search_cache,*archive_search);
    }
  }
  // Clear searches vector
  vector_clear(decode_candidates_buffer->archive_searches);
}
GEM_INLINE void search_stage_decode_candidates_buffer_delete(
    search_stage_decode_candidates_buffer_t* const decode_candidates_buffer,
    archive_search_cache_t* const archive_search_cache) {
  gpu_buffer_fmi_decode_delete(decode_candidates_buffer->gpu_buffer_fmi_decode);
  // Return searches to the cache
  VECTOR_ITERATE(decode_candidates_buffer->archive_searches,archive_search,n,archive_search_t*) {
    archive_search_cache_free(archive_search_cache,*archive_search);
  }
  // Delete searches vector
  vector_delete(decode_candidates_buffer->archive_searches);
}
/*
 * Occupancy
 */
GEM_INLINE bool search_stage_decode_candidates_buffer_fits(
    search_stage_decode_candidates_buffer_t* const decode_candidates_buffer,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
// gpu_buffer_fmi_search_t* const gpu_buffer_fmi_decode = decode_candidates_buffer->gpu_buffer_fmi_decode;
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
GEM_INLINE void search_stage_decode_candidates_buffer_send(
    search_stage_decode_candidates_buffer_t* const decode_candidates_buffer) {
  gpu_buffer_fmi_decode_send(decode_candidates_buffer->gpu_buffer_fmi_decode);
}
GEM_INLINE void search_stage_decode_candidates_buffer_receive(
    search_stage_decode_candidates_buffer_t* const decode_candidates_buffer) {
  gpu_buffer_fmi_decode_receive(decode_candidates_buffer->gpu_buffer_fmi_decode);
}
/*
 * Accessors
 */
GEM_INLINE void search_stage_decode_candidates_buffer_add(
    search_stage_decode_candidates_buffer_t* const decode_candidates_buffer,
    archive_search_t* const archive_search) {
  // Add archive-search
  vector_insert(decode_candidates_buffer->archive_searches,archive_search,archive_search_t*);
  // Copy the candidate-positions (encoded) to the buffer
  gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode = decode_candidates_buffer->gpu_buffer_fmi_decode;
  archive_search_se_stepwise_decode_candidates_copy(archive_search,gpu_buffer_fmi_decode);
}
GEM_INLINE void search_stage_decode_candidates_buffer_retrieve(
    search_stage_decode_candidates_buffer_t* const decode_candidates_buffer,
    const uint64_t search_idx,archive_search_t** const archive_search) {
  // Retrieve archive-search
  *archive_search = *vector_get_elm(decode_candidates_buffer->archive_searches,search_idx,archive_search_t*);
  // Retrieve candidate-positions (decoded) from the buffer
  gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode = decode_candidates_buffer->gpu_buffer_fmi_decode;
  archive_search_se_stepwise_decode_candidates_retrieve(*archive_search,gpu_buffer_fmi_decode);
}


