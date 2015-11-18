/*
 * PROJECT: GEMMapper
 * FILE: search_group_verify_candidates.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_stage_verify_candidates.h"
#include "search_stage_verify_candidates_buffer.h"
#include "archive_search_se_stepwise.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Error Messages
 */
#define GEM_ERROR_SEARCH_STAGE_VC_UNPAIRED_QUERY "Search-Group Stage Buffer. Couldn't retrieve query-pair"

/*
 * Internal Accessors
 */
GEM_INLINE search_stage_verify_candidates_buffer_t* search_stage_vc_get_buffer(
    search_stage_verify_candidates_t* const search_stage_vc,const uint64_t buffer_pos) {
  return *vector_get_elm(search_stage_vc->buffers,buffer_pos,search_stage_verify_candidates_buffer_t*);
}
GEM_INLINE search_stage_verify_candidates_buffer_t* search_stage_vc_get_current_buffer(
    search_stage_verify_candidates_t* const search_stage_vc) {
  return search_stage_vc_get_buffer(search_stage_vc,search_stage_vc->iterator.current_buffer_idx);
}
/*
 * Setup
 */
GEM_INLINE search_stage_verify_candidates_t* search_stage_verify_candidates_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,const uint64_t num_buffers,
    const bool cpu_emulated) {
  // Alloc
  search_stage_verify_candidates_t* const search_stage_vc = mm_alloc(search_stage_verify_candidates_t);
  // Init Buffers
  uint64_t i;
  search_stage_vc->buffers = vector_new(num_buffers,search_stage_verify_candidates_buffer_t*);
  for (i=0;i<num_buffers;++i) {
    search_stage_verify_candidates_buffer_t* const buffer_vc =
        search_stage_verify_candidates_buffer_new(gpu_buffer_collection,buffers_offset+i,cpu_emulated);
    vector_insert(search_stage_vc->buffers,buffer_vc,search_stage_verify_candidates_buffer_t*);
  }
  // Return
  return search_stage_vc;
}
GEM_INLINE void search_stage_verify_candidates_clear(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_cache_t* const archive_search_cache) {
  // Init state
  search_stage_vc->search_stage_mode = search_group_buffer_phase_sending;
  // Clear & Init buffers
  const uint64_t num_buffers = search_stage_vc->iterator.num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    search_stage_verify_candidates_buffer_clear(
        search_stage_vc_get_buffer(search_stage_vc,i),archive_search_cache);
  }
  search_stage_vc->iterator.current_buffer_idx = 0; // Init iterator
}
GEM_INLINE void search_stage_verify_candidates_delete(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_cache_t* const archive_search_cache) {
  // Delete buffers
  const uint64_t num_buffers = search_stage_vc->iterator.num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    search_stage_verify_candidates_buffer_delete(
        search_stage_vc_get_buffer(search_stage_vc,i),archive_search_cache);
  }
  vector_delete(search_stage_vc->buffers); // Delete vector
  mm_free(search_stage_vc); // Free handler
}
/*
 * Accessors
 */
GEM_INLINE bool search_stage_verify_candidates_is_empty(search_stage_verify_candidates_t* const search_stage_vc) {
  // TODO TODO TODO TODO TODO TODO TODO
  return false;
}
/*
 * Send Searches (buffered)
 */
GEM_INLINE bool search_stage_verify_candidates_send_se_search(
    search_stage_verify_candidates_t* const search_stage_vc,archive_search_t* const archive_search) {
  // Check Occupancy (fits in current buffer)
  search_stage_verify_candidates_buffer_t* current_buffer = search_stage_vc_get_current_buffer(search_stage_vc);
  while (!search_stage_verify_candidates_buffer_fits(current_buffer,archive_search,NULL)) {
    // Change group
    const uint64_t last_buffer_idx = search_stage_vc->iterator.num_buffers - 1;
    if (search_stage_vc->iterator.current_buffer_idx < last_buffer_idx) {
      // Send the current group to verification
      search_stage_verify_candidates_buffer_send(current_buffer);
      // Next buffer
      ++(search_stage_vc->iterator.current_buffer_idx);
      current_buffer = search_stage_vc_get_current_buffer(search_stage_vc);
    } else {
      return false;
    }
  }
  // Add SE Search
  search_stage_verify_candidates_buffer_add(current_buffer,archive_search);
  // Return ok
  return true;
}
GEM_INLINE bool search_stage_verify_candidates_send_pe_search(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
  // Check Occupancy (fits in current buffer)
  search_stage_verify_candidates_buffer_t* current_buffer = search_stage_vc_get_current_buffer(search_stage_vc);
  while (!search_stage_verify_candidates_buffer_fits(current_buffer,archive_search_end1,archive_search_end2)) {
    // Change group
    const uint64_t last_buffer_idx = search_stage_vc->iterator.num_buffers - 1;
    if (search_stage_vc->iterator.current_buffer_idx < last_buffer_idx) {
      // Send the current group to verification
      search_stage_verify_candidates_buffer_send(current_buffer);
      // Next buffer
      ++(search_stage_vc->iterator.current_buffer_idx);
      current_buffer = search_stage_vc_get_current_buffer(search_stage_vc);
    } else {
      return false;
    }
  }
  // Add PE Search
  search_stage_verify_candidates_buffer_add(current_buffer,archive_search_end1);
  search_stage_verify_candidates_buffer_add(current_buffer,archive_search_end2);
  // Return ok
  return true;
}
/*
 * Retrieve operators
 */
GEM_INLINE void search_stage_verify_candidates_retrieve_begin(search_stage_verify_candidates_t* const search_stage_vc) {
  // PROFILE
  // FIXME PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_GROUP_BUFFERS_USED,search_group_vc->current_buffer_idx);
  search_stage_verify_candidates_buffer_t* current_buffer;
  // Send the current buffer
  current_buffer = search_stage_vc_get_current_buffer(search_stage_vc);
  search_stage_verify_candidates_buffer_send(current_buffer);
  // Initialize the iterator
  search_stage_iterator_t* const iterator = &search_stage_vc->iterator;
  iterator->current_buffer_idx = 0;
  // Reset searches iterator
  current_buffer = search_stage_vc_get_current_buffer(search_stage_vc);
  iterator->current_search_idx = 0;
  iterator->num_searches = vector_get_used(current_buffer->archive_searches);
  // Fetch first group
  search_stage_verify_candidates_buffer_receive(current_buffer);
}
GEM_INLINE bool search_stage_verify_candidates_retrieve_next(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_t** const archive_search,matches_t* const matches) {
  // Check state
  if (search_stage_vc->search_stage_mode == search_group_buffer_phase_sending) {
    search_stage_verify_candidates_retrieve_begin(search_stage_vc);
    search_stage_vc->search_stage_mode = search_group_buffer_phase_retrieving;
  }
  // Check end-of-iteration
  search_stage_verify_candidates_buffer_t* current_buffer = search_stage_vc_get_current_buffer(search_stage_vc);
  search_stage_iterator_t* const iterator = &search_stage_vc->iterator;
  if (iterator->current_search_idx==iterator->num_searches) {
    // Next buffer
    ++(iterator->current_buffer_idx);
    if (iterator->current_buffer_idx==iterator->num_buffers) return false;
    // Reset searches iterator
    current_buffer = search_stage_vc_get_current_buffer(search_stage_vc);
    iterator->current_search_idx = 0;
    iterator->num_searches = vector_get_used(current_buffer->archive_searches);
    if (iterator->num_searches==0) return false;
    // Receive Buffer
    search_stage_verify_candidates_buffer_receive(current_buffer);
  }
  // Retrieve Search
  search_stage_verify_candidates_buffer_retrieve(current_buffer,iterator->current_search_idx,archive_search,matches);
  ++(iterator->current_search_idx); // Next
  return true;
}
/*
 * Retrieve Searches (buffered)
 */
GEM_INLINE bool search_stage_verify_candidates_retrieve_se_search(
    search_stage_verify_candidates_t* const search_stage_vc,archive_search_t** const archive_search,
    text_collection_t* const text_collection,matches_t* const matches) {
  matches_clear(matches); // Clear Matches
  text_collection_clear(text_collection); // Clear text-collection
  return search_stage_verify_candidates_retrieve_next(search_stage_vc,archive_search,matches);
}
GEM_INLINE bool search_stage_verify_candidates_retrieve_pe_search(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_t** const archive_search_end1,archive_search_t** const archive_search_end2,
    text_collection_t* const text_collection,paired_matches_t* const paired_matches) {
  // Init
  bool success;
  paired_matches_clear(paired_matches);
  text_collection_clear(text_collection); // Clear text-collection
  // Get End/1
  success = search_stage_verify_candidates_retrieve_next(
      search_stage_vc,archive_search_end1,paired_matches->matches_end1);
  if (!success) return false;
  // Get End/2
  success = search_stage_verify_candidates_retrieve_next(
      search_stage_vc,archive_search_end2,paired_matches->matches_end2);
  gem_cond_fatal_error(!success,SEARCH_STAGE_VC_UNPAIRED_QUERY);
  // Return ok
  return true;
}
