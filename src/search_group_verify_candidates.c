/*
 * PROJECT: GEMMapper
 * FILE: search_group_verify_candidates.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_group_verify_candidates.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Setup
 */
GEM_INLINE search_group_verify_candidates_t* search_group_verify_candidates_new(
    bpm_gpu_buffer_t* const bpm_gpu_buffers,const uint64_t num_buffers,
    const bool cpu_emulated) {
  // Alloc
  search_group_verify_candidates_t* const search_group_vc = mm_alloc(search_group_verify_candidates_t);
  // Configure
  search_group_vc->hint_patterns_per_search = 1;
  search_group_vc->cpu_emulated = cpu_emulated;
  // Buffers
  search_group_vc->num_buffers = num_buffers;
  const uint64_t num_initial_searches = bpm_gpu_buffer_get_max_queries(bpm_gpu_buffers);
  search_group_vc->buffers = mm_calloc(num_buffers,search_group_verify_candidates_buffer_t,true);
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    search_group_vc->buffers[i].archive_searches = vector_new(num_initial_searches,archive_search_t*);
    search_group_vc->buffers[i].bpm_gpu_buffer = bpm_gpu_buffers + i;
  }
  search_group_vc->current_buffer_idx = 0;
  search_group_vc->current_buffer = search_group_vc->buffers;
  // Return
  return search_group_vc;
}
GEM_INLINE void search_group_verify_candidates_init(search_group_verify_candidates_t* const search_group_vc) {
  // Initialize BPM-Buffers
  const uint64_t num_buffers = search_group_vc->num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    bpm_gpu_init_buffer(search_group_vc->buffers[i].bpm_gpu_buffer);
  }
}
GEM_INLINE void search_group_verify_candidates_clear(
    search_group_verify_candidates_t* const search_group_vc,
    archive_search_cache_t* const archive_search_cache) {
  // Detach all archive searches & clear BPM-buffers
  const uint64_t num_buffers = search_group_vc->num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    // Clear archive-searches
    vector_t* const archive_searches = search_group_vc->buffers[i].archive_searches;
    if (archive_search_cache!=NULL) {
      VECTOR_ITERATE(archive_searches,archive_search,j,archive_search_t*) {
        archive_search_cache_free(archive_search_cache,*archive_search);
      }
    }
    vector_clear(archive_searches);
    // Clear BPM-buffer
    bpm_gpu_buffer_clear(search_group_vc->buffers[i].bpm_gpu_buffer);
  }
  search_group_vc->current_buffer_idx = 0;
  search_group_vc->current_buffer = search_group_vc->buffers;
}
GEM_INLINE void search_group_verify_candidates_delete(search_group_verify_candidates_t* const search_group_vc) {
  // Free all
  const uint64_t num_buffers = search_group_vc->num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    vector_delete(search_group_vc->buffers[i].archive_searches);
  }
  mm_free(search_group_vc->buffers);
  // Free handler
  mm_free(search_group_vc);
}
/*
 * Accessors
 */
GEM_INLINE bool search_group_verify_candidates_buffer_is_empty(search_group_verify_candidates_t* const search_group_vc) {
  return vector_is_empty(search_group_vc->current_buffer->archive_searches);
}
GEM_INLINE bool search_group_verify_candidates_fits_in_buffer(
    search_group_verify_candidates_t* const search_group_vc,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
  bpm_gpu_buffer_t* const bpm_gpu_buffer = search_group_vc->current_buffer->bpm_gpu_buffer;
  // Compute dimensions
  uint64_t total_entries = 0,total_query_chunks = 0,total_candidate_chunks = 0;
  bpm_gpu_buffer_compute_dimensions(bpm_gpu_buffer,
      &archive_search_end1->forward_search_state.pattern,
      archive_search_get_search_canditates(archive_search_end1),
      &total_entries,&total_query_chunks,&total_candidate_chunks);
  if (archive_search_end2!=NULL) {
    bpm_gpu_buffer_compute_dimensions(bpm_gpu_buffer,
        &archive_search_end2->forward_search_state.pattern,
        archive_search_get_search_canditates(archive_search_end2),
        &total_entries,&total_query_chunks,&total_candidate_chunks);
  }
  // Check if current search fits in buffer
  return bpm_gpu_buffer_fits_in_buffer(bpm_gpu_buffer,total_entries,total_query_chunks,total_candidate_chunks);
}
GEM_INLINE void search_group_verify_candidates_retrieve_begin(search_group_verify_candidates_t* const search_group_vc) {
  // Send the current group to verification
  PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_GROUP_BUFFERS_USED,search_group_vc->current_buffer_idx);
  bpm_gpu_buffer_send(search_group_vc->current_buffer->bpm_gpu_buffer);
  // Initialize the iterator
  search_group_vc->current_buffer_idx = 0;
  search_group_vc->current_buffer = search_group_vc->buffers;
  search_group_vc->current_archive_search_idx = 0;
  search_group_vc->current_archive_search = vector_get_mem(search_group_vc->current_buffer->archive_searches,archive_search_t*);
  search_group_vc->num_archive_searches = vector_get_used(search_group_vc->current_buffer->archive_searches);
  // Fetch first group
  if (!search_group_vc->cpu_emulated && search_group_vc->num_archive_searches > 0) {
    PROFILE_START(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY,PROFILE_LEVEL);
    bpm_gpu_buffer_receive(search_group_vc->current_buffer->bpm_gpu_buffer);
    PROFILE_STOP(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY,PROFILE_LEVEL);
  }
}
/*
 * SE search-group verification
 */
GEM_INLINE bool search_group_verify_candidates_se_add_search(
    search_group_verify_candidates_t* const search_group_vc,archive_search_t* const archive_search) {
  // Check if it fits in current buffer
  while (!search_group_verify_candidates_fits_in_buffer(search_group_vc,archive_search,NULL)) {
    const bool empty_group = vector_is_empty(search_group_vc->current_buffer->archive_searches);
    gem_cond_fatal_error(empty_group,SEARCH_GROUP_VERIFY_CANDIDATES_QUERY_TOO_BIG);
    // Change group
    if (search_group_vc->current_buffer_idx < search_group_vc->num_buffers - 1) {
      // Send the current group to verification
      if (!search_group_vc->cpu_emulated) {
        bpm_gpu_buffer_send(search_group_vc->current_buffer->bpm_gpu_buffer);
      }
      // Next group
      ++(search_group_vc->current_buffer_idx);
      ++(search_group_vc->current_buffer);
    } else {
      return false;
    }
  }
  // Copy the candidates to the buffer
  if (!search_group_vc->cpu_emulated) {
    archive_search_se_stepwise_copy_candidates(archive_search,search_group_vc->current_buffer->bpm_gpu_buffer);
  }
  // Add the search to the current group
  vector_insert(search_group_vc->current_buffer->archive_searches,archive_search,archive_search_t*);
  // Return ok
  return true;
}
GEM_INLINE bool search_group_verify_candidates_retrieve_search(
    search_group_verify_candidates_t* const search_group_vc,
    archive_search_t** const archive_search,matches_t* const matches) {
  // Check end-of-iteration
  if (search_group_vc->current_archive_search_idx==search_group_vc->num_archive_searches) {
    // Next group
    ++(search_group_vc->current_buffer_idx);
    ++(search_group_vc->current_buffer);
    if (search_group_vc->current_buffer_idx==search_group_vc->num_buffers) return false;
    // Clear archive-search iterator
    search_group_vc->current_archive_search_idx = 0;
    search_group_vc->current_archive_search = vector_get_mem(search_group_vc->current_buffer->archive_searches,archive_search_t*);
    search_group_vc->num_archive_searches = vector_get_used(search_group_vc->current_buffer->archive_searches);
    if (search_group_vc->num_archive_searches==0) return false;
    // Fetch group
    if (!search_group_vc->cpu_emulated) {
      PROFILE_START(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY,PROFILE_LEVEL);
      bpm_gpu_buffer_receive(search_group_vc->current_buffer->bpm_gpu_buffer);
      PROFILE_STOP(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY,PROFILE_LEVEL);
    }
  }
  // Return next archive-search
  bpm_gpu_buffer_t* const bpm_gpu_buffer = search_group_vc->current_buffer->bpm_gpu_buffer;
  *archive_search = *search_group_vc->current_archive_search;
  // Retrieve candidates
  if (search_group_vc->cpu_emulated) {
    archive_search_se_stepwise_verify_candidates(*archive_search,matches); // CPU emulation
  } else {
    archive_search_se_stepwise_retrieve_candidates(*archive_search,bpm_gpu_buffer,matches);
  }
  // Next archive-search
  ++(search_group_vc->current_archive_search_idx);
  ++(search_group_vc->current_archive_search);
  // Return ok
  return true;
}
GEM_INLINE bool search_group_verify_candidates_se_get_search(
    search_group_verify_candidates_t* const search_group_vc,archive_search_t** const archive_search,
    text_collection_t* const text_collection,matches_t* const matches) {
  matches_clear(matches); // Clear Matches
  text_collection_clear(text_collection); // Clear text-collection
  return search_group_verify_candidates_retrieve_search(search_group_vc,archive_search,matches);
}
/*
 * PE search-group verification
 */
GEM_INLINE bool search_group_verify_candidates_pe_add_search(
    search_group_verify_candidates_t* const search_group_vc,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
  // Check if it fits in current buffer
  while (!search_group_verify_candidates_fits_in_buffer(search_group_vc,archive_search_end1,archive_search_end2)) {
    const bool empty_group = vector_is_empty(search_group_vc->current_buffer->archive_searches);
    gem_cond_fatal_error(empty_group,SEARCH_GROUP_VERIFY_CANDIDATES_QUERY_TOO_BIG);
    // Change group
    if (search_group_vc->current_buffer_idx < search_group_vc->num_buffers - 1) {
      // Send the current group to verification
      if (!search_group_vc->cpu_emulated) {
        bpm_gpu_buffer_send(search_group_vc->current_buffer->bpm_gpu_buffer);
      }
      // Next group
      ++(search_group_vc->current_buffer_idx);
      ++(search_group_vc->current_buffer);
    } else {
      return false;
    }
  }
  // Copy the candidates to the buffer
  if (!search_group_vc->cpu_emulated) {
    archive_search_se_stepwise_copy_candidates(archive_search_end1,search_group_vc->current_buffer->bpm_gpu_buffer);
    archive_search_se_stepwise_copy_candidates(archive_search_end2,search_group_vc->current_buffer->bpm_gpu_buffer);
  }
  // Add the search to the current group
  vector_insert(search_group_vc->current_buffer->archive_searches,archive_search_end1,archive_search_t*);
  vector_insert(search_group_vc->current_buffer->archive_searches,archive_search_end2,archive_search_t*);
  // Return ok
  return true;
}
GEM_INLINE bool search_group_verify_candidates_pe_get_search(
    search_group_verify_candidates_t* const search_group_vc,archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2,text_collection_t* const text_collection,
    paired_matches_t* const paired_matches) {
  // Init
  bool success;
  paired_matches_clear(paired_matches);
  text_collection_clear(text_collection); // Clear text-collection
  // Get End/1
  success = search_group_verify_candidates_retrieve_search(
      search_group_vc,archive_search_end1,paired_matches->matches_end1);
  if (!success) return false;
  // Get End/2
  success = search_group_verify_candidates_retrieve_search(
      search_group_vc,archive_search_end2,paired_matches->matches_end2);
  gem_cond_fatal_error(!success,SEARCH_GROUP_VERIFY_CANDIDATES_UNPAIRED_QUERY);
  // Return ok
  return true;
}

