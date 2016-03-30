/*
 * PROJECT: GEMMapper
 * FILE: search_state_decode_candidates.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_pipeline/search_stage_decode_candidates.h"
#include "search_pipeline/search_stage_decode_candidates_buffer.h"
#include "archive/archive_search_se_stepwise.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Error Messages
 */
#define GEM_ERROR_SEARCH_STAGE_DC_UNPAIRED_QUERY "Search-Group Stage Buffer. Couldn't retrieve query-pair"

/*
 * Internal Accessors
 */
search_stage_decode_candidates_buffer_t* search_stage_dc_get_buffer(
    search_stage_decode_candidates_t* const search_stage_dc,
    const uint64_t buffer_pos) {
  return *vector_get_elm(search_stage_dc->buffers,buffer_pos,search_stage_decode_candidates_buffer_t*);
}
search_stage_decode_candidates_buffer_t* search_stage_dc_get_current_buffer(
    search_stage_decode_candidates_t* const search_stage_dc) {
  return search_stage_dc_get_buffer(search_stage_dc,search_stage_dc->iterator.current_buffer_idx);
}
/*
 * Setup
 */
search_stage_decode_candidates_t* search_stage_decode_candidates_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const uint64_t num_buffers,
    fm_index_t* const fm_index,
    const bool gpu_decode_sa,
    const bool gpu_decode_text,
    mm_stack_t* const mm_stack) {
  // Alloc
  search_stage_decode_candidates_t* const search_stage_dc = mm_alloc(search_stage_decode_candidates_t);
  // Init Buffers
  uint64_t i;
  search_stage_dc->buffers = vector_new(num_buffers,search_stage_decode_candidates_buffer_t*);
  for (i=0;i<num_buffers;++i) {
    search_stage_decode_candidates_buffer_t* const buffer_vc = search_stage_decode_candidates_buffer_new(
        gpu_buffer_collection,buffers_offset+i,fm_index,gpu_decode_sa,gpu_decode_text);
    vector_insert(search_stage_dc->buffers,buffer_vc,search_stage_decode_candidates_buffer_t*);
  }
  search_stage_dc->iterator.num_buffers = num_buffers;
  search_stage_decode_candidates_clear(search_stage_dc,NULL);
  // Init Support Data Structures
  filtering_candidates_init(&search_stage_dc->filtering_candidates_forward_end1);
  filtering_candidates_init(&search_stage_dc->filtering_candidates_reverse_end1);
  filtering_candidates_init(&search_stage_dc->filtering_candidates_forward_end2);
  filtering_candidates_init(&search_stage_dc->filtering_candidates_reverse_end2);
  search_stage_dc->mm_stack = mm_stack;
  // Return
  return search_stage_dc;
}
void search_stage_decode_candidates_prepare_se_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t* const archive_search) {
  // Clear & Inject Support Data Structures
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_forward_end1);
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_reverse_end1);
  archive_search_inject_filtering_candidates(archive_search,
      &search_stage_dc->filtering_candidates_forward_end1,
      &search_stage_dc->filtering_candidates_reverse_end1,
      NULL, /* Not needed for decoding */
      search_stage_dc->mm_stack);
}
void search_stage_decode_candidates_prepare_pe_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  // Clear & Inject Support Data Structures (End/1)
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_forward_end1);
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_reverse_end1);
  archive_search_inject_filtering_candidates(archive_search_end1,
      &search_stage_dc->filtering_candidates_forward_end1,
      &search_stage_dc->filtering_candidates_reverse_end1,
      NULL, /* Not needed for decoding */
      search_stage_dc->mm_stack);
  // Clear & Inject Support Data Structures (End/2)
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_forward_end2);
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_reverse_end2);
  archive_search_inject_filtering_candidates(archive_search_end2,
      &search_stage_dc->filtering_candidates_forward_end2,
      &search_stage_dc->filtering_candidates_reverse_end2,
      NULL, /* Not needed for decoding */
      search_stage_dc->mm_stack);
}
void search_stage_decode_candidates_clear(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_cache_t* const archive_search_cache) {
  // Init state
  search_stage_dc->search_stage_mode = search_group_buffer_phase_sending;
  // Clear & Init buffers
  const uint64_t num_buffers = search_stage_dc->iterator.num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    search_stage_decode_candidates_buffer_clear(
        search_stage_dc_get_buffer(search_stage_dc,i),archive_search_cache);
  }
  search_stage_dc->iterator.current_buffer_idx = 0; // Init iterator
}
void search_stage_decode_candidates_delete(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_cache_t* const archive_search_cache) {
  // Delete buffers
  const uint64_t num_buffers = search_stage_dc->iterator.num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    search_stage_decode_candidates_buffer_delete(
        search_stage_dc_get_buffer(search_stage_dc,i),archive_search_cache);
  }
  vector_delete(search_stage_dc->buffers); // Delete vector
  // Delete Support Data Structures
  filtering_candidates_destroy(&search_stage_dc->filtering_candidates_forward_end1);
  filtering_candidates_destroy(&search_stage_dc->filtering_candidates_reverse_end1);
  filtering_candidates_destroy(&search_stage_dc->filtering_candidates_forward_end2);
  filtering_candidates_destroy(&search_stage_dc->filtering_candidates_reverse_end2);
  // Free handler
  mm_free(search_stage_dc);
}
/*
 * Send Searches (buffered)
 */
bool search_stage_decode_candidates_send_se_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t* const archive_search) {
  // Check Occupancy (fits in current buffer)
  search_stage_decode_candidates_buffer_t* current_buffer = search_stage_dc_get_current_buffer(search_stage_dc);
  while (!search_stage_decode_candidates_buffer_fits(current_buffer,archive_search,NULL)) {
    // Change group
    const uint64_t last_buffer_idx = search_stage_dc->iterator.num_buffers - 1;
    if (search_stage_dc->iterator.current_buffer_idx < last_buffer_idx) {
      // Send the current group to verification
      search_stage_decode_candidates_buffer_send(current_buffer);
      // Next buffer
      ++(search_stage_dc->iterator.current_buffer_idx);
      current_buffer = search_stage_dc_get_current_buffer(search_stage_dc);
    } else {
      return false;
    }
  }
  // Add SE Search
  search_stage_decode_candidates_buffer_add(current_buffer,archive_search);
  // Copy the candidate-positions (encoded) to the buffer
  archive_search_se_stepwise_decode_candidates_copy(archive_search,current_buffer->gpu_buffer_fmi_decode);
  // Return ok
  return true;
}
bool search_stage_decode_candidates_send_pe_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  // Check Occupancy (fits in current buffer)
  search_stage_decode_candidates_buffer_t* current_buffer = search_stage_dc_get_current_buffer(search_stage_dc);
  while (!search_stage_decode_candidates_buffer_fits(current_buffer,archive_search_end1,archive_search_end2)) {
    // Change group
    const uint64_t last_buffer_idx = search_stage_dc->iterator.num_buffers - 1;
    if (search_stage_dc->iterator.current_buffer_idx < last_buffer_idx) {
      // Send the current group to verification
      search_stage_decode_candidates_buffer_send(current_buffer);
      // Next buffer
      ++(search_stage_dc->iterator.current_buffer_idx);
      current_buffer = search_stage_dc_get_current_buffer(search_stage_dc);
    } else {
      return false;
    }
  }
  // Add PE Search
  search_stage_decode_candidates_buffer_add(current_buffer,archive_search_end1);
  search_stage_decode_candidates_buffer_add(current_buffer,archive_search_end2);
  // Copy the candidate-positions (encoded) to the buffer
  gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode = current_buffer->gpu_buffer_fmi_decode;
  archive_search_se_stepwise_decode_candidates_copy(archive_search_end1,gpu_buffer_fmi_decode);
  archive_search_se_stepwise_decode_candidates_copy(archive_search_end2,gpu_buffer_fmi_decode);
  // Return ok
  return true;
}
/*
 * Retrieve operators
 */
void search_stage_decode_candidates_retrieve_begin(search_stage_decode_candidates_t* const search_stage_dc) {
  search_stage_iterator_t* const iterator = &search_stage_dc->iterator;
  search_stage_decode_candidates_buffer_t* current_buffer;
  // Change mode
  search_stage_dc->search_stage_mode = search_group_buffer_phase_retrieving;
  PROF_ADD_COUNTER(GP_SEARCH_STAGE_DECODE_CANDIDATES_BUFFERS_USED,iterator->current_buffer_idx+1);
  // Send the current buffer
  current_buffer = search_stage_dc_get_current_buffer(search_stage_dc);
  search_stage_decode_candidates_buffer_send(current_buffer);
  // Initialize the iterator
  iterator->current_buffer_idx = 0;
  // Reset searches iterator
  current_buffer = search_stage_dc_get_current_buffer(search_stage_dc);
  iterator->current_search_idx = 0;
  iterator->num_searches = vector_get_used(current_buffer->archive_searches);
  // Fetch first group
  search_stage_decode_candidates_buffer_receive(current_buffer);
}
bool search_stage_decode_candidates_retrieve_finished(
    search_stage_decode_candidates_t* const search_stage_dc) {
  // Mode Sending (Retrieval finished)
  if (search_stage_dc->search_stage_mode==search_group_buffer_phase_sending) return true;
  // Mode Retrieve (Check iterator)
  search_stage_iterator_t* const iterator = &search_stage_dc->iterator;
  return iterator->current_buffer_idx==iterator->num_buffers &&
         iterator->current_search_idx==iterator->num_searches;
}
bool search_stage_decode_candidates_retrieve_next(
    search_stage_decode_candidates_t* const search_stage_dc,
    search_stage_decode_candidates_buffer_t** const current_buffer,
    archive_search_t** const archive_search) {
  // Check state
  if (search_stage_dc->search_stage_mode == search_group_buffer_phase_sending) {
    search_stage_decode_candidates_retrieve_begin(search_stage_dc);
  }
  // Check end-of-iteration
  *current_buffer = search_stage_dc_get_current_buffer(search_stage_dc);
  search_stage_iterator_t* const iterator = &search_stage_dc->iterator;
  if (iterator->current_search_idx==iterator->num_searches) {
    // Next buffer
    ++(iterator->current_buffer_idx);
    if (iterator->current_buffer_idx==iterator->num_buffers) return false;
    // Reset searches iterator
    *current_buffer = search_stage_dc_get_current_buffer(search_stage_dc);
    iterator->current_search_idx = 0;
    iterator->num_searches = vector_get_used((*current_buffer)->archive_searches);
    if (iterator->num_searches==0) return false;
    // Receive Buffer
    search_stage_decode_candidates_buffer_receive(*current_buffer);
  }
  // Retrieve Search
  search_stage_decode_candidates_buffer_retrieve(*current_buffer,iterator->current_search_idx,archive_search);
  ++(iterator->current_search_idx); // Next
  return true;
}
/*
 * Retrieve Searches (buffered)
 */
bool search_stage_decode_candidates_retrieve_se_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t** const archive_search) {
  // Retrieve next
  search_stage_decode_candidates_buffer_t* current_buffer;
  const bool success = search_stage_decode_candidates_retrieve_next(search_stage_dc,&current_buffer,archive_search);
  if (!success) return false;
  // Clear & Inject Support Data Structures
  search_stage_decode_candidates_prepare_se_search(search_stage_dc,*archive_search);
  // Retrieve candidate-positions (decoded) from the buffer
  archive_search_se_stepwise_decode_candidates_retrieve(*archive_search,current_buffer->gpu_buffer_fmi_decode);
  // Return
  return true;
}
bool search_stage_decode_candidates_retrieve_pe_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  search_stage_decode_candidates_buffer_t* current_buffer;
  bool success;
  /*
   * End/1
   */
  // Retrieve next (End/1)
  success = search_stage_decode_candidates_retrieve_next(search_stage_dc,&current_buffer,archive_search_end1);
  if (!success) return false;
  // Clear & Inject Support Data Structures (End/1)
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_forward_end1);
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_reverse_end1);
  archive_search_inject_filtering_candidates(*archive_search_end1,
      &search_stage_dc->filtering_candidates_forward_end1,
      &search_stage_dc->filtering_candidates_reverse_end1,
      NULL, /* Not needed for decoding */
      search_stage_dc->mm_stack);
  // Retrieve candidate-positions (decoded) from the buffer (End/1)
  archive_search_se_stepwise_decode_candidates_retrieve(*archive_search_end1,current_buffer->gpu_buffer_fmi_decode);
  /*
   * End/2
   */
  // Retrieve next (End/2)
  success = search_stage_decode_candidates_retrieve_next(search_stage_dc,&current_buffer,archive_search_end2);
  gem_cond_fatal_error(!success,SEARCH_STAGE_DC_UNPAIRED_QUERY);
  // Clear & Inject Support Data Structures (End/2)
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_forward_end2);
  filtering_candidates_clear(&search_stage_dc->filtering_candidates_reverse_end2);
  archive_search_inject_filtering_candidates(*archive_search_end2,
      &search_stage_dc->filtering_candidates_forward_end2,
      &search_stage_dc->filtering_candidates_reverse_end2,
      NULL, /* Not needed for decoding */
      search_stage_dc->mm_stack);
  // Retrieve candidate-positions (decoded) from the buffer (End/2)
  archive_search_se_stepwise_decode_candidates_retrieve(*archive_search_end2,current_buffer->gpu_buffer_fmi_decode);
  // Return ok
  return true;
}

