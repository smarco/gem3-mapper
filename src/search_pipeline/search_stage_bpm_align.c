/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_pipeline/search_stage_bpm_align.h"
#include "search_pipeline/search_stage_bpm_align_buffer.h"
#include "archive/search/archive_search_se_stepwise.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Error Messages
 */
#define GEM_ERROR_SEARCH_STAGE_BPM_ALIGN_UNPAIRED_QUERY \
  "Search-stage BPM-Align Buffer. Couldn't retrieve query-pair"

/*
 * Internal Accessors
 */
search_stage_bpm_align_buffer_t* search_stage_bpm_align_get_buffer(
    search_stage_bpm_align_t* const search_stage,
    const uint64_t buffer_pos) {
  return *vector_get_elm(search_stage->buffers,buffer_pos,search_stage_bpm_align_buffer_t*);
}
search_stage_bpm_align_buffer_t* search_stage_bpm_align_get_current_buffer(
    search_stage_bpm_align_t* const search_stage) {
  return search_stage_bpm_align_get_buffer(search_stage,search_stage->iterator.current_buffer_idx);
}
/*
 * Setup
 */
search_stage_bpm_align_t* search_stage_bpm_align_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const uint64_t num_buffers,
    const bool paired_end,
    const bool bpm_align_enabled,
    search_pipeline_handlers_t* const search_pipeline_handlers) {
  // Alloc
  search_stage_bpm_align_t* const search_stage = mm_alloc(search_stage_bpm_align_t);
  search_stage->paired_end = paired_end;
  // Init Support Data Structures
  search_stage->search_pipeline_handlers = search_pipeline_handlers;
  if (paired_end) {
    search_stage->paired_matches = paired_matches_new();
  } else {
    search_stage->matches = matches_new();
  }
  // Init Buffers
  uint64_t i;
  search_stage->buffers = vector_new(num_buffers,search_stage_bpm_align_buffer_t*);
  for (i=0;i<num_buffers;++i) {
    search_stage_bpm_align_buffer_t* const bpm_align_buffer =
        search_stage_bpm_align_buffer_new(gpu_buffer_collection,buffers_offset+i,bpm_align_enabled);
    vector_insert(search_stage->buffers,bpm_align_buffer,search_stage_bpm_align_buffer_t*);
  }
  search_stage->iterator.num_buffers = num_buffers;
  search_stage_bpm_align_clear(search_stage); // Clear buffers
  // Return
  return search_stage;
}
void search_stage_bpm_align_clear(
    search_stage_bpm_align_t* const search_stage) {
  // Init state
  search_stage->search_stage_mode = search_group_buffer_phase_sending;
  // Clear & Init buffers
  const uint64_t num_buffers = search_stage->iterator.num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    search_stage_bpm_align_buffer_clear(
        search_stage_bpm_align_get_buffer(search_stage,i));
  }
  search_stage->iterator.current_buffer_idx = 0; // Init iterator
}
void search_stage_bpm_align_delete(
    search_stage_bpm_align_t* const search_stage) {
  // Delete buffers
  const uint64_t num_buffers = search_stage->iterator.num_buffers;
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    search_stage_bpm_align_buffer_delete(
        search_stage_bpm_align_get_buffer(search_stage,i));
  }
  vector_delete(search_stage->buffers); // Delete vector
  // Delete Support Data Structures
  if (search_stage->paired_end) {
    paired_matches_delete(search_stage->paired_matches);
  } else {
    matches_delete(search_stage->matches);
  }
  // Free handler
  mm_free(search_stage);
}
/*
 * Send Searches (buffered)
 */
bool search_stage_bpm_align_send_se_search(
    search_stage_bpm_align_t* const search_stage,
    archive_search_t* const archive_search) {
  // Check Occupancy (fits in current buffer)
  search_stage_bpm_align_buffer_t* current_buffer =
      search_stage_bpm_align_get_current_buffer(search_stage);
  while (!search_stage_bpm_align_buffer_fits(current_buffer,archive_search,NULL)) {
    // Change group
    const uint64_t last_buffer_idx = search_stage->iterator.num_buffers - 1;
    if (search_stage->iterator.current_buffer_idx < last_buffer_idx) {
      // Send the current group to align
      search_stage_bpm_align_buffer_send(current_buffer);
      // Next buffer
      ++(search_stage->iterator.current_buffer_idx);
      current_buffer = search_stage_bpm_align_get_current_buffer(search_stage);
    } else {
      return false;
    }
  }
  // Add SE Search
  search_stage_bpm_align_buffer_add(current_buffer,archive_search);
  // Copy candidates to the buffer
  archive_search_se_stepwise_bpm_align_copy(archive_search,current_buffer->gpu_buffer_bpm_align);
  // Return ok
  return true;
}
bool search_stage_bpm_align_send_pe_search(
    search_stage_bpm_align_t* const search_stage,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  // Check Occupancy (fits in current buffer)
  search_stage_bpm_align_buffer_t* current_buffer =
      search_stage_bpm_align_get_current_buffer(search_stage);
  while (!search_stage_bpm_align_buffer_fits(current_buffer,archive_search_end1,archive_search_end2)) {
    // Change group
    const uint64_t last_buffer_idx = search_stage->iterator.num_buffers - 1;
    if (search_stage->iterator.current_buffer_idx < last_buffer_idx) {
      // Send the current group to verification
      search_stage_bpm_align_buffer_send(current_buffer);
      // Next buffer
      ++(search_stage->iterator.current_buffer_idx);
      current_buffer = search_stage_bpm_align_get_current_buffer(search_stage);
    } else {
      return false;
    }
  }
  // Add PE Search
  search_stage_bpm_align_buffer_add(current_buffer,archive_search_end1);
  search_stage_bpm_align_buffer_add(current_buffer,archive_search_end2);
  // Copy candidates to the buffer
  archive_search_se_stepwise_bpm_align_copy(archive_search_end1,current_buffer->gpu_buffer_bpm_align);
  archive_search_se_stepwise_bpm_align_copy(archive_search_end2,current_buffer->gpu_buffer_bpm_align);
  // Return ok
  return true;
}
/*
 * Retrieve operators
 */
void search_stage_bpm_align_retrieve_begin(search_stage_bpm_align_t* const search_stage) {
  search_stage_iterator_t* const iterator = &search_stage->iterator;
  search_stage_bpm_align_buffer_t* current_buffer;
  // Change mode
  search_stage->search_stage_mode = search_group_buffer_phase_retrieving;
  PROF_ADD_COUNTER(GP_SSTAGE_BPM_ALIGN_BUFFERS,iterator->current_buffer_idx+1);
  // Send the current buffer
  current_buffer = search_stage_bpm_align_get_current_buffer(search_stage);
  search_stage_bpm_align_buffer_send(current_buffer);
  // Initialize the iterator
  iterator->current_buffer_idx = 0;
  // Reset searches iterator
  current_buffer = search_stage_bpm_align_get_current_buffer(search_stage);
  iterator->current_search_idx = 0;
  iterator->num_searches = vector_get_used(current_buffer->archive_searches);
  // Fetch first group
  search_stage_bpm_align_buffer_receive(current_buffer);
}
bool search_stage_bpm_align_retrieve_finished(
    search_stage_bpm_align_t* const search_stage) {
  // Mode Sending (Retrieval finished)
  if (search_stage->search_stage_mode==search_group_buffer_phase_sending) return true;
  // Mode Retrieve (Check iterator)
  search_stage_iterator_t* const iterator = &search_stage->iterator;
  return iterator->current_buffer_idx==iterator->num_buffers &&
         iterator->current_search_idx==iterator->num_searches;
}
bool search_stage_bpm_align_retrieve_next(
    search_stage_bpm_align_t* const search_stage,
    search_stage_bpm_align_buffer_t** const current_buffer,
    archive_search_t** const archive_search) {
  // Check state
  if (search_stage->search_stage_mode == search_group_buffer_phase_sending) {
    search_stage_bpm_align_retrieve_begin(search_stage);
  }
  // Check end-of-iteration
  *current_buffer = search_stage_bpm_align_get_current_buffer(search_stage);
  search_stage_iterator_t* const iterator = &search_stage->iterator;
  if (iterator->current_search_idx==iterator->num_searches) {
    // Next buffer
    ++(iterator->current_buffer_idx);
    if (iterator->current_buffer_idx==iterator->num_buffers) return false;
    // Reset searches iterator
    *current_buffer = search_stage_bpm_align_get_current_buffer(search_stage);
    iterator->current_search_idx = 0;
    iterator->num_searches = vector_get_used((*current_buffer)->archive_searches);
    if (iterator->num_searches==0) return false;
    // Receive Buffer
    search_stage_bpm_align_buffer_receive(*current_buffer);
  }
  // Retrieve Search
  search_stage_bpm_align_buffer_retrieve(*current_buffer,iterator->current_search_idx,archive_search);
  ++(iterator->current_search_idx); // Next
  return true;
}
/*
 * Retrieve Searches (buffered)
 */
void search_stage_bpm_align_prepare(
    archive_search_t* const archive_search,
    search_pipeline_handlers_t* const search_pipeline_handlers,
    filtering_candidates_t* const filtering_candidates) {
  // Parameters
  approximate_search_t* const search = &archive_search->approximate_search;
  // Prepare Support Data Structures
  filtering_candidates_clear(filtering_candidates,true);
  search->filtering_candidates = filtering_candidates;
  filtering_candidates_inject_handlers(
      filtering_candidates,archive_search->archive,
      &archive_search->search_parameters,
      search_pipeline_handlers->mm_allocator);
}
bool search_stage_bpm_align_retrieve_se_search(
    search_stage_bpm_align_t* const search_stage,
    archive_search_t** const archive_search) {
  // Retrieve next
  search_stage_bpm_align_buffer_t* current_buffer;
  const bool success = search_stage_bpm_align_retrieve_next(search_stage,&current_buffer,archive_search);
  if (!success) return false;
  // Prepare Archive Search
  search_pipeline_handlers_t* const search_pipeline_handlers = search_stage->search_pipeline_handlers;
  matches_clear(search_stage->matches);
  search_stage_bpm_align_prepare(
      *archive_search,search_pipeline_handlers,
      &search_pipeline_handlers->fc_bpm_align_end1);
  // Retrieve candidates from the buffer
  archive_search_se_stepwise_bpm_align_retrieve(*archive_search,
      current_buffer->gpu_buffer_bpm_align,search_stage->matches);
  // Return
  return true;
}
bool search_stage_bpm_align_retrieve_pe_search(
    search_stage_bpm_align_t* const search_stage,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  search_stage_bpm_align_buffer_t* current_buffer;
  bool success;
  /*
   * End/1
   */
  // Retrieve next (End/1)
  success = search_stage_bpm_align_retrieve_next(search_stage,&current_buffer,archive_search_end1);
  if (!success) return false;
  // Clear paired-matches
  paired_matches_clear(search_stage->paired_matches,true);
  // Prepare Archive Search
  search_pipeline_handlers_t* const search_pipeline_handlers = search_stage->search_pipeline_handlers;
  search_stage_bpm_align_prepare(
      *archive_search_end1,search_pipeline_handlers,
      &search_pipeline_handlers->fc_bpm_align_end1);
  // Retrieve candidates from the buffer (End/1)
  archive_search_se_stepwise_bpm_align_retrieve(
      *archive_search_end1,current_buffer->gpu_buffer_bpm_align,
      search_stage->paired_matches->matches_end1);
  /*
   * End/2
   */
  // Retrieve next (End/2)
  success = search_stage_bpm_align_retrieve_next(search_stage,&current_buffer,archive_search_end2);
  gem_cond_fatal_error(!success,SEARCH_STAGE_BPM_ALIGN_UNPAIRED_QUERY);
  // Prepare Archive Search
  search_stage_bpm_align_prepare(
      *archive_search_end2,search_pipeline_handlers,
      &search_pipeline_handlers->fc_bpm_align_end2);
  // Retrieve candidates from the buffer (End/2)
  archive_search_se_stepwise_bpm_align_retrieve(
      *archive_search_end2,current_buffer->gpu_buffer_bpm_align,
      search_stage->paired_matches->matches_end2);
  // Return ok
  return true;
}
