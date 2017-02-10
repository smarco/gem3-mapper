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

#include "search_pipeline/search_stage_verify_candidates_buffer.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Setup
 */
search_stage_verify_candidates_buffer_t* search_stage_verify_candidates_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool verify_candidates_enabled) {
  // Alloc
  search_stage_verify_candidates_buffer_t* const verify_candidates_buffer = mm_alloc(search_stage_verify_candidates_buffer_t);
  // Init
  verify_candidates_buffer->gpu_buffer_align_bpm =
      gpu_buffer_align_bpm_new(gpu_buffer_collection,buffer_no,verify_candidates_enabled);
  const uint64_t max_queries = gpu_buffer_align_bpm_get_max_queries(verify_candidates_buffer->gpu_buffer_align_bpm);
  verify_candidates_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return verify_candidates_buffer;
}
void search_stage_verify_candidates_buffer_clear(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    archive_search_cache_t* const archive_search_cache) {
  gpu_buffer_align_bpm_clear(verify_candidates_buffer->gpu_buffer_align_bpm);
  // Return searches to the cache
  if (archive_search_cache!=NULL) {
    VECTOR_ITERATE(verify_candidates_buffer->archive_searches,archive_search,n,archive_search_t*) {
      archive_search_cache_free(archive_search_cache,*archive_search);
    }
  }
  // Clear searches vector
  vector_clear(verify_candidates_buffer->archive_searches);
}
void search_stage_verify_candidates_buffer_delete(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    archive_search_cache_t* const archive_search_cache) {
  gpu_buffer_align_bpm_delete(verify_candidates_buffer->gpu_buffer_align_bpm);
  // Return searches to the cache
  VECTOR_ITERATE(verify_candidates_buffer->archive_searches,archive_search,n,archive_search_t*) {
    archive_search_cache_free(archive_search_cache,*archive_search);
  }
  // Delete searches vector
  vector_delete(verify_candidates_buffer->archive_searches);
  mm_free(verify_candidates_buffer);
}
/*
 * Occupancy
 */
bool search_stage_verify_candidates_buffer_fits(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm = verify_candidates_buffer->gpu_buffer_align_bpm;
  // Compute dimensions
  uint64_t total_entries = 0,total_queries = 0,total_candidates = 0;
  pattern_t* const pattern_end1 = &archive_search_end1->approximate_search.pattern;
  gpu_buffer_align_bpm_compute_dimensions(gpu_buffer_align_bpm,
      pattern_end1->key_length,pattern_end1->pattern_tiled.num_tiles,
      archive_search_get_num_verify_candidates(archive_search_end1),
      &total_entries,&total_queries,&total_candidates);
  if (archive_search_end2!=NULL) {
    pattern_t* const pattern_end2 = &archive_search_end2->approximate_search.pattern;
    gpu_buffer_align_bpm_compute_dimensions(gpu_buffer_align_bpm,
        pattern_end2->key_length,pattern_end2->pattern_tiled.num_tiles,
        archive_search_get_num_verify_candidates(archive_search_end2),
        &total_entries,&total_queries,&total_candidates);
  }
  // Return if current search fits in buffer
  return gpu_buffer_align_bpm_fits_in_buffer(gpu_buffer_align_bpm,
      total_entries,total_queries,total_candidates);
}
/*
 * Send/Receive
 */
void search_stage_verify_candidates_buffer_send(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer) {
  PROF_ADD_COUNTER(GP_SEARCH_STAGE_VERIFY_CANDIDATES_SEARCHES_IN_BUFFER,
      vector_get_used(verify_candidates_buffer->archive_searches));
  gpu_buffer_align_bpm_send(verify_candidates_buffer->gpu_buffer_align_bpm);
}
void search_stage_verify_candidates_buffer_receive(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer) {
  gpu_buffer_align_bpm_receive(verify_candidates_buffer->gpu_buffer_align_bpm);
}
/*
 * Accessors
 */
void search_stage_verify_candidates_buffer_add(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    archive_search_t* const archive_search) {
  // Add archive-search
  vector_insert(verify_candidates_buffer->archive_searches,archive_search,archive_search_t*);
}
void search_stage_verify_candidates_buffer_retrieve(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search) {
  // Retrieve archive-search
  *archive_search = *vector_get_elm(verify_candidates_buffer->archive_searches,search_idx,archive_search_t*);
}

