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

#include "search_pipeline/search_stage_bpm_distance_buffer.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Setup
 */
search_stage_bpm_distance_buffer_t* search_stage_bpm_distance_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool bpm_distance_enabled) {
  // Alloc
  search_stage_bpm_distance_buffer_t* const bpm_distance_buffer =
      mm_alloc(search_stage_bpm_distance_buffer_t);
  // Init
  bpm_distance_buffer->gpu_buffer_bpm_distance =
      gpu_buffer_bpm_distance_new(gpu_buffer_collection,buffer_no,bpm_distance_enabled);
  const uint64_t max_queries =
      gpu_buffer_bpm_distance_get_max_queries(bpm_distance_buffer->gpu_buffer_bpm_distance);
  bpm_distance_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return bpm_distance_buffer;
}
void search_stage_bpm_distance_buffer_clear(
    search_stage_bpm_distance_buffer_t* const bpm_distance_buffer) {
  gpu_buffer_bpm_distance_clear(bpm_distance_buffer->gpu_buffer_bpm_distance);
  vector_clear(bpm_distance_buffer->archive_searches);
}
void search_stage_bpm_distance_buffer_delete(
    search_stage_bpm_distance_buffer_t* const bpm_distance_buffer) {
  gpu_buffer_bpm_distance_delete(bpm_distance_buffer->gpu_buffer_bpm_distance);
  vector_delete(bpm_distance_buffer->archive_searches);
  mm_free(bpm_distance_buffer);
}
/*
 * Occupancy
 */
bool search_stage_bpm_distance_buffer_fits(
    search_stage_bpm_distance_buffer_t* const bpm_distance_buffer,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance =
      bpm_distance_buffer->gpu_buffer_bpm_distance;
  // Compute dimensions
  uint64_t total_queries = 0,total_query_entries = 0,total_candidates = 0;
  pattern_t* const pattern_end1 = &archive_search_end1->approximate_search.pattern;
  gpu_buffer_bpm_distance_compute_dimensions(
      gpu_buffer_bpm_distance,pattern_end1,
      archive_search_get_num_bpm_distance_candidates(archive_search_end1),
      &total_queries,&total_query_entries,&total_candidates);
  if (archive_search_end2!=NULL) {
    pattern_t* const pattern_end2 = &archive_search_end2->approximate_search.pattern;
    gpu_buffer_bpm_distance_compute_dimensions(
        gpu_buffer_bpm_distance,pattern_end2,
        archive_search_get_num_bpm_distance_candidates(archive_search_end2),
        &total_queries,&total_query_entries,&total_candidates);
  }
  // Return if current search fits in buffer
  return gpu_buffer_bpm_distance_fits_in_buffer(
      gpu_buffer_bpm_distance,total_queries,
      total_query_entries,total_candidates);
}
/*
 * Send/Receive
 */
void search_stage_bpm_distance_buffer_send(
    search_stage_bpm_distance_buffer_t* const bpm_distance_buffer) {
  PROF_ADD_COUNTER(GP_SSTAGE_BPM_DISTANCE_BUFFER_SEARCHES,
      vector_get_used(bpm_distance_buffer->archive_searches));
  gpu_buffer_bpm_distance_send(bpm_distance_buffer->gpu_buffer_bpm_distance);
}
void search_stage_bpm_distance_buffer_receive(
    search_stage_bpm_distance_buffer_t* const bpm_distance_buffer) {
  gpu_buffer_bpm_distance_receive(bpm_distance_buffer->gpu_buffer_bpm_distance);
}
/*
 * Accessors
 */
void search_stage_bpm_distance_buffer_add(
    search_stage_bpm_distance_buffer_t* const bpm_distance_buffer,
    archive_search_t* const archive_search) {
  vector_insert(bpm_distance_buffer->archive_searches,archive_search,archive_search_t*);
}
void search_stage_bpm_distance_buffer_retrieve(
    search_stage_bpm_distance_buffer_t* const bpm_distance_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search) {
  *archive_search = *vector_get_elm(bpm_distance_buffer->archive_searches,search_idx,archive_search_t*);
}

