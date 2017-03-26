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

#include "search_pipeline/search_stage_decode_buffer.h"

/*
 * Setup
 */
search_stage_decode_buffer_t* search_stage_decode_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const uint32_t sampling_rate,
    const bool decode_sa_enabled,
    const bool decode_text_enabled) {
  // Alloc
  search_stage_decode_buffer_t* const decode_buffer = mm_alloc(search_stage_decode_buffer_t);
  // Init
  decode_buffer->gpu_buffer_fmi_decode = gpu_buffer_fmi_decode_new(
      gpu_buffer_collection,buffer_no,sampling_rate,decode_sa_enabled,decode_text_enabled);
  const uint64_t max_queries = gpu_buffer_fmi_decode_get_max_queries(decode_buffer->gpu_buffer_fmi_decode);
  decode_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return decode_buffer;
}
void search_stage_decode_buffer_clear(
    search_stage_decode_buffer_t* const decode_buffer) {
  gpu_buffer_fmi_decode_clear(decode_buffer->gpu_buffer_fmi_decode);
  // Clear searches vector
  vector_clear(decode_buffer->archive_searches);
}
void search_stage_decode_buffer_delete(
    search_stage_decode_buffer_t* const decode_buffer) {
  gpu_buffer_fmi_decode_delete(decode_buffer->gpu_buffer_fmi_decode);
  // Delete searches vector
  vector_delete(decode_buffer->archive_searches);
  mm_free(decode_buffer);
}
/*
 * Occupancy
 */
bool search_stage_decode_buffer_fits(
    search_stage_decode_buffer_t* const decode_buffer,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  // Compute dimensions (Total number of candidates to decode)
  uint64_t num_decode = archive_search_get_num_decode_candidates(archive_search_end1);
  if (archive_search_end2!=NULL) {
    num_decode += archive_search_get_num_decode_candidates(archive_search_end2);
  }
  // Return if current search fits in buffer
  return gpu_buffer_fmi_decode_fits_in_buffer(decode_buffer->gpu_buffer_fmi_decode,num_decode);
}
/*
 * Send/Receive
 */
void search_stage_decode_buffer_send(
    search_stage_decode_buffer_t* const decode_buffer) {
  PROF_ADD_COUNTER(GP_SSTAGE_DECODE_BUFFER_SEARCHES,
      vector_get_used(decode_buffer->archive_searches));
  gpu_buffer_fmi_decode_send(decode_buffer->gpu_buffer_fmi_decode);
}
void search_stage_decode_buffer_receive(
    search_stage_decode_buffer_t* const decode_buffer) {
  gpu_buffer_fmi_decode_receive(decode_buffer->gpu_buffer_fmi_decode);
}
/*
 * Accessors
 */
void search_stage_decode_buffer_add(
    search_stage_decode_buffer_t* const decode_buffer,
    archive_search_t* const archive_search) {
  vector_insert(decode_buffer->archive_searches,archive_search,archive_search_t*);
}
void search_stage_decode_buffer_retrieve(
    search_stage_decode_buffer_t* const decode_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search) {
  *archive_search = *vector_get_elm(decode_buffer->archive_searches,search_idx,archive_search_t*);
}


