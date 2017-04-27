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

#include "search_pipeline/search_stage_kmer_filter_buffer.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Setup
 */
search_stage_kmer_filter_buffer_t* search_stage_kmer_filter_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool kmer_filter_enabled) {
  // Alloc
  search_stage_kmer_filter_buffer_t* const kmer_filter_buffer =
      mm_alloc(search_stage_kmer_filter_buffer_t);
  // Init
  kmer_filter_buffer->gpu_buffer_kmer_filter =
      gpu_buffer_kmer_filter_new(gpu_buffer_collection,buffer_no,kmer_filter_enabled);
  const uint64_t max_queries = gpu_buffer_kmer_filter_get_max_queries(kmer_filter_buffer->gpu_buffer_kmer_filter);
  kmer_filter_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return kmer_filter_buffer;
}
void search_stage_kmer_filter_buffer_clear(
    search_stage_kmer_filter_buffer_t* const kmer_filter_buffer) {
  gpu_buffer_kmer_filter_clear(kmer_filter_buffer->gpu_buffer_kmer_filter);
  // Clear searches vector
  vector_clear(kmer_filter_buffer->archive_searches);
}
void search_stage_kmer_filter_buffer_delete(
    search_stage_kmer_filter_buffer_t* const kmer_filter_buffer) {
  gpu_buffer_kmer_filter_delete(kmer_filter_buffer->gpu_buffer_kmer_filter);
  // Delete searches vector
  vector_delete(kmer_filter_buffer->archive_searches);
  mm_free(kmer_filter_buffer);
}
/*
 * Occupancy
 */
bool search_stage_kmer_filter_buffer_fits(
    search_stage_kmer_filter_buffer_t* const kmer_filter_buffer,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter = kmer_filter_buffer->gpu_buffer_kmer_filter;
  // Compute dimensions
  uint64_t total_queries = 0, total_queries_length = 0, total_candidates = 0;
  pattern_t* const pattern_end1 = &archive_search_end1->approximate_search.pattern;
  kmer_counting_nway_t* const kmer_counting_nway_end1 = &pattern_end1->pattern_tiled.kmer_filter_nway;
  gpu_buffer_kmer_filter_compute_dimensions(
      gpu_buffer_kmer_filter,kmer_counting_nway_end1,
      archive_search_get_num_kmer_filter_candidates(archive_search_end1),
      &total_queries,&total_queries_length,&total_candidates);
  if (archive_search_end2!=NULL) {
    pattern_t* const pattern_end2 = &archive_search_end2->approximate_search.pattern;
    kmer_counting_nway_t* const kmer_counting_nway_end2 = &pattern_end2->pattern_tiled.kmer_filter_nway;
    gpu_buffer_kmer_filter_compute_dimensions(
        gpu_buffer_kmer_filter,kmer_counting_nway_end2,
        archive_search_get_num_kmer_filter_candidates(archive_search_end2),
        &total_queries,&total_queries_length,&total_candidates);
  }
  // Return if current search fits in buffer
  return gpu_buffer_kmer_filter_fits_in_buffer(
      gpu_buffer_kmer_filter,total_queries,total_queries_length,total_candidates);
}
/*
 * Send/Receive
 */
void search_stage_kmer_filter_buffer_send(
    search_stage_kmer_filter_buffer_t* const kmer_filter_buffer) {
  PROF_ADD_COUNTER(GP_SSTAGE_KMER_FILTER_BUFFER_SEARCHES,
      vector_get_used(kmer_filter_buffer->archive_searches));
  gpu_buffer_kmer_filter_send(kmer_filter_buffer->gpu_buffer_kmer_filter);
}
void search_stage_kmer_filter_buffer_receive(
    search_stage_kmer_filter_buffer_t* const kmer_filter_buffer) {
  gpu_buffer_kmer_filter_receive(kmer_filter_buffer->gpu_buffer_kmer_filter);
}
/*
 * Accessors
 */
void search_stage_kmer_filter_buffer_add(
    search_stage_kmer_filter_buffer_t* const kmer_filter_buffer,
    archive_search_t* const archive_search) {
  vector_insert(kmer_filter_buffer->archive_searches,archive_search,archive_search_t*);
}
void search_stage_kmer_filter_buffer_retrieve(
    search_stage_kmer_filter_buffer_t* const kmer_filter_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search) {
  *archive_search = *vector_get_elm(kmer_filter_buffer->archive_searches,search_idx,archive_search_t*);
}

