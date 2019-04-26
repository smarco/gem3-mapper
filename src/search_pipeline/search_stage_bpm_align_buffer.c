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

#include "search_pipeline/search_stage_bpm_align_buffer.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Setup
 */
search_stage_bpm_align_buffer_t* search_stage_bpm_align_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool align_enabled) {
  // Alloc
  search_stage_bpm_align_buffer_t* const bpm_align_buffer = mm_alloc(search_stage_bpm_align_buffer_t);
  // Init
  bpm_align_buffer->gpu_buffer_bpm_align =
      gpu_buffer_bpm_align_new(gpu_buffer_collection,buffer_no,align_enabled);
  const uint64_t max_queries =
      gpu_buffer_bpm_align_get_max_queries(bpm_align_buffer->gpu_buffer_bpm_align);
  bpm_align_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return bpm_align_buffer;
}
void search_stage_bpm_align_buffer_clear(
    search_stage_bpm_align_buffer_t* const bpm_align_buffer) {
  gpu_buffer_bpm_align_clear(bpm_align_buffer->gpu_buffer_bpm_align);
  vector_clear(bpm_align_buffer->archive_searches);
}
void search_stage_bpm_align_buffer_delete(
    search_stage_bpm_align_buffer_t* const bpm_align_buffer) {
  gpu_buffer_bpm_align_delete(bpm_align_buffer->gpu_buffer_bpm_align);
  vector_delete(bpm_align_buffer->archive_searches);
  mm_free(bpm_align_buffer);
}
/*
 * Occupancy
 */
//TODO: Refactorize and implement proper interfaces
uint64_t search_stage_bpm_align_buffer_fits_max_candidates(
         archive_search_t* const archive_search){
 filtering_candidates_t* const filtering_candidates = archive_search->approximate_search.filtering_candidates;
 search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
 select_parameters_t* const select_parameters = &search_parameters->select_parameters;
 const uint64_t id_stats_histo = MIN(archive_search_get_num_bpm_align_canonical_candidates(archive_search), GEM_HIST_CAND_ALIGNED-1);
 const uint64_t num_samples_realigned = COUNTER_GET_NUM_SAMPLES(&filtering_candidates->candidates_aligned_histo[id_stats_histo]);
 const uint64_t average_samples_realigned = (uint64_t) ceil(COUNTER_GET_MEAN(&filtering_candidates->candidates_aligned_histo[id_stats_histo]));
 const uint64_t registered_samples_realigned = (num_samples_realigned==0) ? select_parameters->max_reported_matches : average_samples_realigned;
 const uint64_t max_buffered_candidates_aligned = MIN(registered_samples_realigned,archive_search_get_num_bpm_align_canonical_candidates(archive_search));
 return (max_buffered_candidates_aligned);
}

bool search_stage_bpm_align_buffer_fits(
    search_stage_bpm_align_buffer_t* const bpm_align_buffer,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align = bpm_align_buffer->gpu_buffer_bpm_align;
  // Compute dimensions
  uint64_t total_queries = 0, total_query_entries = 0, total_query_length = 0;
  uint64_t total_candidates = 0;
  // Updating the canonical candidates
  if (archive_search_end2!=NULL) {
    archive_search_pe_stepwise_bpm_align_update(archive_search_end1,archive_search_end2);
  }else{
    archive_search_se_stepwise_bpm_align_update(archive_search_end1);
  }
  // Calculate the minimum data structures to process the query on the device
  pattern_t* const pattern_end1 = &archive_search_end1->approximate_search.pattern;
  uint64_t num_canonical_regions = MIN(search_stage_bpm_align_buffer_fits_max_candidates(archive_search_end1), GEM_HIST_CAND_ALIGNED-1);
	//printf("idThread=%u, HISTO CANDIDATES=%d \n", (uint32_t)pthread_self(), num_canonical_regions);
	//fflush(stdout);
  gpu_buffer_bpm_align_compute_dimensions(
      gpu_buffer_bpm_align,pattern_end1,num_canonical_regions,
      &total_queries,&total_query_entries,&total_query_length,
      &total_candidates);
  if (archive_search_end2!=NULL) {
    pattern_t* const pattern_end2 = &archive_search_end2->approximate_search.pattern;
    uint64_t num_canonical_regions_end2 = MIN(search_stage_bpm_align_buffer_fits_max_candidates(archive_search_end2), GEM_HIST_CAND_ALIGNED-1);
    gpu_buffer_bpm_align_compute_dimensions(
        gpu_buffer_bpm_align,pattern_end2,num_canonical_regions_end2,
        &total_queries,&total_query_entries,&total_query_length,
        &total_candidates);
  }
  // Return if current search fits in buffer
  return gpu_buffer_bpm_align_fits_in_buffer(
      gpu_buffer_bpm_align,total_queries,
      total_query_entries,total_query_length,
      total_candidates);
}
/*
 * Send/Receive
 */
void search_stage_bpm_align_buffer_send(
    search_stage_bpm_align_buffer_t* const bpm_align_buffer) {
  PROF_ADD_COUNTER(GP_SSTAGE_BPM_ALIGN_BUFFER_SEARCHES,
      vector_get_used(bpm_align_buffer->archive_searches));
  gpu_buffer_bpm_align_send(bpm_align_buffer->gpu_buffer_bpm_align);
}
void search_stage_bpm_align_buffer_receive(
    search_stage_bpm_align_buffer_t* const bpm_align_buffer) {
  gpu_buffer_bpm_align_receive(bpm_align_buffer->gpu_buffer_bpm_align);
}
/*
 * Accessors
 */
void search_stage_bpm_align_buffer_add(
    search_stage_bpm_align_buffer_t* const bpm_align_buffer,
    archive_search_t* const archive_search) {
  vector_insert(bpm_align_buffer->archive_searches,archive_search,archive_search_t*);
}
void search_stage_bpm_align_buffer_retrieve(
    search_stage_bpm_align_buffer_t* const bpm_align_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search) {
  *archive_search = *vector_get_elm(bpm_align_buffer->archive_searches,search_idx,archive_search_t*);
}

