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

#include "search_pipeline/search_stage_region_profile_buffer.h"
#include "archive/search/archive_search_se_stepwise.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Setup
 */
search_stage_region_profile_buffer_t* search_stage_region_profile_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    const bool region_profile_enabled,
    const uint32_t occ_min_threshold,
    const uint32_t extra_search_steps,
    const uint32_t alphabet_size) {
  // Alloc
  search_stage_region_profile_buffer_t* const region_profile_buffer =
      mm_alloc(search_stage_region_profile_buffer_t);
  // Init gpu-buffer
#ifdef GPU_REGION_PROFILE_ADAPTIVE
  region_profile_buffer->gpu_buffer_fmi_asearch = gpu_buffer_fmi_asearch_new(
      gpu_buffer_collection,buffer_no,region_profile_enabled,
      occ_min_threshold,extra_search_steps,alphabet_size);
  const uint64_t max_queries =
      gpu_buffer_fmi_asearch_get_max_queries(region_profile_buffer->gpu_buffer_fmi_asearch);
#else
  region_profile_buffer->gpu_buffer_fmi_ssearch = gpu_buffer_fmi_ssearch_new(
      gpu_buffer_collection,buffer_no,region_profile_enabled);
  const uint64_t max_queries =
      gpu_buffer_fmi_ssearch_get_max_queries(region_profile_buffer->gpu_buffer_fmi_ssearch);
#endif
  // Allocate searches-buffer
  region_profile_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return region_profile_buffer;
}
void search_stage_region_profile_buffer_clear(
    search_stage_region_profile_buffer_t* const region_profile_buffer) {
#ifdef GPU_REGION_PROFILE_ADAPTIVE
  gpu_buffer_fmi_asearch_clear(region_profile_buffer->gpu_buffer_fmi_asearch);
#else
  gpu_buffer_fmi_ssearch_clear(region_profile_buffer->gpu_buffer_fmi_ssearch);
#endif
  // Clear searches vector
  vector_clear(region_profile_buffer->archive_searches);
}
void search_stage_region_profile_buffer_delete(
    search_stage_region_profile_buffer_t* const region_profile_buffer) {
#ifdef GPU_REGION_PROFILE_ADAPTIVE
  gpu_buffer_fmi_asearch_delete(region_profile_buffer->gpu_buffer_fmi_asearch);
#else
  gpu_buffer_fmi_ssearch_delete(region_profile_buffer->gpu_buffer_fmi_ssearch);
#endif
  // Delete searches vector
  vector_delete(region_profile_buffer->archive_searches);
  mm_free(region_profile_buffer);
}
/*
 * Occupancy
 */
bool search_stage_region_profile_buffer_fits(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
#ifdef GPU_REGION_PROFILE_ADAPTIVE
  // Compute dimensions
  uint64_t total_queries = 1;
  uint64_t total_bases = archive_search_end1->approximate_search.pattern.key_length;
  uint64_t total_regions = archive_search_end1->approximate_search.region_profile.max_expected_regions;
  if (archive_search_end2!=NULL) {
    ++total_queries;
    total_bases += archive_search_end2->approximate_search.pattern.key_length;
    total_regions += archive_search_end2->approximate_search.region_profile.max_expected_regions;
  }
  // Return if current search fits in buffer
  return gpu_buffer_fmi_asearch_fits_in_buffer(
      region_profile_buffer->gpu_buffer_fmi_asearch,
      total_queries,total_bases,total_regions);
#else
  // Compute dimensions
  uint64_t num_regions_profile;
  num_regions_profile = archive_search_get_num_regions_profile(archive_search_end1);
  if (archive_search_end2!=NULL) {
    num_regions_profile += archive_search_get_num_regions_profile(archive_search_end2);
  }
  // Return if current search fits in buffer
  return gpu_buffer_fmi_ssearch_fits_in_buffer(
      region_profile_buffer->gpu_buffer_fmi_ssearch,num_regions_profile);
#endif
}
/*
 * Send/Receive
 */
void search_stage_region_profile_buffer_send(
    search_stage_region_profile_buffer_t* const region_profile_buffer) {
  PROF_ADD_COUNTER(GP_SSTAGE_REGION_PROFILE_BUFFER_SEARCHES,
      vector_get_used(region_profile_buffer->archive_searches));
#ifdef GPU_REGION_PROFILE_ADAPTIVE
  gpu_buffer_fmi_asearch_send(region_profile_buffer->gpu_buffer_fmi_asearch);
#else
  gpu_buffer_fmi_ssearch_send(region_profile_buffer->gpu_buffer_fmi_ssearch);
#endif
}
void search_stage_region_profile_buffer_receive(
    search_stage_region_profile_buffer_t* const region_profile_buffer) {
#ifdef GPU_REGION_PROFILE_ADAPTIVE
  gpu_buffer_fmi_asearch_receive(region_profile_buffer->gpu_buffer_fmi_asearch);
#else
  gpu_buffer_fmi_ssearch_receive(region_profile_buffer->gpu_buffer_fmi_ssearch);
#endif
}
/*
 * Accessors
 */
void search_stage_region_profile_buffer_add(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_t* const archive_search) {
  vector_insert(region_profile_buffer->archive_searches,archive_search,archive_search_t*);
}
void search_stage_region_profile_buffer_retrieve(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search) {
  *archive_search = *vector_get_elm(region_profile_buffer->archive_searches,search_idx,archive_search_t*);
}


