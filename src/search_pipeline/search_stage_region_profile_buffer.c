/*
 * PROJECT: GEMMapper
 * FILE: search_stage_region_profile_buffer.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_pipeline/search_stage_region_profile_buffer.h"
#include "archive/archive_search_se_stepwise.h"

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
    fm_index_t* const fm_index,
    const bool cpu_emulated) {
  // Alloc
  search_stage_region_profile_buffer_t* const region_profile_buffer = mm_alloc(search_stage_region_profile_buffer_t);
  // Init
  region_profile_buffer->gpu_buffer_fmi_search =
      gpu_buffer_fmi_search_new(gpu_buffer_collection,buffer_no,fm_index);
  if (cpu_emulated) gpu_buffer_fmi_search_set_device_cpu(region_profile_buffer->gpu_buffer_fmi_search);
  const uint64_t max_queries = gpu_buffer_fmi_search_get_max_queries(region_profile_buffer->gpu_buffer_fmi_search);
  region_profile_buffer->archive_searches = vector_new(max_queries,archive_search_t*);
  // Return
  return region_profile_buffer;
}
void search_stage_region_profile_buffer_clear(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_cache_t* const archive_search_cache) {
  gpu_buffer_fmi_search_clear(region_profile_buffer->gpu_buffer_fmi_search);
  // Return searches to the cache
  if (archive_search_cache!=NULL) {
    VECTOR_ITERATE(region_profile_buffer->archive_searches,archive_search,n,archive_search_t*) {
      archive_search_cache_free(archive_search_cache,*archive_search);
    }
  }
  // Clear searches vector
  vector_clear(region_profile_buffer->archive_searches);
}
void search_stage_region_profile_buffer_delete(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_cache_t* const archive_search_cache) {
  gpu_buffer_fmi_search_delete(region_profile_buffer->gpu_buffer_fmi_search);
  // Return searches to the cache
  VECTOR_ITERATE(region_profile_buffer->archive_searches,archive_search,n,archive_search_t*) {
    archive_search_cache_free(archive_search_cache,*archive_search);
  }
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
  // Compute dimensions (Total number of regions to profile)
  uint64_t num_regions_profile;
  num_regions_profile = archive_search_get_num_regions_profile(archive_search_end1);
  if (archive_search_end2!=NULL) {
    num_regions_profile += archive_search_get_num_regions_profile(archive_search_end2);
  }
  // Return if current search fits in buffer
  return gpu_buffer_fmi_search_fits_in_buffer(region_profile_buffer->gpu_buffer_fmi_search,num_regions_profile);
}
/*
 * Send/Receive
 */
void search_stage_region_profile_buffer_send(
    search_stage_region_profile_buffer_t* const region_profile_buffer) {
  PROF_ADD_COUNTER(GP_SEARCH_STAGE_REGION_PROFILE_SEARCHES_IN_BUFFER,
      vector_get_used(region_profile_buffer->archive_searches));
  gpu_buffer_fmi_search_send(region_profile_buffer->gpu_buffer_fmi_search);
}
void search_stage_region_profile_buffer_receive(
    search_stage_region_profile_buffer_t* const region_profile_buffer) {
  gpu_buffer_fmi_search_receive(region_profile_buffer->gpu_buffer_fmi_search);
}
/*
 * Accessors
 */
void search_stage_region_profile_buffer_add(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_t* const archive_search) {
  // Add archive-search
  vector_insert(region_profile_buffer->archive_searches,archive_search,archive_search_t*);
}
void search_stage_region_profile_buffer_retrieve(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search) {
  // Retrieve archive-search
  *archive_search = *vector_get_elm(region_profile_buffer->archive_searches,search_idx,archive_search_t*);
}


