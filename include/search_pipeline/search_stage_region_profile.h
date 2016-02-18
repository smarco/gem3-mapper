/*
 * PROJECT: GEMMapper
 * FILE: search_state_region_profile.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_STATE_REGION_PROFILE_H_
#define SEARCH_STATE_REGION_PROFILE_H_

#include "utils/essentials.h"
#include "archive/archive_search.h"
#include "archive/archive_search_cache.h"
#include "search_pipeline/search_stage_state.h"

/*
 * Search-Stage Region-Profile
 */
typedef struct {
  // Configuration
  search_stage_mode_t search_stage_mode;  // Stage Mode (Sending/Receiving)
  // Region-Profile Buffers
  vector_t* buffers;                      // Region-Profile Buffers (search_stage_region_profile_buffer_t*)
  search_stage_iterator_t iterator;       // Buffers Iterator
} search_stage_region_profile_t;

/*
 * Setup
 */
search_stage_region_profile_t* search_stage_region_profile_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const uint64_t num_buffers,
    fm_index_t* const fm_index,
    const bool cpu_emulated);
void search_stage_region_profile_clear(
    search_stage_region_profile_t* const search_stage_rp,
    archive_search_cache_t* const archive_search_cache);
void search_stage_region_profile_delete(
    search_stage_region_profile_t* const search_stage_rp,
    archive_search_cache_t* const archive_search_cache);

/*
 * Accessors
 */
bool search_stage_region_profile_is_empty(search_stage_region_profile_t* const search_stage_rp);

/*
 * Send Searches (buffered)
 */
bool search_stage_region_profile_send_se_search(
    search_stage_region_profile_t* const search_stage_rp,
    archive_search_t* const archive_search);
bool search_stage_region_profile_send_pe_search(
    search_stage_region_profile_t* const search_stage_rp,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Retrieve Searches (buffered)
 */
bool search_stage_region_profile_retrieve_finished(
    search_stage_region_profile_t* const search_stage_rp);
bool search_stage_region_profile_retrieve_se_search(
    search_stage_region_profile_t* const search_stage_rp,
    archive_search_t** const archive_search);
bool search_stage_region_profile_retrieve_pe_search(
    search_stage_region_profile_t* const search_stage_rp,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);

#endif /* SEARCH_STATE_REGION_PROFILE_H_ */
