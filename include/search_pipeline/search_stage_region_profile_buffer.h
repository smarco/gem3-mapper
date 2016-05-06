/*
 * PROJECT: GEMMapper
 * FILE: search_state_region_profile_buffer.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_STATE_REGION_PROFILE_BUFFER_H_
#define SEARCH_STATE_REGION_PROFILE_BUFFER_H_

#include "utils/essentials.h"
#include "archive/archive_search.h"
#include "archive/archive_search_cache.h"
#include "gpu/gpu_buffer_fmi_bsearch.h"

/*
 * Search-Stage Region Profile Buffer
 */
typedef struct {
  gpu_buffer_fmi_search_t* gpu_buffer_fmi_search; // GPU-BPM-Buffer
  vector_t* archive_searches;                     // Vector of archive-searches (archive_search_t*)
} search_stage_region_profile_buffer_t;

/*
 * Setup
 */
search_stage_region_profile_buffer_t* search_stage_region_profile_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,
    fm_index_t* const fm_index,
    const bool cpu_emulated);
void search_stage_region_profile_buffer_clear(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_cache_t* const archive_search_cache);
void search_stage_region_profile_buffer_delete(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_cache_t* const archive_search_cache);

/*
 * Occupancy
 */
bool search_stage_region_profile_buffer_fits(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Send/Receive
 */
void search_stage_region_profile_buffer_send(
    search_stage_region_profile_buffer_t* const region_profile_buffer);
void search_stage_region_profile_buffer_receive(
    search_stage_region_profile_buffer_t* const region_profile_buffer);

/*
 * Accessors
 */
void search_stage_region_profile_buffer_add(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    archive_search_t* const archive_search);
void search_stage_region_profile_buffer_retrieve(
    search_stage_region_profile_buffer_t* const region_profile_buffer,
    const uint64_t search_idx,
    archive_search_t** const archive_search);

#endif /* SEARCH_STATE_REGION_PROFILE_BUFFER_H_ */
