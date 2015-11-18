/*
 * PROJECT: GEMMapper
 * FILE: search_group_verify_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_STAGE_VERIFY_CANDIDATES_BUFFER_H_
#define SEARCH_STAGE_VERIFY_CANDIDATES_BUFFER_H_

#include "essentials.h"
#include "archive_search.h"
#include "archive_search_cache.h"
#include "gpu_buffer_align_bpm.h"

/*
 * Search-Stage Verify Candidates Buffer
 */
typedef struct {
  gpu_buffer_align_bpm_t* gpu_buffer_align_bpm; // GPU-BPM-Buffers
  vector_t* archive_searches;                   // Vector of archive-searches (archive_search_t*)
} search_stage_verify_candidates_buffer_t;

/*
 * Setup
 */
GEM_INLINE search_stage_verify_candidates_buffer_t* search_stage_verify_candidates_buffer_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffer_no,const bool cpu_emulated);
GEM_INLINE void search_stage_verify_candidates_buffer_clear(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    archive_search_cache_t* const archive_search_cache);
GEM_INLINE void search_stage_verify_candidates_buffer_delete(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    archive_search_cache_t* const archive_search_cache);

/*
 * Occupancy
 */
GEM_INLINE bool search_stage_verify_candidates_buffer_fits(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2);

/*
 * Send/Receive
 */
GEM_INLINE void search_stage_verify_candidates_buffer_send(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer);
GEM_INLINE void search_stage_verify_candidates_buffer_receive(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer);

/*
 * Accessors
 */
GEM_INLINE void search_stage_verify_candidates_buffer_add(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    archive_search_t* const archive_search);
GEM_INLINE void search_stage_verify_candidates_buffer_retrieve(
    search_stage_verify_candidates_buffer_t* const verify_candidates_buffer,
    const uint64_t search_idx,archive_search_t** const archive_search,
    matches_t* const matches);

#endif /* SEARCH_STAGE_VERIFY_CANDIDATES_BUFFER_H_ */
