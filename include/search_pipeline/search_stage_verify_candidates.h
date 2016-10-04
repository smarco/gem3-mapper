/*
 * PROJECT: GEMMapper
 * FILE: search_group_verify_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_STAGE_VERIFY_CANDIDATES_H_
#define SEARCH_STAGE_VERIFY_CANDIDATES_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_cache.h"
#include "search_pipeline/search_stage_state.h"
#include "search_pipeline/search_pipeline_handlers.h"

/*
 * Search-group Verify Candidates
 */
typedef struct {
  // Configuration
  bool paired_end;                                          // Paired-end search
  search_stage_mode_t search_stage_mode;                    // Stage Mode (Sending/Receiving)
  // Verify Candidates Buffers
  vector_t* buffers;                                        // Verify Candidates Buffers (search_group_buffer_vc_t*)
  search_stage_iterator_t iterator;                         // Buffers State
  /* Support Data Structures */
  matches_t* matches;                                       // Matches
  paired_matches_t* paired_matches;                         // Paired-Matches
  search_pipeline_handlers_t* search_pipeline_handlers;
} search_stage_verify_candidates_t;

/*
 * Setup
 */
search_stage_verify_candidates_t* search_stage_verify_candidates_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const uint64_t num_buffers,
    const bool paired_end,
    const bool verify_candidates_enabled,
    search_pipeline_handlers_t* const search_pipeline_handlers);
void search_stage_verify_candidates_clear(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_cache_t* const archive_search_cache);
void search_stage_verify_candidates_delete(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_cache_t* const archive_search_cache);

/*
 * Send Searches (buffered)
 */
bool search_stage_verify_candidates_send_se_search(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_t* const archive_search);
bool search_stage_verify_candidates_send_pe_search(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Retrieve Searches (buffered)
 */
bool search_stage_verify_candidates_retrieve_finished(
    search_stage_verify_candidates_t* const search_stage_vc);
bool search_stage_verify_candidates_retrieve_se_search(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_t** const archive_search);
bool search_stage_verify_candidates_retrieve_pe_search(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);

#endif /* SEARCH_STAGE_VERIFY_CANDIDATES_H_ */
