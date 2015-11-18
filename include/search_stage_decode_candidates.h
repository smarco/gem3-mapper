/*
 * PROJECT: GEMMapper
 * FILE: search_state_decode_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_STATE_DECODE_CANDIDATES_H_
#define SEARCH_STATE_DECODE_CANDIDATES_H_


#include "essentials.h"
#include "archive_search.h"
#include "archive_search_cache.h"
#include "search_stage_state.h"

/*
 * Search-group Decode Candidates
 */
typedef struct {
  // Configuration
  search_stage_mode_t search_stage_mode;  // Stage Mode (Sending/Receiving)
  // Decode Candidates Buffers
  vector_t* buffers;                      // Verify Candidates Buffers (search_stage_decode_candidates_buffer_t*)
  search_stage_iterator_t iterator;       // Buffers State
} search_stage_decode_candidates_t;

/*
 * Setup
 */
search_stage_decode_candidates_t* search_stage_decode_candidates_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,const uint64_t num_buffers,
    const bool cpu_emulated);
void search_stage_decode_candidates_clear(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_cache_t* const archive_search_cache);
void search_stage_decode_candidates_delete(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_cache_t* const archive_search_cache);

/*
 * Accessors
 */
bool search_stage_decode_candidates_is_empty(search_stage_decode_candidates_t* const search_stage_dc);

/*
 * Send Searches (buffered)
 */
bool search_stage_decode_candidates_send_se_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t* const archive_search);
bool search_stage_decode_candidates_send_pe_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Retrieve Searches (buffered)
 */
bool search_stage_decode_candidates_retrieve_se_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t** const archive_search);
bool search_stage_decode_candidates_retrieve_pe_search(
    search_stage_decode_candidates_t* const search_stage_dc,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);


#endif /* SEARCH_STATE_DECODE_CANDIDATES_H_ */
