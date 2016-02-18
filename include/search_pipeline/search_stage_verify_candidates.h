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
#include "archive/archive_search.h"
#include "archive/archive_search_cache.h"
#include "search_pipeline/search_stage_state.h"

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
  filtering_candidates_t filtering_candidates_forward_end1; // Filtering Candidates (end/1:F)
  filtering_candidates_t filtering_candidates_reverse_end1; // Filtering Candidates (end/1:R)
  filtering_candidates_t filtering_candidates_forward_end2; // Filtering Candidates (end/2:F)
  filtering_candidates_t filtering_candidates_reverse_end2; // Filtering Candidates (end/2:R)
  text_collection_t text_collection;                        // Stores text-traces
  mm_stack_t* mm_stack;                                     // MM-Stack
} search_stage_verify_candidates_t;

/*
 * Setup
 */
search_stage_verify_candidates_t* search_stage_verify_candidates_new(
    const gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const uint64_t num_buffers,
    const bool paired_end,
    const bool cpu_emulated,
    archive_text_t* const archive_text,
    mm_stack_t* const mm_stack);
void search_stage_verify_candidates_clear(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_cache_t* const archive_search_cache);
void search_stage_verify_candidates_delete(
    search_stage_verify_candidates_t* const search_stage_vc,
    archive_search_cache_t* const archive_search_cache);

/*
 * Accessors
 */
bool search_stage_verify_candidates_is_empty(search_stage_verify_candidates_t* const search_stage_vc);

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
