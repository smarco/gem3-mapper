/*
 * PROJECT: GEMMapper
 * FILE: search_group_verify_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_GROUP_VERIFY_CANDIDATES_H_
#define SEARCH_GROUP_VERIFY_CANDIDATES_H_

#include "essentials.h"
#include "archive_search.h"
#include "archive_search_cache.h"
#include "archive_search_se_stepwise.h"

typedef struct {
  /* BPM-GPU candidates buffer */
  bpm_gpu_buffer_t* bpm_gpu_buffer;      // BPM-Buffer
  /* Archive searches */
  vector_t* archive_searches;            // Vector of archive-searches (archive_search_t*)
} search_group_verify_candidates_buffer_t;

typedef struct {
  /* Configuration */
  bool cpu_emulated;                                       // Emulate verification using CPU
  uint64_t hint_patterns_per_search;                       // Number of patterns per search-group
  /* Verification buffers */
  search_group_verify_candidates_buffer_t* buffers;        // Search-group buffers
  uint64_t num_buffers;                                    // Total number of search-groups
  uint64_t num_archive_searches;                           // Total number of archive-searches
  /* Verification buffers iterator */
  uint64_t current_buffer_idx;                             // Current search-group index
  search_group_verify_candidates_buffer_t* current_buffer; // Current search-group
  /* Archive-Search iterator */
  uint64_t current_archive_search_idx;                     // Current archive-search index
  archive_search_t** current_archive_search;               // Current archive-search
} search_group_verify_candidates_t;

/*
 * Setup
 */
search_group_verify_candidates_t* search_group_verify_candidates_new(
    bpm_gpu_buffer_t* const bpm_gpu_buffers,const uint64_t num_buffers,
    const bool cpu_emulated);
void search_group_verify_candidates_init(search_group_verify_candidates_t* const search_group_vc);
void search_group_verify_candidates_clear(
    search_group_verify_candidates_t* const search_group_vc,
    archive_search_cache_t* const archive_search_cache);
void search_group_verify_candidates_delete(search_group_verify_candidates_t* const search_group_vc);

/*
 * Accessors
 */
bool search_group_verify_candidates_buffer_is_empty(search_group_verify_candidates_t* const search_group_vc);
void search_group_verify_candidates_retrieve_begin(search_group_verify_candidates_t* const search_group_vc);

/*
 * SE search-group verification
 */
bool search_group_verify_candidates_se_add_search(
    search_group_verify_candidates_t* const search_group_vc,archive_search_t* const archive_search);
bool search_group_verify_candidates_se_get_search(
    search_group_verify_candidates_t* const search_group_vc,archive_search_t** const archive_search,
    text_collection_t* const text_collection,matches_t* const matches);

/*
 * PE search-group verification
 */
bool search_group_verify_candidates_pe_add_search(
    search_group_verify_candidates_t* const search_group_vc,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2);
bool search_group_verify_candidates_pe_get_search(
    search_group_verify_candidates_t* const search_group_vc,archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2,text_collection_t* const text_collection,
    paired_matches_t* const paired_matches);

/*
 * Error Messages
 */
#define GEM_ERROR_SEARCH_GROUP_VERIFY_CANDIDATES_QUERY_TOO_BIG  "Search-Group verify candidates. Couldn't copy query to BPM-buffer (Query too big)"
#define GEM_ERROR_SEARCH_GROUP_VERIFY_CANDIDATES_UNPAIRED_QUERY "Search-Group verify candidates. Couldn't retrieve query-pair"

#endif /* SEARCH_GROUP_VERIFY_CANDIDATES_H_ */
