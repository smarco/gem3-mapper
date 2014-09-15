/*
 * PROJECT: GEMMapper
 * FILE: archive_search_group.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_GROUP_H_
#define ARCHIVE_SEARCH_GROUP_H_

#include "essentials.h"
#include "archive_search.h"
#include "bpm_align_gpu.h"

/*
 * Checker
 */
// TODO

/*
 * Mapper Search Dispatcher (Gathers all BPM-GPU buffers and related mapper searches)
 */
typedef struct {
  /* Search */
  archive_search_t* archive_search;     // Archive search
  /* Candidates verified */
  uint64_t results_buffer_offset;       // Offset in the results vector
} archive_search_member_t;
typedef enum {
  archive_search_group_free,
  archive_search_group_generating_candidates,
  archive_search_group_verifying_candidates,
  archive_search_group_selecting_candidates
} archive_search_group_state_t;
typedef struct {
  /* State */
  uint32_t mayor_group_id; // TODO Group mayor ID (for input block)
  uint32_t minor_group_id; // TODO Group minor ID (for segmented groups)
  archive_search_group_state_t state;                // State of the search group
  /* BPM-GPU candidates buffer */
  bpm_gpu_buffer_t* bpm_gpu_buffer;                  // BPM-Buffer
  /* Archive searches */
  vector_t* archive_searches;                        // Vector of search members (archive_search_member_t)
  vector_t* archive_searches_attr;                   // TODO Vector of generic attr of the members (E.g related buffered_ouput)
  /* MM */
  mm_search_t* mm_search_end1;                       // Memory Managed Search
  mm_search_t* mm_search_end2;                       // Memory Managed Search
} archive_search_group_t;
typedef struct {
  /* Dispatcher State */
  uint64_t num_groups;                  // Total number of search-groups allocated
  uint64_t num_groups_free;             // Free groups (ready to be filled & used)
  uint64_t num_groups_generating;       // Groups dispatched for candidate generation
  uint64_t num_groups_verifying;        // Groups being verified (BPM-CUDA)
  uint64_t num_groups_selecting;        // Groups dispatched for candidate selection
  /* Dispatcher Record */
  uint64_t num_threads_generating;
  /* Dispatcher Search Groups */
  archive_search_group_t* search_group; // Search Groups
  /* BPM-GPU Buffer*/
  bpm_gpu_buffer_collection_t* bpm_gpu_buffer_collection; // BPM Buffers
  /* Mutex/CV */
  pthread_mutex_t dispatcher_mutex;
  pthread_cond_t groups_free_cond;
  pthread_cond_t groups_verifying_cond;
} archive_search_group_dispatcher_t;

/*
 * Archive-search group
 */
GEM_INLINE void archive_search_group_add(
    archive_search_group_t* const archive_search_group,
    archive_search_t* const archive_search,const uint64_t results_buffer_offset);

/*
 * Dispatcher
 */
GEM_INLINE archive_search_group_dispatcher_t* archive_search_group_dispatcher_new();
GEM_INLINE void archive_search_group_dispatcher_delete(archive_search_group_dispatcher_t* const dispatcher);

GEM_INLINE void archive_search_group_dispatcher_register_generating(
    archive_search_group_dispatcher_t* const dispatcher,const uint64_t num_threads);
GEM_INLINE void archive_search_group_dispatcher_deregister_generating(
    archive_search_group_dispatcher_t* const dispatcher,const uint64_t num_threads);

GEM_INLINE archive_search_group_t* archive_search_group_dispatcher_request_generating(
    archive_search_group_dispatcher_t* const dispatcher);
GEM_INLINE void archive_search_group_dispatcher_return_generating(
    archive_search_group_dispatcher_t* const dispatcher,archive_search_group_t* const search_group);
GEM_INLINE archive_search_group_t* archive_search_group_dispatcher_request_selecting(
    archive_search_group_dispatcher_t* const dispatcher);
GEM_INLINE void archive_search_group_dispatcher_return_selecting(
    archive_search_group_dispatcher_t* const dispatcher,archive_search_group_t* const search_group);

/*
 * Step-wise SE-Search
 */
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search);
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer);
GEM_INLINE void archive_search_select_candidates(
    archive_search_t* const archive_search,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t results_buffer_offset);

/*
 * Error Messages
 */
#define GEM_ERROR_ARCHIVE_SEARCH_GROUP_MAPPING_MODE_NOT_SUPPORTED "Archive search-group. Mapping mode not supported (adaptive | fixed)"

#endif /* ARCHIVE_SEARCH_GROUP_H_ */
