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
#include "buffered_output_file.h"
#include "mapper.h"

typedef struct _search_group_t search_group_t;
typedef struct _search_group_dispatcher_t search_group_dispatcher_t;

/*
 * Archive-search group
 */
GEM_INLINE void search_group_set_incomplete(search_group_t* const search_group);
GEM_INLINE bool search_group_is_incomplete(search_group_t* const search_group);

GEM_INLINE uint64_t search_group_get_group_id(search_group_t* const search_group);
GEM_INLINE bpm_gpu_buffer_t* search_group_get_bpm_buffer(search_group_t* const search_group);
GEM_INLINE uint64_t search_group_get_num_searches(search_group_t* const search_group);
GEM_INLINE void search_group_get_search(
    search_group_t* const search_group,const uint64_t position,
    archive_search_t** const archive_search,uint64_t* const results_buffer_offset);

GEM_INLINE void search_group_add_search(
    search_group_t* const search_group,
    archive_search_t* const archive_search,const uint64_t results_buffer_offset);

// Archive Search Group Allocator (Cache)
GEM_INLINE archive_search_t* search_group_alloc(search_group_t* const search_group);
GEM_INLINE void search_group_release(search_group_t* const search_group,archive_search_t* const archive_search);
GEM_INLINE void search_group_configure(search_group_t* const search_group,archive_search_t* const archive_search);

/*
 * Dispatcher
 */
GEM_INLINE search_group_dispatcher_t* search_group_dispatcher_new(
    mapper_parameters_t* const mapper_parameters,archive_t* const archive,
    const uint64_t num_search_groups,const uint64_t bpm_buffer_size,
    const uint64_t average_query_size,const uint64_t candidates_per_query);
GEM_INLINE void search_group_dispatcher_delete(search_group_dispatcher_t* const dispatcher);

// Register/Deregister Threads
GEM_INLINE void search_group_dispatcher_register_generating(
    search_group_dispatcher_t* const dispatcher,const uint64_t num_threads);
GEM_INLINE void search_group_dispatcher_deregister_generating(
    search_group_dispatcher_t* const dispatcher,const uint64_t num_threads);

// Generating group
GEM_INLINE search_group_t* search_group_dispatcher_request_generating(
    search_group_dispatcher_t* const dispatcher,const uint32_t request_id);
GEM_INLINE search_group_t* search_group_dispatcher_request_generating_extension(
    search_group_dispatcher_t* const dispatcher,search_group_t* const search_group);
GEM_INLINE void search_group_dispatcher_return_generating(
    search_group_dispatcher_t* const dispatcher,search_group_t* const search_group);

// Selecting group
GEM_INLINE search_group_t* search_group_dispatcher_request_selecting(
    search_group_dispatcher_t* const dispatcher);
GEM_INLINE search_group_t* search_group_dispatcher_request_selecting_next(
    search_group_dispatcher_t* const dispatcher,search_group_t* search_group);
GEM_INLINE void search_group_dispatcher_return_selecting(
    search_group_dispatcher_t* const dispatcher,search_group_t* const search_group);

/*
 * Step-wise SE-Search
 */
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search);
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer);
GEM_INLINE void archive_search_retrieve_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t results_buffer_offset,matches_t* const matches);

/*
 * Error Messages
 */
#define GEM_ERROR_ARCHIVE_SEARCH_GROUP_MAPPING_MODE_NOT_SUPPORTED "Archive search-group. Mapping mode not supported (adaptive | fixed)"

#endif /* ARCHIVE_SEARCH_GROUP_H_ */
