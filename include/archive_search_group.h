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

/*
 * Archive-Search groups
 */
typedef struct _archive_search_group_t archive_search_group_t;

/*
 * Archive-Search group
 */
GEM_INLINE archive_search_group_t* archive_search_group_new(
    mapper_parameters_t* const mapper_parameters,bpm_gpu_buffer_t* const bpm_gpu_buffers,const uint64_t total_search_groups);
GEM_INLINE void archive_search_group_clear(archive_search_group_t* const archive_search_group);
GEM_INLINE void archive_search_group_delete(archive_search_group_t* const archive_search_group);

GEM_INLINE archive_search_t* archive_search_group_alloc(archive_search_group_t* const archive_search_group);

GEM_INLINE bool archive_search_group_is_empty(archive_search_group_t* const archive_search_group);
GEM_INLINE bool archive_search_group_fits_in_buffer(
    archive_search_group_t* const archive_search_group,archive_search_t* const archive_search);

GEM_INLINE bool archive_search_group_add_search(
    archive_search_group_t* const archive_search_group,archive_search_t* const archive_search);
GEM_INLINE void archive_search_group_retrieve_begin(archive_search_group_t* const archive_search_group);
GEM_INLINE bool archive_search_group_get_search(
    archive_search_group_t* const archive_search_group,
    archive_search_t** const archive_search,bpm_gpu_buffer_t** const bpm_gpu_buffer);

/*
 * Error Messages
 */
#define GEM_ERROR_ARCHIVE_SEARCH_GROUP_QUERY_TOO_BIG "Archive-Search group. Couldn't copy query to BPM-buffer (Query too big)"

#endif /* ARCHIVE_SEARCH_GROUP_H_ */
