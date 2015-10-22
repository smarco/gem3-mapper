/*
 * PROJECT: GEMMapper
 * FILE: search_group.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_GROUP_H_
#define SEARCH_GROUP_H_

#include "essentials.h"
#include "align_bpm_gpu.h"
#include "search_group_verify_candidates.h"
#include "archive_search_cache.h"
#include "archive_search.h"
#include "mapper.h"

/*
 * Archive-Search groups
 */
typedef struct {
  /* Search-group buffer to verify candidates */
  search_group_verify_candidates_t* search_group_vc;
  /* Archive-search cache */
  archive_search_cache_t* archive_search_cache;
  /* MM */
  mm_search_t* mm_search;
} search_group_t;

/*
 * Setup
 */
search_group_t* search_group_new(
    mapper_parameters_t* const mapper_parameters,bpm_gpu_buffer_t* const bpm_gpu_buffers,
    const uint64_t num_verify_candidates_buffers);
void search_group_init(search_group_t* const search_group);
void search_group_clear(search_group_t* const search_group);
void search_group_delete(search_group_t* const search_group);
/*
 * Archive-Search allocation
 */
archive_search_t* search_group_allocate_se(search_group_t* const search_group);
void search_group_allocate_pe(
    search_group_t* const search_group,
    archive_search_t** const archive_search_end1,archive_search_t** const archive_search_end2);
/*
 * Accessors
 */
text_collection_t* search_group_get_text_collection(search_group_t* const search_group);

#endif /* SEARCH_GROUP_H_ */
