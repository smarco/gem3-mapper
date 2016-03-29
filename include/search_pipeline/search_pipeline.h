/*
 * PROJECT: GEMMapper
 * FILE: search_pipeline.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_PIPELINE_H_
#define SEARCH_PIPELINE_H_

#include "utils/essentials.h"
#include "search_pipeline/search_stage_region_profile.h"
#include "search_pipeline/search_stage_decode_candidates.h"
#include "search_pipeline/search_stage_verify_candidates.h"
#include "archive/archive_search.h"
#include "archive/archive_search_cache.h"
#include "mapper/mapper.h"

/*
 * Archive-Search groups
 */
typedef struct {
  /* Search-Stages buffer to verify candidates */
  search_stage_region_profile_t* stage_region_profile;
  search_stage_decode_candidates_t* stage_decode_candidates;
  search_stage_verify_candidates_t* stage_verify_candidates;
  /* Archive-search cache */
  archive_search_cache_t* archive_search_cache;
  /* Support Data Structures */
  mm_stack_t* mm_stack;                 // Memory-Stack allocator
  mapper_stats_t* mapper_stats;         // Mapping Statistics
  interval_set_t interval_set;          // Interval-Set
} search_pipeline_t;

/*
 * Setup
 */
search_pipeline_t* search_pipeline_new(
    mapper_parameters_t* const mapper_parameters,
    gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const bool paired_end);
void search_pipeline_clear(search_pipeline_t* const search_pipeline);
void search_pipeline_delete(search_pipeline_t* const search_pipeline);

/*
 * Archive-Search allocation
 */
void search_pipeline_allocate_se(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search);
void search_pipeline_allocate_pe(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);

#endif /* SEARCH_PIPELINE_H_ */
