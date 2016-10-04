/*
 * PROJECT: GEMMapper
 * FILE: search_pipeline_handlers.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_PIPELINE_HANDLERS_H_
#define SEARCH_PIPELINE_HANDLERS_H_

#include "utils/essentials.h"
#include "archive/archive.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_se_parameters.h"
#include "filtering/candidates/filtering_candidates_mm.h"
#include "mapper/mapper_stats.h"

/*
 * Search Pipeline Handlers
 */
typedef struct {
  /* Mapping Statistics */
  mapper_stats_t* mapper_stats;         // Mapping Statistics
  /* Filtering Candidates */
  filtering_candidates_t fc_decode_end1;
  filtering_candidates_t fc_decode_end2;
  filtering_candidates_mm_t fc_decode_mm;
  filtering_candidates_t fc_verify_end1;
  filtering_candidates_t fc_verify_end2;
  filtering_candidates_mm_t fc_verify_mm;
  filtering_candidates_buffered_mm_t fc_buffered_mm;
  /* Neighborhood Search */
  nsearch_schedule_t nsearch_schedule;
  /* MM */
  mm_slab_t* mm_slab;                   // MM-Slab
  mm_stack_t* mm_stack;                 // MM-Stack
} search_pipeline_handlers_t;

/*
 * Setup
 */
search_pipeline_handlers_t* search_pipeline_handlers_new(void);
void search_pipeline_handlers_clear(search_pipeline_handlers_t* const search_pipeline_handlers);
void search_pipeline_handlers_clear_positions(search_pipeline_handlers_t* const search_pipeline_handlers);
void search_pipeline_handlers_clear_regions(search_pipeline_handlers_t* const search_pipeline_handlers);
void search_pipeline_handlers_delete(search_pipeline_handlers_t* const search_pipeline_handlers);

/*
 * Injection (Support Data Structures)
 */
void search_pipeline_handlers_inject_se(
    archive_search_t* const archive_search,
    search_pipeline_handlers_t* const search_pipeline_handlers);
void search_pipeline_handlers_inject_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    search_pipeline_handlers_t* const search_pipeline_handlers);

#endif /* SEARCH_PIPELINE_HANDLERS_H_ */
