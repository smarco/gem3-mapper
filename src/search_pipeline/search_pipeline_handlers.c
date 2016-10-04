/*
 * PROJECT: GEMMapper
 * FILE: search_pipeline_handlers.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "search_pipeline/search_pipeline_handlers.h"

/*
 * Setup
 */
search_pipeline_handlers_t* search_pipeline_handlers_new(void) {
  // Allocate
  search_pipeline_handlers_t* const search_pipeline_handlers = mm_alloc(search_pipeline_handlers_t);
  // Filtering Candidates
  filtering_candidates_init(&search_pipeline_handlers->fc_decode_end1);
  filtering_candidates_init(&search_pipeline_handlers->fc_decode_end2);
  filtering_candidates_mm_init(&search_pipeline_handlers->fc_decode_mm);
  filtering_candidates_init(&search_pipeline_handlers->fc_verify_end1);
  filtering_candidates_init(&search_pipeline_handlers->fc_verify_end2);
  filtering_candidates_mm_init(&search_pipeline_handlers->fc_verify_mm);
  // Filtering Candidates Buffered
  filtering_candidates_buffered_mm_init(&search_pipeline_handlers->fc_buffered_mm);
  filtering_candidates_mm_inject_buffered_mm(
      &search_pipeline_handlers->fc_verify_mm,
      &search_pipeline_handlers->fc_buffered_mm);
  filtering_candidates_mm_inject_buffered_mm(
      &search_pipeline_handlers->fc_decode_mm,
      &search_pipeline_handlers->fc_buffered_mm);
  // Stats
  search_pipeline_handlers->mapper_stats = mapper_stats_new();
  // MM
  search_pipeline_handlers->mm_slab = mm_slab_new_(
      BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM,"archive_search_handlers.8MB");
  search_pipeline_handlers->mm_stack = mm_stack_new(search_pipeline_handlers->mm_slab);
  // Return
  return search_pipeline_handlers;
}
void search_pipeline_handlers_clear(search_pipeline_handlers_t* const search_pipeline_handlers) {
  filtering_candidates_mm_clear(&search_pipeline_handlers->fc_decode_mm);
  filtering_candidates_mm_clear(&search_pipeline_handlers->fc_verify_mm);
  filtering_candidates_buffered_mm_clear(&search_pipeline_handlers->fc_buffered_mm);
  mm_stack_clear(search_pipeline_handlers->mm_stack);
}
void search_pipeline_handlers_delete(search_pipeline_handlers_t* const search_pipeline_handlers) {
  filtering_candidates_buffered_mm_destroy(&search_pipeline_handlers->fc_buffered_mm);
  filtering_candidates_mm_destroy(&search_pipeline_handlers->fc_decode_mm);
  filtering_candidates_mm_destroy(&search_pipeline_handlers->fc_verify_mm);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_decode_end1);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_decode_end2);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_verify_end1);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_verify_end2);
  mapper_stats_delete(search_pipeline_handlers->mapper_stats);
  mm_stack_delete(search_pipeline_handlers->mm_stack);
  mm_slab_delete(search_pipeline_handlers->mm_slab);
  mm_free(search_pipeline_handlers);
}
/*
 * Injection (Support Data Structures)
 */
void search_pipeline_handlers_inject_se(
    archive_search_t* const archive_search,
    search_pipeline_handlers_t* const search_pipeline_handlers) {
  // Archive Search
  archive_search->mapper_stats = search_pipeline_handlers->mapper_stats;
  archive_search->mm_stack = search_pipeline_handlers->mm_stack;
  // Approximate Search
  approximate_search_t* const search = &archive_search->approximate_search;
  region_profile_inject_mm(&search->region_profile,search_pipeline_handlers->mm_stack);
  search->nsearch_schedule = &search_pipeline_handlers->nsearch_schedule;
  nsearch_schedule_inject_mm(search->nsearch_schedule,search_pipeline_handlers->mm_stack);
}
void search_pipeline_handlers_inject_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    search_pipeline_handlers_t* const search_pipeline_handlers) {
  // Archive Search
  archive_search_end1->mapper_stats = search_pipeline_handlers->mapper_stats;
  archive_search_end2->mapper_stats = search_pipeline_handlers->mapper_stats;
  archive_search_end1->mm_stack = search_pipeline_handlers->mm_stack;
  archive_search_end2->mm_stack = search_pipeline_handlers->mm_stack;
  // Approximate Search
  approximate_search_t* const search_end1 = &archive_search_end1->approximate_search;
  approximate_search_t* const search_end2 = &archive_search_end2->approximate_search;
  region_profile_inject_mm(&search_end1->region_profile,search_pipeline_handlers->mm_stack);
  region_profile_inject_mm(&search_end2->region_profile,search_pipeline_handlers->mm_stack);
  search_end1->nsearch_schedule = &search_pipeline_handlers->nsearch_schedule;
  nsearch_schedule_inject_mm(search_end1->nsearch_schedule,search_pipeline_handlers->mm_stack);
  search_end2->nsearch_schedule = &search_pipeline_handlers->nsearch_schedule;
  nsearch_schedule_inject_mm(search_end2->nsearch_schedule,search_pipeline_handlers->mm_stack);
}
