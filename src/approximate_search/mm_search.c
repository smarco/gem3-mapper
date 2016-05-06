/*
 * PROJECT: GEMMapper
 * FILE: mm_search.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *
 */

#include "approximate_search/mm_search.h"

/*
 * MM Search
 */
mm_search_t* mm_search_new(mm_slab_t* const mm_slab) {
  // Allocate
  mm_search_t* const mm_search = mm_alloc(mm_search_t);
  // Filtering candidates
  filtering_candidates_init(&mm_search->filtering_candidates_forward_end1);
  filtering_candidates_init(&mm_search->filtering_candidates_reverse_end1);
  filtering_candidates_init(&mm_search->filtering_candidates_forward_end2);
  filtering_candidates_init(&mm_search->filtering_candidates_reverse_end2);
  // Text-Collection Buffer
  text_collection_init(&mm_search->text_collection);
  // Interval Set
  interval_set_init(&mm_search->interval_set);
  // Stats
  mm_search->mapper_stats = mapper_stats_new();
  // MM-Stack
  mm_search->mm_stack = mm_stack_new(mm_slab);
  // Return
  return mm_search;
}
void mm_search_clear(mm_search_t* const mm_search) {
  filtering_candidates_clear(&mm_search->filtering_candidates_forward_end1);
  filtering_candidates_clear(&mm_search->filtering_candidates_reverse_end1);
  filtering_candidates_clear(&mm_search->filtering_candidates_forward_end2);
  filtering_candidates_clear(&mm_search->filtering_candidates_reverse_end2);
  text_collection_clear(&mm_search->text_collection);
  interval_set_clear(&mm_search->interval_set);
  mm_stack_clear(mm_search->mm_stack);
}
void mm_search_delete(mm_search_t* const mm_search) {
  filtering_candidates_destroy(&mm_search->filtering_candidates_forward_end1);
  filtering_candidates_destroy(&mm_search->filtering_candidates_reverse_end1);
  filtering_candidates_destroy(&mm_search->filtering_candidates_forward_end2);
  filtering_candidates_destroy(&mm_search->filtering_candidates_reverse_end2);
  text_collection_destroy(&mm_search->text_collection);
  interval_set_destroy(&mm_search->interval_set);
  mapper_stats_delete(mm_search->mapper_stats);
  mm_stack_delete(mm_search->mm_stack);
  mm_free(mm_search);
}
