/*
 * PROJECT: GEMMapper
 * FILE: mm_search.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *
 */

#include "mm_search.h"

/*
 * MM Search
 */
GEM_INLINE mm_search_t* mm_search_new(mm_slab_t* const mm_slab) {
  // Allocate
  mm_search_t* const mm_search = mm_alloc(mm_search_t);
  // MM-Stack
  mm_search->mm_stack = mm_stack_new(mm_slab);
  // Text-Collection Buffer
  text_collection_init(&mm_search->text_collection);
  // Filtering Candidates Buffers
  filtering_candidates_init(&mm_search->filtering_candidates_forward);
  filtering_candidates_init(&mm_search->filtering_candidates_reverse);
  // Interval Set
  interval_set_init(&mm_search->interval_set);
  // Return
  return mm_search;
}
GEM_INLINE void mm_search_clear(mm_search_t* const mm_search) {
  mm_stack_free(mm_search->mm_stack);
  text_collection_clear(&mm_search->text_collection);
  filtering_candidates_clear(&mm_search->filtering_candidates_forward);
  filtering_candidates_clear(&mm_search->filtering_candidates_reverse);
  interval_set_clear(&mm_search->interval_set);
}
GEM_INLINE void mm_search_delete(mm_search_t* const mm_search) {
  mm_stack_delete(mm_search->mm_stack);
  text_collection_destroy(&mm_search->text_collection);
  filtering_candidates_destroy(&mm_search->filtering_candidates_forward);
  filtering_candidates_destroy(&mm_search->filtering_candidates_reverse);
  interval_set_destroy(&mm_search->interval_set);
  mm_free(mm_search);
}
