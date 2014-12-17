/*
 * PROJECT: GEMMapper
 * FILE: mm_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *
 */

#ifndef MM_SEARCH_H_
#define MM_SEARCH_H_

#include "essentials.h"
#include "text_collection.h"
#include "filtering_candidates.h"

/*
 * Checkers
 */
#define MM_SEARCH_CHECK(mm_search) GEM_CHECK_NULL(mm_search)

/*
 * MM Search
 */
typedef struct {
  /* MM-Stack */
  mm_stack_t* mm_stack;               // Memory-Stack allocator
  /* Text-Collection Buffer */
  text_collection_t text_collection;  // Stores text-traces (candidates/matches/regions/...)
  /* Interval Set */
  interval_set_t interval_set;
} mm_search_t;

/*
 * Setup
 */
GEM_INLINE mm_search_t* mm_search_new();
GEM_INLINE void mm_search_clear(mm_search_t* const mm_search);
GEM_INLINE void mm_search_delete(mm_search_t* const mm_search);

/*
 * Errors
 */
//#define GEM_ERROR_MM_SEARCH_ ""

#endif /* MM_STACK_H_ */
