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

#include "utils/essentials.h"
#include "data_structures/text_collection.h"
#include "filtering/filtering_candidates.h"
#include "mapper/mapper_stats.h"

/*
 * MM Search
 */
typedef struct {
  /* Filtering candidates */
  filtering_candidates_t filtering_candidates_forward_end1;
  filtering_candidates_t filtering_candidates_reverse_end1;
  filtering_candidates_t filtering_candidates_forward_end2;
  filtering_candidates_t filtering_candidates_reverse_end2;
  /* Text-Collection Buffer */
  text_collection_t text_collection;    // Stores text-traces (candidates/matches/regions/...)
  /* Interval Set */
  interval_set_t interval_set;          // Interval-Set
  /* Stats */
  mapper_stats_t* mapper_stats;         // Mapping Statistics
  /* MM-Stack */
  mm_stack_t* mm_stack;                 // Memory-Stack allocator
} mm_search_t;

/*
 * Setup
 */
mm_search_t* mm_search_new();
void mm_search_clear(mm_search_t* const mm_search);
void mm_search_delete(mm_search_t* const mm_search);

#endif /* MM_STACK_H_ */
