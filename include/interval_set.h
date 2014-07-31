/*
 * PROJECT: GEMMapper
 * FILE: interval_set.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef INTERVAL_SET_H_
#define INTERVAL_SET_H_

#include "essentials.h"

typedef struct {
  /* HI/LO Interval */
  uint64_t lo;
  uint64_t hi;
  /* Levenshtein Error */
  uint64_t error;
} interval_t;

typedef struct {
  /* HI/LO Interval */
  uint64_t lo;
  uint64_t hi;
//  /* DP Column */
//  interval_search_t* prev_....
//  uint64_t* dp;
} interval_search_t;


typedef struct {
  svector_t* intervals; /* (interval_t) */

  // TODO
  svector_iterator_t pending_candidates_iterator; // Writing Iterator (always appending)
  /* MM */ // TODO
  mm_slab_t* mm_slab;
  mm_stack_t* mm_stack;
} interval_set_t;

GEM_INLINE void interval_set_new(interval_set_t* const interval_set);
GEM_INLINE void interval_set_clear(interval_set_t* const interval_set);
GEM_INLINE void interval_set_delete(interval_set_t* const interval_set);

GEM_INLINE uint64_t interval_set_count_intervals(interval_set_t* const interval_set);
GEM_INLINE uint64_t interval_set_count_intervals_length(interval_set_t* const interval_set);
GEM_INLINE uint64_t interval_set_count_intervals_length_thresholded(
    interval_set_t* const interval_set,const uint64_t max_error);


#endif /* INTERVAL_SET_H_ */
