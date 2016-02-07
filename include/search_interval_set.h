/*
 * PROJECT: GEMMapper
 * FILE: search_interval_set.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef SEARCH_INTERVAL_SET_H_
#define SEARCH_INTERVAL_SET_H_

#include "essentials.h"

typedef struct {
  /* BWT Interval */
  uint64_t lo;
  uint64_t hi;
  /* Matching Properties */
  uint64_t distance;
  uint64_t length;
} search_interval_t;

typedef struct {
  vector_t* intervals; // (interval_t*)
} search_interval_set_t;

/*
 * Setup
 */
void search_interval_set_init(search_interval_set_t* const search_interval_set);
void search_interval_set_clear(search_interval_set_t* const search_interval_set);
void search_interval_set_destroy(search_interval_set_t* const search_interval_set);

/*
 * Counting
 */
uint64_t search_interval_set_count_intervals(search_interval_set_t* const search_interval_set);
uint64_t search_interval_set_count_intervals_length(search_interval_set_t* const search_interval_set);
uint64_t search_interval_set_count_intervals_length_thresholded(
    search_interval_set_t* const search_interval_set,const uint64_t max_error);

/*
 * Adding
 */
void search_interval_set_add(
    search_interval_set_t* const search_interval_set,const uint64_t lo,const uint64_t hi,
    const uint64_t distance,const uint64_t length);

/*
 * Set Operators
 */
void search_interval_set_union(search_interval_set_t* const search_interval_set_a,search_interval_set_t* const search_interval_set_b);
void search_interval_set_subtract(search_interval_set_t* const result_set,search_interval_set_t* const exclusion_set);

/*
 * Macro iterator
 */
#define SEARCH_INTERVAL_SET_ITERATE(search_interval_set,element) \
  VECTOR_ITERATE(search_interval_set->intervals,element,search_interval_set_##position,search_interval_t)

#endif /* SEARCH_INTERVAL_SET_H_ */
