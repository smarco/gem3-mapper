/*
 * PROJECT: GEMMapper
 * FILE: interval_set.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef INTERVAL_SET_H_
#define INTERVAL_SET_H_

#include "utils/essentials.h"

typedef struct {
  /* BWT Interval */
  uint64_t lo;
  uint64_t hi;
  /* Matching Properties */
  uint64_t distance;
  uint64_t length;
} interval_t;

typedef struct {
  vector_t* intervals; // (interval_t*)
} interval_set_t;

/*
 * Setup
 */
void interval_set_init(interval_set_t* const interval_set);
void interval_set_clear(interval_set_t* const interval_set);
void interval_set_destroy(interval_set_t* const interval_set);

/*
 * Counting
 */
uint64_t interval_set_count_intervals(interval_set_t* const interval_set);
uint64_t interval_set_count_intervals_length(interval_set_t* const interval_set);
uint64_t interval_set_count_intervals_length_thresholded(
    interval_set_t* const interval_set,
    const uint64_t max_error);

/*
 * Adding
 */
void interval_set_add(
    interval_set_t* const interval_set,
    const uint64_t lo,
    const uint64_t hi,
    const uint64_t distance,
    const uint64_t length);

/*
 * Set Operators
 */
void interval_set_union(
    interval_set_t* const interval_set_a,
    interval_set_t* const interval_set_b);
void interval_set_subtract(
    interval_set_t* const result_set,
    interval_set_t* const exclusion_set);

/*
 * Macro iterator
 */
#define INTERVAL_SET_ITERATE(interval_set,element) \
  VECTOR_ITERATE(interval_set->intervals,element,interval_set_##position,interval_t)


#endif /* INTERVAL_SET_H_ */
