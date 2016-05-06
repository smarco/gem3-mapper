/*
 * PROJECT: GEMMapper
 * FILE: region_set.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef REGION_SET_H_
#define REGION_SET_H_

#include "utils/essentials.h"

/*
 * Region Set
 */
typedef struct {
  uint64_t begin_position;
  uint64_t end_position;
} region_interval_t;
typedef struct {
  vector_t* region_intervals; // (region_interval_t*)
} region_set_t;

/*
 * Setup
 */
void region_set_init(region_set_t* const region_set);
void region_set_clear(region_set_t* const region_set);
void region_set_destroy(region_set_t* const region_set);

/*
 * Adding
 */
void region_set_add(
    region_set_t* const region_set,
    const uint64_t begin_position,
    const uint64_t end_position);

/*
 * Sort
 */
void region_set_sort(region_set_t* const region_set);

/*
 * Set operators
 */
bool region_set_is_contained(
    region_set_t* const region_set,
    const uint64_t begin_position,
    const uint64_t end_position);

/*
 * Macro iterator
 */
#define REGION_SET_ITERATE(region_set,element) \
  VECTOR_ITERATE(region_set->region_intervals,element,region_set_##position,region_interval_t)


#endif /* REGION_SET_H_ */
