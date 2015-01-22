/*
 * PROJECT: GEMMapper
 * FILE: neighborhood_search.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NEIGHBORHOOD_SEARCH_H_
#define NEIGHBORHOOD_SEARCH_H_

#include "essentials.h"
#include "fm_index.h"
#include "interval_set.h"

/*
 * TLS related structures
 */
//typedef struct {
//  idx_t lo;
//  idx_t hi;
//  idx_t misms;
//} interval_t;
//typedef struct {
//  idx_t lo;
//  idx_t hi;
//} basic_interval_t;
//typedef struct {
//  idx_t begin;
//  idx_t count;
//  idx_t state;
//} group_intervals_t;

// GEM_INLINE uint64_t neighborhood_search();

GEM_INLINE void neighborhood_search(
    const fm_index_t* const fm_index,const uint8_t* const key,const uint64_t length,
    const uint64_t max_error,interval_set_t* const intervals_result,mm_stack_t* const mm_stack);

#endif /* NEIGHBORHOOD_SEARCH_H_ */
