/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_cache.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_REGION_CACHE_H_
#define FILTERING_REGION_CACHE_H_

#include "filtering_region.h"

typedef struct {
  /* Footprint Checksum */
  ihash_t* footprint_hash; // Footprint Checksum (filtering_region_t*)
} filtering_region_cache_t;

/*
 * Setup
 */
void filtering_region_cache_init(filtering_region_cache_t* const filtering_region_cache);
void filtering_region_cache_clear(filtering_region_cache_t* const filtering_region_cache);
void filtering_region_cache_destroy(filtering_region_cache_t* const filtering_region_cache);

/*
 * Accessors
 */
void filtering_region_cache_compute_footprint(
    filtering_region_t* const filtering_region,text_collection_t* const text_collection);
void filtering_region_cache_add(
    filtering_region_cache_t* const filtering_region_cache,filtering_region_t* const filtering_region,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack);

/*
 * Search
 */
match_trace_t* filtering_region_cache_search(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,text_collection_t* const text_collection);

#endif /* FILTERING_REGION_CACHE_H_ */
