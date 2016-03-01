/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_cache.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_REGION_CACHE_H_
#define FILTERING_REGION_CACHE_H_

#include "filtering/filtering_region.h"
#include "matches/matches.h"

typedef struct {
  filtering_region_t* filtering_region;
  uint64_t* match_trace_offset;
} filtering_region_cache_element_t;
typedef struct {
  /* Footprint Checksum */
  ihash_t* footprint_hash; // Footprint Checksum (filtering_region_t*)
  /* Last Region Aligned */
  filtering_region_cache_element_t last_aligned;
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
bool filtering_region_transient_cache_is_empty(
    filtering_region_cache_t* const filtering_region_cache);
void filtering_region_transient_cache_add(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    uint64_t* const match_trace_offset,
    mm_stack_t* const mm_stack);

/*
 * Search
 */
match_trace_t* filtering_region_transient_cache_search(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    text_collection_t* const text_collection,
    matches_t* const matches);

#endif /* FILTERING_REGION_CACHE_H_ */
