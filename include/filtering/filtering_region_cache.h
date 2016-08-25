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
  match_trace_t* match_trace;
} filtering_region_cache_t;

/*
 * Setup
 */
void filtering_region_cache_init(filtering_region_cache_t* const filtering_region_cache);
void filtering_region_cache_clear(filtering_region_cache_t* const filtering_region_cache);

/*
 * Accessors
 */
bool filtering_region_transient_cache_is_empty(
    filtering_region_cache_t* const filtering_region_cache);
void filtering_region_transient_cache_add(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    match_trace_t* const match_trace);

/*
 * Search
 */
match_trace_t* filtering_region_transient_cache_search(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    text_collection_t* const text_collection);

#endif /* FILTERING_REGION_CACHE_H_ */
