/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_cache.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_region_cache.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Constants
 */
#define FILTERING_REGION_CACHE_TEXT_FOOTPRINT_LENGTH 100

/*
 * Setup
 */
void filtering_region_cache_init(filtering_region_cache_t* const filtering_region_cache) {
  filtering_region_cache->filtering_region = NULL;
  filtering_region_cache->match_trace = NULL;
}
void filtering_region_cache_clear(filtering_region_cache_t* const filtering_region_cache) {
  filtering_region_cache->filtering_region = NULL;
  filtering_region_cache->match_trace = NULL;
}
/*
 * Accessors
 */
bool filtering_region_transient_cache_is_empty(
    filtering_region_cache_t* const filtering_region_cache) {
  return filtering_region_cache->filtering_region==NULL;
}
void filtering_region_transient_cache_add(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    match_trace_t* const match_trace) {
  filtering_region_cache->filtering_region = filtering_region;
  filtering_region_cache->match_trace = match_trace;
}
/*
 * Search
 */
bool filtering_region_cache_cmp_regions(
    filtering_region_t* const filtering_region_a,
    filtering_region_t* const filtering_region_b,
    text_collection_t* const text_collection) {
  // Cmp read
  const text_trace_t* const text_trace_a = text_collection_get_trace(text_collection,filtering_region_a->text_trace_offset);
  const text_trace_t* const text_trace_b = text_collection_get_trace(text_collection,filtering_region_b->text_trace_offset);
  if (text_trace_a->text_length != text_trace_b->text_length) return false;
  // TODO Shrink to text-boundaries offsets
  return (memcmp(text_trace_a->text,text_trace_b->text,text_trace_a->text_length)==0);
}
match_trace_t* filtering_region_transient_cache_search(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    text_collection_t* const text_collection) {
  // Basic compare (null & distance)
  if (filtering_region_cache->filtering_region==NULL) return NULL;
  const uint64_t cache_distance_min_bound = filtering_region_cache->filtering_region->region_alignment.distance_min_bound;
  const uint64_t region_distance_min_bound = filtering_region->region_alignment.distance_min_bound;
  if (cache_distance_min_bound != region_distance_min_bound) return NULL;
  // Compare filtering regions
  if (filtering_region_cache_cmp_regions(filtering_region,filtering_region_cache->filtering_region,text_collection)) {
    PROF_INC_COUNTER(GP_FC_CACHE_SEARCH_HIT);
    return filtering_region_cache->match_trace;
  } else {
    return NULL;
  }
}

