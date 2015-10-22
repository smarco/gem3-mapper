/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_cache.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_region_cache.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Constants
 */
#define FILTERING_REGION_CACHE_TEXT_FOOTPRINT_LENGTH 100

/*
 * Cached
 */
typedef struct {
  filtering_region_t* filtering_region;
  match_trace_t *match_trace;
} filtering_region_cache_element_t;

/*
 * Setup
 */
GEM_INLINE void filtering_region_cache_init(filtering_region_cache_t* const filtering_region_cache) {
  filtering_region_cache->footprint_hash = ihash_new();
}
GEM_INLINE void filtering_region_cache_clear(filtering_region_cache_t* const filtering_region_cache) {
  ihash_clear(filtering_region_cache->footprint_hash);
}
GEM_INLINE void filtering_region_cache_destroy(filtering_region_cache_t* const filtering_region_cache) {
  ihash_delete(filtering_region_cache->footprint_hash);
}
/*
 * Accessors
 */
GEM_INLINE void filtering_region_cache_compute_footprint(
    filtering_region_t* const filtering_region,text_collection_t* const text_collection) {
  PROFILE_START(GP_FC_CACHE_COMPUTE_FOOTPRINT,PROFILE_LEVEL);
  const text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  const uint64_t text_length = filtering_region->end_position-filtering_region->begin_position;
  const uint8_t* const text = text_trace->text;
  const uint64_t footprint_text_length = MIN(text_length,FILTERING_REGION_CACHE_TEXT_FOOTPRINT_LENGTH);
  // Compute footprint as a sort of circular checksum
  uint64_t footprint = 0, i;
  for (i=0;i<footprint_text_length;++i) {
    footprint ^= ((uint64_t)text[i] << ((i%16)*4));
  }
  // Set footprint
  filtering_region->footprint = footprint;
  PROFILE_STOP(GP_FC_CACHE_COMPUTE_FOOTPRINT,PROFILE_LEVEL);
}
GEM_INLINE void filtering_region_cache_add(
    filtering_region_cache_t* const filtering_region_cache,filtering_region_t* const filtering_region,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack) {
  filtering_region_cache_element_t* const cache_element = mm_stack_alloc(mm_stack,filtering_region_cache_element_t);
  cache_element->filtering_region = filtering_region;
  cache_element->match_trace = match_trace;
  ihash_insert(filtering_region_cache->footprint_hash,filtering_region->footprint,cache_element);
}
/*
 * Search
 */
GEM_INLINE bool filtering_region_cache_cmp_regions(
    filtering_region_t* const filtering_region_a,filtering_region_t* const filtering_region_b,
    text_collection_t* const text_collection) {
  // Cmp text-lengths
  const uint64_t text_length_a = filtering_region_a->end_position-filtering_region_a->begin_position;
  const uint64_t text_length_b = filtering_region_b->end_position-filtering_region_b->begin_position;
  if (text_length_a != text_length_b) return false;
  // Cmp read
  const text_trace_t* const text_trace_a = text_collection_get_trace(text_collection,filtering_region_a->text_trace_offset);
  const uint8_t* const text_a = text_trace_a->text;
  const text_trace_t* const text_trace_b = text_collection_get_trace(text_collection,filtering_region_b->text_trace_offset);
  const uint8_t* const text_b = text_trace_b->text;
  return (memcmp(text_a,text_b,text_length_b)==0);
}
GEM_INLINE match_trace_t* filtering_region_cache_search(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,text_collection_t* const text_collection) {
  PROFILE_START(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
  // Search for filtering-region footprint (checksum)
  filtering_region_cache_element_t* const cache_element =
      ihash_get(filtering_region_cache->footprint_hash,filtering_region->footprint,filtering_region_cache_element_t);
  if (cache_element==NULL) return NULL;
  // Compare filtering regions
  if (filtering_region_cache_cmp_regions(filtering_region,cache_element->filtering_region,text_collection)) {
    PROF_INC_COUNTER(GP_FC_CACHE_SEARCH_HIT);
    PROFILE_STOP(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
    return cache_element->match_trace;
  } else {
    PROFILE_STOP(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
    return NULL;
  }
}

