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
  filtering_region_cache->footprint_hash = ihash_new();
  filtering_region_cache->last_aligned.filtering_region = NULL;
}
void filtering_region_cache_clear(filtering_region_cache_t* const filtering_region_cache) {
  ihash_clear(filtering_region_cache->footprint_hash);
  filtering_region_cache->last_aligned.filtering_region = NULL;
}
void filtering_region_cache_destroy(filtering_region_cache_t* const filtering_region_cache) {
  ihash_delete(filtering_region_cache->footprint_hash);
}
/*
 * Accessors
 */
void filtering_region_cache_compute_footprint(
    filtering_region_t* const filtering_region,
    text_collection_t* const text_collection) {
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
bool filtering_region_transient_cache_is_empty(
    filtering_region_cache_t* const filtering_region_cache) {
  return filtering_region_cache->last_aligned.filtering_region==NULL;
}
void filtering_region_transient_cache_add(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    uint64_t* const match_trace_offset,
    mm_stack_t* const mm_stack) {
  filtering_region_cache->last_aligned.filtering_region = filtering_region;
  filtering_region_cache->last_aligned.match_trace_offset = match_trace_offset;
}
void filtering_region_permanent_cache_add(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    uint64_t* const match_trace_offset,
    mm_stack_t* const mm_stack) {
  filtering_region_cache_element_t* const cache_element = mm_stack_alloc(mm_stack,filtering_region_cache_element_t);
  cache_element->filtering_region = filtering_region;
  cache_element->match_trace_offset = match_trace_offset;
  ihash_insert(filtering_region_cache->footprint_hash,filtering_region->footprint,cache_element);
}
/*
 * Search
 */
bool filtering_region_cache_cmp_regions(
    filtering_region_t* const filtering_region_a,
    filtering_region_t* const filtering_region_b,
    text_collection_t* const text_collection) {
  // Cmp text-lengths
  // const uint64_t text_length_a = filtering_region_a->end_position-filtering_region_a->begin_position;
  // const uint64_t text_length_b = filtering_region_b->end_position-filtering_region_b->begin_position;
  // if (text_length_a != text_length_b) return false;
  // Cmp read
  const text_trace_t* const text_trace_a = text_collection_get_trace(text_collection,filtering_region_a->text_trace_offset);
  const text_trace_t* const text_trace_b = text_collection_get_trace(text_collection,filtering_region_b->text_trace_offset);
  if (text_trace_a->regular_text_length != text_trace_b->regular_text_length) return false;
  return (memcmp(text_trace_a->regular_text,text_trace_b->regular_text,text_trace_a->regular_text_length)==0);
}
match_trace_t* filtering_region_transient_cache_search(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    text_collection_t* const text_collection,
    matches_t* const matches) {
  PROFILE_START(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
  // Basic compare (null & distance)
  filtering_region_cache_element_t* const last_aligned = &filtering_region_cache->last_aligned;
  if (last_aligned->filtering_region==NULL ||
      last_aligned->filtering_region->align_distance != filtering_region->align_distance) {
    PROFILE_STOP(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
    return NULL;
  }
  // Compare filtering regions
  if (filtering_region_cache_cmp_regions(filtering_region,last_aligned->filtering_region,text_collection)) {
    PROF_INC_COUNTER(GP_FC_CACHE_SEARCH_HIT);
    PROFILE_STOP(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
    return matches_get_match_trace(matches,*(last_aligned->match_trace_offset));
  } else {
    PROFILE_STOP(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
    return NULL;
  }
}
match_trace_t* filtering_region_permanent_cache_search(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region,
    text_collection_t* const text_collection,
    matches_t* const matches) {
  PROFILE_START(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
  // Search for filtering-region footprint (checksum)
  filtering_region_cache_element_t* const cache_element =
      ihash_get(filtering_region_cache->footprint_hash,filtering_region->footprint,filtering_region_cache_element_t);
  if (cache_element==NULL) {
    PROFILE_STOP(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
    return NULL;
  }
  // Compare filtering regions
  if (filtering_region_cache_cmp_regions(filtering_region,cache_element->filtering_region,text_collection)) {
    PROF_INC_COUNTER(GP_FC_CACHE_SEARCH_HIT);
    PROFILE_STOP(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
    return matches_get_match_trace(matches,*(cache_element->match_trace_offset));
  } else {
    PROFILE_STOP(GP_FC_CACHE_SEARCH,PROFILE_LEVEL);
    return NULL;
  }
}


