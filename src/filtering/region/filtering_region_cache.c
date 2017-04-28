/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Filtering region module provides a cache of already processed
 *   alignments as to avoid the overhead of (re)computing identical
 *   alignments
 */

#include "filtering/region/filtering_region_cache.h"

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
    filtering_region_t* const filtering_region_b) {
  // Compare lengths
  const text_trace_t* const text_trace_a = &filtering_region_a->text_trace;
  const text_trace_t* const text_trace_b = &filtering_region_b->text_trace;
  if (text_trace_a->text_length != text_trace_b->text_length) return false;
  // Compare strand (due to left-gap alignment, same text on diff strand can give diff alignments)
  if (text_trace_a->strand != text_trace_b->strand) return false;
  // Compare read // TODO Shrink to text-boundaries offsets
  return (memcmp(text_trace_a->text,text_trace_b->text,text_trace_a->text_length)==0);
}
match_trace_t* filtering_region_transient_cache_search(
    filtering_region_cache_t* const filtering_region_cache,
    filtering_region_t* const filtering_region) {
  // Basic compare (null & distance)
  if (filtering_region_cache->filtering_region==NULL) return NULL;
  const uint64_t cache_distance_min_bound =
      filtering_region_cache->filtering_region->alignment.distance_min_bound;
  const uint64_t region_distance_min_bound = filtering_region->alignment.distance_min_bound;
  if (cache_distance_min_bound != region_distance_min_bound) return NULL;
  // Compare filtering regions
  if (filtering_region_cache_cmp_regions(
      filtering_region,filtering_region_cache->filtering_region)) {
    PROF_INC_COUNTER(GP_FC_CACHE_SEARCH_HIT);
    return filtering_region_cache->match_trace;
  } else {
    return NULL;
  }
}

