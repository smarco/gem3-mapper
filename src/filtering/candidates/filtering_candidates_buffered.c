/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 * DESCRIPTION:
 *   Filtering module provides a simple data structure to store
 *   filtering-regions & filtering-positions for batch processing
 */

#include "filtering/candidates/filtering_candidates_buffered.h"

/*
 * Setup
 */
void filtering_candidates_buffered_inject_handlers(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    mm_allocator_t* const mm_allocator) {
  // Init structures
  filtering_candidates_buffered->positions = NULL;
  filtering_candidates_buffered->num_positions = 0;
  filtering_candidates_buffered->regions = NULL;
  filtering_candidates_buffered->num_regions = 0;
  filtering_candidates_buffered->num_canonical_regions = 0;
  filtering_candidates_buffered->discarded_regions = NULL;
  filtering_candidates_buffered->num_discarded_regions = 0;
  // MM
  filtering_candidates_buffered->mm_allocator = mm_allocator;
}
/*
 * Allocators
 */
void filtering_candidates_buffered_allocate_positions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_filtering_positions) {
  filtering_candidates_buffered->positions = mm_allocator_calloc(
      filtering_candidates_buffered->mm_allocator,num_filtering_positions,filtering_position_t*,false);
}
void filtering_candidates_buffered_free_positions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered) {
  mm_allocator_free(filtering_candidates_buffered->mm_allocator,filtering_candidates_buffered->positions);
  filtering_candidates_buffered->positions = NULL;
}
void filtering_candidates_buffered_allocate_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_filtering_regions) {
  filtering_candidates_buffered->regions = mm_allocator_calloc(
      filtering_candidates_buffered->mm_allocator,num_filtering_regions,filtering_region_t*,false);
}
void filtering_candidates_buffered_free_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered) {
  mm_allocator_free(filtering_candidates_buffered->mm_allocator,filtering_candidates_buffered->regions);
  filtering_candidates_buffered->regions = NULL;
}
void filtering_candidates_buffered_allocate_discarded_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_discarded_regions) {
  filtering_candidates_buffered->discarded_regions = mm_allocator_calloc(
      filtering_candidates_buffered->mm_allocator,num_discarded_regions,filtering_region_t*,false);
}
void filtering_candidates_buffered_free_discarded_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered) {
  mm_allocator_free(filtering_candidates_buffered->mm_allocator,filtering_candidates_buffered->discarded_regions);
  filtering_candidates_buffered->discarded_regions = NULL;
}
/*
 * Accessors
 */
uint64_t filtering_candidates_buffered_get_num_regions(
    const filtering_candidates_buffered_t* const filtering_candidates_buffered) {
  return filtering_candidates_buffered->num_regions;
}
uint64_t filtering_candidates_buffered_get_num_canonical_regions(
    const filtering_candidates_buffered_t* const filtering_candidates_buffered) {
  return filtering_candidates_buffered->num_canonical_regions;
}
uint64_t filtering_candidates_buffered_get_region_tiles_length(
    const filtering_candidates_buffered_t* const filtering_candidates_buffered) {
  uint64_t i, region_tiles_length = 0;
  for (i=0;i<filtering_candidates_buffered->num_regions;++i) {
    filtering_region_t* const filtering_region = filtering_candidates_buffered->regions[i];
    region_tiles_length += filtering_region->text_end_position - filtering_region->text_begin_position;
  }
  return region_tiles_length;
}
