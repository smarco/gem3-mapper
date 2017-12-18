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

#ifndef FILTERING_CANDIDATES_BUFFERED_H_
#define FILTERING_CANDIDATES_BUFFERED_H_

#include "utils/essentials.h"
#include "filtering/candidates/filtering_candidates.h"
#include "matches/align/match_alignment.h"

/*
 * Filtering Candidates Buffered
 */
typedef struct {
  /* Filtering Positions */
  filtering_position_t** positions;
  uint64_t num_positions;
  /* Filtering Regions */
  filtering_region_t** regions;
  uint64_t num_regions;
  /* Canonical: Number of regions processed by GPU */
  uint64_t num_canonical_regions;
  /* Discarded Regions */
  filtering_region_t** discarded_regions;
  uint64_t num_discarded_regions;
  /* MM */
  mm_allocator_t* mm_allocator;
} filtering_candidates_buffered_t;

/*
 * Setup
 */
void filtering_candidates_buffered_inject_handlers(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    mm_allocator_t* const mm_allocator);

/*
 * Allocators
 */
void filtering_candidates_buffered_allocate_positions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_filtering_positions);
void filtering_candidates_buffered_free_positions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered);

void filtering_candidates_buffered_allocate_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_filtering_regions);
void filtering_candidates_buffered_free_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered);

void filtering_candidates_buffered_allocate_discarded_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_discarded_regions);
void filtering_candidates_buffered_free_discarded_regions(
    filtering_candidates_buffered_t* const filtering_candidates_buffered);

/*
 * Accessors
 */
uint64_t filtering_candidates_buffered_get_num_regions(
    const filtering_candidates_buffered_t* const filtering_candidates_buffered);
uint64_t filtering_candidates_buffered_get_num_canonical_regions(
    const filtering_candidates_buffered_t* const filtering_candidates_buffered);
uint64_t filtering_candidates_buffered_get_region_tiles_length(
    const filtering_candidates_buffered_t* const filtering_candidates_buffered);

#endif /* FILTERING_CANDIDATES_BUFFERED_H_ */
