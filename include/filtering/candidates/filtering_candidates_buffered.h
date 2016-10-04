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
 *   Filtering module provides a simple data structure to store
 *   filtering-regions & filtering-positions for batch processing
 */

#ifndef FILTERING_CANDIDATES_BUFFERED_H_
#define FILTERING_CANDIDATES_BUFFERED_H_

#include "utils/essentials.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/candidates/filtering_candidates_mm.h"
#include "matches/align/match_alignment.h"

/*
 * Filtering position/region buffered
 */
typedef struct {
  // Source Region
  uint64_t source_region_begin;        // Source-region Begin
  uint64_t source_region_end;          // Source-region End
} filtering_position_buffered_t;
typedef struct {
  /* Source Region Offset */
  uint32_t text_source_region_offset;          // Text-Offset to the begin of the source-region
  uint32_t key_source_region_offset;           // Key-Offset to the begin of the source-region
  /* Text */
  uint64_t text_begin_position;                // Region effective begin position (adjusted to error boundaries)
  uint64_t text_end_position;                  // Region effective end position (adjusted to error boundaries)
  /* Alignment */
  alignment_t alignment;                       // Alignment
  match_alignment_region_t* alignment_regions; // Match alignment-regions
  uint64_t num_alignment_regions;
  uint64_t scaffold_coverage;
} filtering_region_buffered_t;

/*
 * Filtering Candidates Buffered
 */
typedef struct {
  // Filtering Positions
  filtering_position_buffered_t* positions_buffered;
  uint64_t num_positions;
  // Filtering Regions
  filtering_region_buffered_t* regions_buffered;
  uint64_t num_regions;
} filtering_candidates_buffered_t;

/*
 * Setup
 */
void filtering_candidates_buffered_clear(
    filtering_candidates_buffered_t* const filtering_candidates_buffered);

/*
 * Allocators
 */
void filtering_candidates_buffered_allocate_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_filtering_positions);
void filtering_candidates_buffered_allocate_regions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_filtering_regions);


#endif /* FILTERING_CANDIDATES_BUFFERED_H_ */
