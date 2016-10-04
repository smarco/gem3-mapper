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

#include "filtering/candidates/filtering_candidates_buffered.h"

/*
 * Setup
 */
void filtering_candidates_buffered_clear(
    filtering_candidates_buffered_t* const filtering_candidates_buffered) {
  filtering_candidates_buffered->positions_buffered = NULL;
  filtering_candidates_buffered->num_positions = 0;
  filtering_candidates_buffered->regions_buffered = NULL;
  filtering_candidates_buffered->num_regions = 0;
}
/*
 * Allocators
 */
void filtering_candidates_buffered_allocate_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_filtering_positions) {
  filtering_candidates_buffered->positions_buffered =
      mm_stack_calloc(filtering_candidates->buffered_mm->mm_positions,
          num_filtering_positions,filtering_position_buffered_t,false);
  filtering_candidates_buffered->num_positions = num_filtering_positions;
}
void filtering_candidates_buffered_allocate_regions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t num_filtering_regions) {
  filtering_candidates_buffered->regions_buffered =
      mm_stack_calloc(filtering_candidates->buffered_mm->mm_regions,
          num_filtering_regions,filtering_region_buffered_t,false);
  filtering_candidates_buffered->num_regions = num_filtering_regions;
}
