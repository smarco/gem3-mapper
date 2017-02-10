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
 *   Filtering candidates module provides functions to process all
 *   bwt-encoded positions (originated from a candidate generation process)
 *   and compose them into filtering-regions (using decoded text coordinates)
 */

#ifndef FILTERING_CANDIDATES_PROCESS_H_
#define FILTERING_CANDIDATES_PROCESS_H_

#include "filtering/candidates/filtering_candidates.h"
#include "align/pattern/pattern.h"

/*
 * Constants
 */
#define DECODE_NUM_POSITIONS_PREFETCHED          10

/*
 * Batch decode
 */
typedef struct {
  uint64_t vector_rank;
  uint64_t index_position;
  uint64_t distance;
  uint64_t used_slot;
  bwt_block_locator_t bwt_block_locator;
} fc_batch_decode_candidate;

/*
 * Adjust the filtering-position and compute the coordinates or the candidate text
 */
void filtering_candidates_compute_text_coordinates(
    filtering_candidates_t* const filtering_candidates,
    filtering_position_t* const filtering_position,
    pattern_t* const pattern);

/*
 * Compose filtering regions
 */
uint64_t filtering_candidates_compose_filtering_regions(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool compose_alignment_regions);

/*
 * Process Candidates
 */
uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool compose_alignment_regions);

#endif /* FILTERING_CANDIDATES_PROCESS_H_ */
