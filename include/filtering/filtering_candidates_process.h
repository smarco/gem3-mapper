/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_process.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_CANDIDATES_PROCESS_H_
#define FILTERING_CANDIDATES_PROCESS_H_

#include "filtering/filtering_candidates.h"
#include "data_structures/pattern.h"

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
    pattern_t* const pattern,
    const uint64_t begin_offset,
    const uint64_t end_offset);

/*
 * Compose filtering regions
 */
uint64_t filtering_candidates_compose_filtering_regions(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t key_length,
    const uint64_t max_delta_difference,
    const bool compose_region_chaining);

/*
 * Process Candidates
 */
uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool keep_matching_regions);

#endif /* FILTERING_CANDIDATES_PROCESS_H_ */
