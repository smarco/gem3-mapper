/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_process.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_CANDIDATES_PROCESS_H_
#define FILTERING_CANDIDATES_PROCESS_H_

#include "filtering_candidates.h"
#include "archive_search_parameters.h"
#include "archive.h"
#include "pattern.h"

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
 * Retrieve all candidates(text) from the index
 */
void filtering_candidates_retrieve_filtering_regions(
    filtering_candidates_t* const filtering_candidates,archive_text_t* const archive_text,
    text_collection_t* const text_collection,mm_stack_t* const mm_stack);

/*
 * Filtering adjustment of the position wrt region/seed on which the candidate is based
 */
void filtering_candidates_adjust_filtering_position(
    filtering_position_t* const filtering_position,archive_text_t* const archive_text,
    const uint64_t begin_offset,const uint64_t end_offset,const uint64_t boundary_error);

/*
 * Compose filtering regions
 */
uint64_t filtering_candidates_compose_filtering_regions(
    filtering_candidates_t* const filtering_candidates,const uint64_t key_length,
    const uint64_t max_delta_difference,const bool compose_region_chaining,
    mm_stack_t* const mm_stack);

/*
 * Process Candidates
 */
uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    const bool compose_region_chaining,mm_stack_t* const mm_stack);

#endif /* FILTERING_CANDIDATES_PROCESS_H_ */
