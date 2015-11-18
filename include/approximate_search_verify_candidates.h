/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_verify_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_VERIFY_CANDIDATES_H_
#define APPROXIMATE_SEARCH_VERIFY_CANDIDATES_H_

#include "approximate_search.h"
#include "gpu_buffer_align_bpm.h"

/*
 * Verify Candidates
 */
void approximate_search_verify_candidates(approximate_search_t* const search,matches_t* const matches);

/*
 * Verify Candidates Buffered
 */
void approximate_search_verify_candidates_buffered_copy(
    approximate_search_t* const search,gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void approximate_search_verify_candidates_buffered_retrieve(
    approximate_search_t* const search,gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_VERIFY_CANDIDATES_H_ */
