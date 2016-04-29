/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_verify_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_VERIFY_CANDIDATES_H_
#define APPROXIMATE_SEARCH_VERIFY_CANDIDATES_H_

#include "approximate_search/approximate_search.h"
#include "gpu/gpu_buffer_align_bpm.h"
#include "matches/paired_matches.h"
#include "mapper/mapper_stats.h"

/*
 * Verify Candidates
 */
void approximate_search_verify_candidates(
    approximate_search_t* const restrict search,
    matches_t* const restrict matches);

/*
 * Verify Extend Candidate
 */
uint64_t approximate_search_verify_extend_candidate(
    filtering_candidates_t* const restrict filtering_candidates,
    pattern_t* const restrict candidate_pattern,
    const match_trace_t* const restrict extended_match,
    mapper_stats_t* const restrict mapper_stats,
    paired_matches_t* const restrict paired_matches,
    const sequence_end_t candidate_end);

/*
 * Verify Candidates Buffered
 */
void approximate_search_verify_candidates_buffered_copy(
    approximate_search_t* const restrict search,
    gpu_buffer_align_bpm_t* const restrict gpu_buffer_align_bpm);
void approximate_search_verify_candidates_buffered_retrieve(
    approximate_search_t* const restrict search,
    gpu_buffer_align_bpm_t* const restrict gpu_buffer_align_bpm,
    matches_t* const restrict matches);

#endif /* APPROXIMATE_SEARCH_VERIFY_CANDIDATES_H_ */
