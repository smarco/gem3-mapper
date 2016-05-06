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
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Verify Extend Candidate
 */
uint64_t approximate_search_verify_extend_candidate(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const candidate_pattern,
    const match_trace_t* const extended_match,
    mapper_stats_t* const mapper_stats,
    paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end);

/*
 * Verify Candidates Buffered
 */
void approximate_search_verify_candidates_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void approximate_search_verify_candidates_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_VERIFY_CANDIDATES_H_ */
