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
#include "mapper_stats.h"
#include "paired_matches.h"

/*
 * Verify Candidates
 */
void approximate_search_verify_candidates(approximate_search_t* const search,matches_t* const matches);

/*
 * Verify Extend Candidate
 */
uint64_t approximate_search_verify_extend_candidate(
    filtering_candidates_t* const filtering_candidates,archive_text_t* const archive_text,
    const locator_t* const locator,text_collection_t* const text_collection,
    const match_trace_t* const extended_match,pattern_t* const candidate_pattern,
    const as_parameters_t* const candidate_actual_parameters,mapper_stats_t* const mapper_stats,
    paired_matches_t* const paired_matches,const sequence_end_t candidate_end,
    mm_stack_t* const mm_stack);

/*
 * Verify Candidates Buffered
 */
void approximate_search_verify_candidates_buffered_copy(
    approximate_search_t* const search,gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void approximate_search_verify_candidates_buffered_retrieve(
    approximate_search_t* const search,gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_VERIFY_CANDIDATES_H_ */
