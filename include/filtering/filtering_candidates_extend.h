/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_extend.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates.h"
#include "data_structures/pattern.h"
#include "matches/paired_matches.h"

/*
 * Pair Extension
 */
uint64_t filtering_candidates_extend_match(
    filtering_candidates_t* const restrict filtering_candidates,
    pattern_t* const restrict candidate_pattern,
    const match_trace_t* const restrict extended_match,
    paired_matches_t* const restrict paired_matches,
    const sequence_end_t candidate_end,
    mapper_stats_t* const restrict mapper_stats);
void filtering_candidates_extend_generate_candidates(
    filtering_candidates_t* const restrict extended_filtering_candidates,
    filtering_candidates_t* const restrict candidate_filtering_candidates,
    const pattern_t* const restrict extended_pattern,
    const pattern_t* const restrict candidate_pattern,
    paired_matches_t* const restrict paired_matches,
    mapper_stats_t* const restrict mapper_stats);
