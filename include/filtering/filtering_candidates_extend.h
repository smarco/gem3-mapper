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
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const candidate_pattern,
    const match_trace_t* const extended_match,
    paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end,
    mapper_stats_t* const mapper_stats);
void filtering_candidates_extend_generate_candidates(
    filtering_candidates_t* const extended_filtering_candidates,
    filtering_candidates_t* const candidate_filtering_candidates,
    const pattern_t* const extended_pattern,
    const pattern_t* const candidate_pattern,
    paired_matches_t* const paired_matches,
    mapper_stats_t* const mapper_stats);
