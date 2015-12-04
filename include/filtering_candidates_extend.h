/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_extend.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "archive_search_parameters.h"
#include "archive.h"
#include "text_collection.h"
#include "pattern.h"
#include "paired_matches.h"

/*
 * Pair Extension
 */
uint64_t filtering_candidates_extend_match(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,const match_trace_t* const extended_match,
    pattern_t* const candidate_pattern,const as_parameters_t* const candidate_actual_parameters,
    mapper_stats_t* const mapper_stats,paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end,mm_stack_t* const mm_stack);
void filtering_candidates_extend_generate_candidates(
    filtering_candidates_t* const extended_filtering_candidates,
    filtering_candidates_t* const candidate_filtering_candidates,
    archive_t* const archive,text_collection_t* const text_collection,
    const pattern_t* const extended_pattern,const pattern_t* const candidate_pattern,
    const search_parameters_t* const search_parameters,mapper_stats_t* const mapper_stats,
    paired_matches_t* const paired_matches,mm_stack_t* const mm_stack);
