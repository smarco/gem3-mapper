/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_align.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "archive_text.h"
#include "text_collection.h"
#include "pattern.h"

/*
 * Alignment Accepted Candidate
 */
uint64_t filtering_candidates_align_accepted_regions(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,const bool extended_match,
    matches_t* const matches,mm_stack_t* const mm_stack);

/*
 * Align Candidates
 */
uint64_t filtering_candidates_align_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    const bool approximated_distance,const bool extended_match,
    matches_t* const matches,mm_stack_t* const mm_stack);
uint64_t filtering_candidates_align_unbounded(
    filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,const locator_t* const locator,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const bool emulated_rc_search,const as_parameters_t* const as_parameters,
    matches_t* const matches,mm_stack_t* const mm_stack);
