/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_verify.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "archive_search_parameters.h"
#include "archive.h"
#include "text_collection.h"
#include "pattern.h"

/*
 * Verify Candidates
 */
uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,archive_t* const archive,
    text_collection_t* const text_collection,pattern_t* const pattern,
    const as_parameters_t* const as_parameters,matches_t* const matches,
    mm_stack_t* const mm_stack);
