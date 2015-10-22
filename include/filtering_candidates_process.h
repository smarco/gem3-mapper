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
#include "archive.h"
#include "pattern.h"
#include "search_parameters.h"

/*
 * Retrieve all candidates(text) from the index
 */
void filtering_candidates_retrieve_filtering_regions(
    filtering_candidates_t* const filtering_candidates,archive_text_t* const archive_text,
    text_collection_t* const text_collection,mm_stack_t* const mm_stack);

/*
 * Process Candidates
 */
uint64_t filtering_candidates_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    const bool compose_region_chaining,mm_stack_t* const mm_stack);

#endif /* FILTERING_CANDIDATES_PROCESS_H_ */
