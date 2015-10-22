/*
 * PROJECT: GEMMapper
 * FILE: archive_search_pe_stepwise.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_PE_STEPWISE_H_
#define ARCHIVE_SEARCH_PE_STEPWISE_H_

#include "essentials.h"
#include "archive_search.h"
#include "paired_matches.h"

/*
 * Archive Search PE Stepwise
 */
void archive_search_pe_stepwise_generate_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_search_pe_stepwise_finish_search(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

#endif /* ARCHIVE_SEARCH_PE_STEPWISE_H_ */
