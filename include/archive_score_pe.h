/*
 * PROJECT: GEMMapper
 * FILE: archive_score_pe.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SCORE_PE_H_
#define ARCHIVE_SCORE_PE_H_

#include "archive_search.h"
#include "matches.h"

/*
 * PE Score Categories
 */
uint8_t archive_score_matches_pe_default_ties(matches_predictors_t* const predictors);
uint8_t archive_score_matches_pe_default_mmap(matches_predictors_t* const predictors);
uint8_t archive_score_matches_pe_default_unique(matches_predictors_t* const predictors);

/*
 * Archive Scoring PE
 */
void archive_score_matches_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

#endif /* ARCHIVE_SCORE_PE_H_ */
