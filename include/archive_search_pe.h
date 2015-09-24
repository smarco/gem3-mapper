/*
 * PROJECT: GEMMapper
 * FILE: archive_search_pe.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_PE_H_
#define ARCHIVE_SEARCH_PE_H_

#include "essentials.h"
#include "archive_search.h"
#include "paired_matches.h"

/*
 * Setup
 */
void archive_search_paired_end_configure(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    mm_search_t* const mm_search);

/*
 * PE Archive Search building blocks
 */
void archive_search_pe_generate_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_search_pe_finish_search(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

/*
 * Paired-End Indexed Search (PE Online Approximate String Search)
 */
void archive_search_paired_end(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

/*
 * Compute Predictors
 */
void archive_search_paired_end_compute_predictors(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,matches_predictors_t* const predictors);

/*
 * Display
 */
GEM_INLINE void archive_search_pe_print(
    FILE* const stream,archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,paired_matches_t* const paired_matches);

/*
 * Errors
 */
#define GEM_ERROR_ARCHIVE_SEARCH_INDEX_COMPLEMENT_REQUIRED "Archive Search. Explicit indexed complement required"

#endif /* ARCHIVE_SEARCH_PE_H_ */
