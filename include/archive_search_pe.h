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
 * Memory Injection (Support Data Structures)
 */
void archive_search_pe_inject_mm(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    mm_search_t* const mm_search);

/*
 * PE Extension Control
 */
bool archive_search_pe_is_extension_feasible(archive_search_t* const archive_search);
bool archive_search_pe_use_shortcut_extension(archive_search_t* const archive_search,matches_t* const matches);
bool archive_search_pe_use_recovery_extension(archive_search_t* const archive_search,matches_t* const matches);

/*
 * Archive Search PE Continue Search
 */
void archive_search_pe_continue(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

/*
 * Paired-End Indexed Search (PE Online Approximate String Search)
 */
void archive_search_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

/*
 * Compute Predictors
 */
void archive_search_pe_compute_predictors(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,matches_predictors_t* const predictors);

/*
 * Display
 */
void archive_search_pe_print(
    FILE* const stream,archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,paired_matches_t* const paired_matches);

#endif /* ARCHIVE_SEARCH_PE_H_ */
