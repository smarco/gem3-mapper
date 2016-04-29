/*
 * PROJECT: GEMMapper
 * FILE: archive_select.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SELECT_H_
#define ARCHIVE_SELECT_H_

#include "archive/archive_search.h"
#include "archive/archive_select_parameters.h"

/*
 * Setup
 */
void archive_select_configure_se(archive_search_t* const restrict archive_search);
void archive_select_configure_pe(archive_search_t* const restrict archive_search);

/*
 * Select Paired-Matches
 */
void archive_select_se_matches(
    archive_search_t* const restrict archive_search,
    select_parameters_t* const restrict select_parameters,
    matches_t* const restrict matches);
void archive_select_pe_matches(
    archive_search_t* const restrict archive_search_end1,
    archive_search_t* const restrict archive_search_end2,
    select_parameters_t* const restrict select_parameters,
    paired_matches_t* const restrict paired_matches);

#endif /* ARCHIVE_SELECT_H_ */
