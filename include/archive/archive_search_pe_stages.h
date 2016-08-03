/*
 * PROJECT: GEMMapper
 * FILE: archive_search_pe_stages.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_PE_STAGES_H_
#define ARCHIVE_SEARCH_PE_STAGES_H_

#include "utils/essentials.h"
#include "archive/archive_search.h"
#include "matches/paired_matches.h"

/*
 * Archive Search PE Stages
 */
void archive_search_pe_begin(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_search_pe_search_end1(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_search_pe_search_end2(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_search_pe_recovery(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_search_pe_find_pairs(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_search_pe_end(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

#endif /* ARCHIVE_SEARCH_PE_STAGES_H_ */
