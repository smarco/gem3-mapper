/*
 * PROJECT: GEMMapper
 * FILE: archive_select.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SELECT_H_
#define ARCHIVE_SELECT_H_

#include "archive_search.h"
#include "archive_select_parameters.h"

/*
 * Decoding Matches (Retrieving & Processing matches)
 */
void archive_select_decode_trace_matches(
    archive_search_t* const archive_search,matches_t* const matches,
    match_trace_t* const match_trace,const uint64_t num_match_traces,
    const bool force_strand,const strand_t strand);
void archive_select_decode_trace_matches_all(
    archive_search_t* const archive_search,matches_t* const matches,
    const bool force_strand,const strand_t strand);

/*
 * Select Paired-Matches
 */
void archive_select_se_matches(
    archive_search_t* const archive_search,
    const bool paired_mapping,matches_t* const matches);
void archive_select_pe_matches(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

#endif /* ARCHIVE_SELECT_H_ */
