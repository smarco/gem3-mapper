/*
 * PROJECT: GEMMapper
 * FILE: archive_check.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_CHECK_H_
#define ARCHIVE_CHECK_H_

#include "archive/archive.h"
#include "archive/archive_search_parameters.h"
#include "data_structures/sequence.h"
#include "matches/match_align.h"
#include "matches/paired_matches.h"

/*
 * Check Matches
 */
void archive_check_se_matches(
    archive_t* const archive,
    const alignment_model_t alignment_model,
    swg_penalties_t* swg_penalties,
    sequence_t* const sequence,
    matches_t* const matches,
    const archive_check_type check_type,
    mm_stack_t* const mm_stack);
void archive_check_pe_matches(
    archive_t* const archive,
    const alignment_model_t alignment_model,
    swg_penalties_t* swg_penalties,
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,
    paired_matches_t* const paired_matches,
    const archive_check_type check_type,
    mm_stack_t* const mm_stack);

#endif /* ARCHIVE_CHECK_H_ */
