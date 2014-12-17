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
 * Select Matches (Retrieving & Processing matches)
 *   - 1. Process CIGAR. (Eg Transform-reverse matches)
 *   - 2. Expand interval-matches (compacted)
 *   - 3. Sort matches wrt distance
 */
GEM_INLINE void archive_select_matches(archive_search_t* const archive_search,matches_t* const matches);

/*
 * Check Matches
 */
GEM_INLINE void archive_check_matches_correctness(
    archive_t* const archive,sequence_t* const sequence,
    matches_t* const matches,mm_stack_t* const mm_stack);
GEM_INLINE void archive_check_matches_optimum(
    archive_t* const archive,sequence_t* const sequence,
    matches_t* const matches,mm_stack_t* const mm_stack);
GEM_INLINE void archive_check_matches_completness(
    archive_t* const archive,sequence_t* const sequence,
    matches_t* const matches,mm_stack_t* const mm_stack);

/*
 * Error Messages
 */
//#define GEM_ERROR_ARCHIVE_SELECT_

#endif /* ARCHIVE_SELECT_H_ */
