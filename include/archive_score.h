/*
 * PROJECT: GEMMapper
 * FILE: archive_score.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SCORE_H_
#define ARCHIVE_SCORE_H_

#include "archive_search.h"
#include "matches.h"

/*
 * GEM Score
 */
GEM_INLINE uint8_t archive_score_matches_gem_se_ties(matches_predictors_t* const predictors);
GEM_INLINE uint8_t archive_score_matches_gem_se_mmap(matches_predictors_t* const predictors);
GEM_INLINE uint8_t archive_score_matches_gem_se_unique(matches_predictors_t* const predictors);

/*
 * SE Scoring
 */
GEM_INLINE void archive_score_matches_se(
    archive_search_t* const archive_search,
    const bool paired_mapping,matches_t* const matches);

/*
 * PE Scoring
 */
GEM_INLINE void archive_score_matches_pe(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

/*
 * Error Messages
 */
//#define GEM_ERROR_ARCHIVE_SCORE_

#endif /* ARCHIVE_SCORE_H_ */
