/*
 * PROJECT: GEMMapper
 * FILE: archive_search_se.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_SE_H_
#define ARCHIVE_SEARCH_SE_H_

#include "essentials.h"
#include "archive_search.h"

/*
 * Setup
 */
void archive_search_single_end_configure(archive_search_t* const archive_search,mm_search_t* const mm_search);

/*
 * SE Archive Search building blocks
 */
void archive_search_generate_candidates(archive_search_t* const archive_search);
void archive_search_verify_candidates(archive_search_t* const archive_search,matches_t* const matches);
void archive_search_finish_search(archive_search_t* const archive_search,matches_t* const matches);
void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer);
void archive_search_retrieve_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer,matches_t* const matches);

/*
 * Single-End Indexed Search (SE Online Approximate String Search)
 */
void archive_search_single_end(archive_search_t* const archive_search,matches_t* const matches);

/*
 * Compute Predictors
 */
void archive_search_compute_predictors(
    archive_search_t* const archive_search,matches_t* const matches,
    matches_predictors_t* const predictors);

/*
 * Errors
 */
#define GEM_ERROR_ARCHIVE_SEARCH_INDEX_COMPLEMENT_REQUIRED "Archive Search. Explicit indexed complement required"

#endif /* ARCHIVE_SEARCH_SE_H_ */
