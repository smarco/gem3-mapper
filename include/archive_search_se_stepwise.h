/*
 * PROJECT: GEMMapper
 * FILE: archive_search_se_stepwise.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_SE_STEPWISE_H_
#define ARCHIVE_SEARCH_SE_STEPWISE_H_

#include "essentials.h"
#include "archive_search.h"
#include "align_bpm_gpu.h"
#include "matches.h"

/*
 * SE Archive Search Stepwise: Init Search
 */
void archive_search_se_stepwise_init_search(archive_search_t* const archive_search);

/*
 * SE Archive Search Stepwise: Region-Profile Generation
 */
void archive_search_se_stepwise_generate_region_profile_partition(
    archive_search_t* const archive_search);
void archive_search_se_stepwise_generate_region_profile_copy_partition(
    archive_search_t* const archive_search);
void archive_search_se_stepwise_generate_region_profile_retrieve_partition(
    archive_search_t* const archive_search);

/*
 * SE Archive Search Stepwise: Candidate Verification
 */
void archive_search_se_stepwise_generate_candidates(archive_search_t* const archive_search);
void archive_search_se_stepwise_verify_candidates(
    archive_search_t* const archive_search,matches_t* const matches);
void archive_search_se_stepwise_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer);
void archive_search_se_stepwise_retrieve_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    matches_t* const matches);

/*
 * SE Archive Search Stepwise: Finish Search
 */
void archive_search_se_stepwise_finish_search(archive_search_t* const archive_search,matches_t* const matches);


#endif /* ARCHIVE_SEARCH_SE_STEPWISE_H_ */
