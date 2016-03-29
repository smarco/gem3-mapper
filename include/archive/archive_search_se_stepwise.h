/*
 * PROJECT: GEMMapper
 * FILE: archive_search_se_stepwise.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_SE_STEPWISE_H_
#define ARCHIVE_SEARCH_SE_STEPWISE_H_

#include "utils/essentials.h"
#include "archive/archive_search.h"
#include "gpu/gpu_buffer_fmi_bsearch.h"
#include "gpu/gpu_buffer_fmi_decode.h"
#include "gpu/gpu_buffer_align_bpm.h"
#include "matches/matches.h"

/*
 * Stepwise: Init Search
 */
void archive_search_se_stepwise_init_search(archive_search_t* const archive_search);

/*
 * Stepwise: Region-Profile
 */
void archive_search_se_stepwise_region_profile_generate_static(
    archive_search_t* const archive_search);
void archive_search_se_stepwise_region_profile_generate_adaptive(
    archive_search_t* const archive_search);
void archive_search_se_stepwise_region_profile_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);
void archive_search_se_stepwise_region_profile_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);

/*
 * Stepwise: Decode-Candidates
 */
void archive_search_se_stepwise_decode_candidates_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void archive_search_se_stepwise_decode_candidates_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

/*
 * Stepwise: Verify-Candidates
 */
void archive_search_se_stepwise_verify_candidates_copy(
    archive_search_t* const archive_search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void archive_search_se_stepwise_verify_candidates_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches);

/*
 * Stepwise: Finish Search
 */
void archive_search_se_stepwise_finish_search(
    archive_search_t* const archive_search,
    matches_t* const matches);

#endif /* ARCHIVE_SEARCH_SE_STEPWISE_H_ */
