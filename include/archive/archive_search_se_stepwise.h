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
#include "gpu/gpu_buffer_fmi_asearch.h"
#include "gpu/gpu_buffer_fmi_ssearch.h"
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
void archive_search_se_stepwise_region_profile_static_generate(
    archive_search_t* const archive_search);
void archive_search_se_stepwise_region_profile_static_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void archive_search_se_stepwise_region_profile_static_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);

void archive_search_se_stepwise_region_profile_adaptive_generate(
    archive_search_t* const archive_search);
void archive_search_se_stepwise_region_profile_adaptive_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
void archive_search_se_stepwise_region_profile_adaptive_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

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
    matches_t* const matches,
    const bool paired_end_search);

#endif /* ARCHIVE_SEARCH_SE_STEPWISE_H_ */
