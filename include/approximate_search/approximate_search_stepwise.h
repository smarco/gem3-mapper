/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering_adaptive_stepwise.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_STEPWISE_H_
#define APPROXIMATE_SEARCH_STEPWISE_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"
#include "gpu/gpu_buffer_fmi_bsearch.h"
#include "gpu/gpu_buffer_fmi_decode.h"
#include "gpu/gpu_buffer_align_bpm.h"

/*
 * AM Stepwise :: Region Profile
 */
void approximate_search_stepwise_region_profile_generate_static(approximate_search_t* const search);
void approximate_search_stepwise_region_profile_generate_adaptive(approximate_search_t* const search);
void approximate_search_stepwise_region_profile_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);
void approximate_search_stepwise_region_profile_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search);

/*
 * AM Stepwise :: Decode Candidates
 */
void approximate_search_stepwise_decode_candidates_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void approximate_search_stepwise_decode_candidates_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

/*
 * AM Stepwise :: Verify Candidates
 */
void approximate_search_stepwise_verify_candidates_copy(
    approximate_search_t* const search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);
void approximate_search_stepwise_verify_candidates_retrieve(
    approximate_search_t* const search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches);

/*
 * AM Stepwise :: Finish Search
 */
void approximate_search_stepwise_finish(
    approximate_search_t* const search,
    matches_t* const matches);

#endif /* APPROXIMATE_SEARCH_STEPWISE_H_ */
