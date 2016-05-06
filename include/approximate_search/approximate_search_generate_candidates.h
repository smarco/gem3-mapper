/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_generate_candidates.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_GENERATE_CANDIDATES_H_
#define APPROXIMATE_SEARCH_GENERATE_CANDIDATES_H_

#include "approximate_search/approximate_search.h"
#include "gpu/gpu_buffer_fmi_decode.h"

/*
 * Generate Candidates from Region-Profile
 */
void approximate_search_generate_candidates_exact(
    approximate_search_t* const search,
    matches_t* const matches);
void approximate_search_generate_candidates_inexact(
    approximate_search_t* const search,
    const bool dynamic_scheduling,
    const bool verify_ahead,
    matches_t* const matches);

/*
 * Buffered Copy/Retrieve
 */
void approximate_search_generate_candidates_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void approximate_search_generate_candidates_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

#endif /* APPROXIMATE_SEARCH_GENERATE_CANDIDATES_H_ */
