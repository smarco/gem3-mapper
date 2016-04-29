/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_process_buffered.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_CANDIDATES_PROCESS_BUFFERED_H_
#define FILTERING_CANDIDATES_PROCESS_BUFFERED_H_

#include "filtering/filtering_candidates.h"
#include "data_structures/pattern.h"
#include "gpu/gpu_buffer_fmi_decode.h"

/*
 * Decode Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_decode_sa_filtering_positions_buffered(
    filtering_candidates_t* const restrict filtering_candidates,
    pattern_t* const restrict pattern,
    region_search_t* const restrict region_search,
    filtering_position_buffered_t* const restrict gpu_filtering_positions,
    gpu_buffer_fmi_decode_t* const restrict gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin);
void filtering_candidates_decode_text_filtering_positions_buffered(
    filtering_candidates_t* const restrict filtering_candidates,
    pattern_t* const restrict pattern,
    region_search_t* const restrict region_search,
    filtering_position_buffered_t* const restrict gpu_filtering_positions,
    gpu_buffer_fmi_decode_t* const restrict gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin);

/*
 * Process Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_process_candidates_buffered(
    filtering_candidates_t* const restrict filtering_candidates,
    pattern_t* const restrict pattern,
    const bool compose_region_chaining);

#endif /* FILTERING_CANDIDATES_PROCESS_BUFFERED_H_ */
