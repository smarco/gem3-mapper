/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_process_buffered.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_CANDIDATES_PROCESS_BUFFERED_H_
#define FILTERING_CANDIDATES_PROCESS_BUFFERED_H_

#include "filtering_candidates.h"
#include "archive_search_parameters.h"
#include "archive.h"
#include "pattern.h"
#include "gpu_buffer_fmi_decode.h"

/*
 * Decode Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_decode_filtering_positions_buffered(
    filtering_candidates_t* const filtering_candidates,locator_t* const locator,
    archive_text_t* const archive_text,fm_index_t* const fm_index,
    region_search_t* const region_search,pattern_t* const pattern,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin,mm_stack_t* const mm_stack);
void filtering_candidates_decode_filtering_positions_buffered_prefetched(
    filtering_candidates_t* const filtering_candidates,locator_t* const locator,
    archive_text_t* const archive_text,fm_index_t* const fm_index,
    region_search_t* const region_search,pattern_t* const pattern,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin,mm_stack_t* const mm_stack);

/*
 * Process Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_process_candidates_buffered(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    const bool compose_region_chaining,mm_stack_t* const mm_stack);

#endif /* FILTERING_CANDIDATES_PROCESS_BUFFERED_H_ */
