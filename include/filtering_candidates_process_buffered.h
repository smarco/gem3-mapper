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
#include "archive.h"
#include "pattern.h"
#include "search_parameters.h"
#include "gpu_buffer_fmi_decode.h"

/*
 * Process Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_decode_filtering_positions_buffered(
    filtering_candidates_t* const filtering_candidates,const locator_t* const locator,
    archive_text_t* const archive_text,const fm_index_t* const fm_index,
    const uint64_t region_begin,const uint64_t region_end,
    const uint64_t region_lo,const uint64_t region_hi,
    const uint64_t key_length,const uint64_t boundary_error,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin);
void filtering_candidates_process_candidates_buffered(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    const bool compose_region_chaining,mm_stack_t* const mm_stack);

#endif /* FILTERING_CANDIDATES_PROCESS_BUFFERED_H_ */
