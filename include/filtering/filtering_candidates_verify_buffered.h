/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_verify_buffered.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates.h"
#include "data_structures/pattern.h"
#include "gpu/gpu_buffer_align_bpm.h"

/*
 * BPM-Buffered Add (Candidates Verification)
 */
void filtering_candidates_verify_buffered_add(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    uint64_t* const gpu_buffer_align_offset,
    filtering_region_buffered_t** const filtering_region_buffered,
    uint64_t* const gpu_num_filtering_regions);

/*
 * BPM-Buffered Retrieve (Candidates Verification)
 */
void filtering_candidates_verify_buffered_retrieve(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset_begin,
    filtering_region_buffered_t* const filtering_region_buffered,
    uint64_t const num_filtering_regions);

/*
 * Display/Benchmark
 */
void filtering_candidates_verify_buffered_print_benchmark(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern);

