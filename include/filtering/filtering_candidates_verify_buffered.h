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
uint64_t filtering_candidates_verify_buffered_add(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

/*
 * BPM-Buffered Retrieve (Candidates Verification)
 */
uint64_t filtering_candidates_verify_buffered_retrieve(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    const uint64_t candidate_offset_begin,
    const uint64_t candidate_offset_end);

/*
 * Display/Benchmark
 */
void filtering_candidates_verify_buffered_print_benchmark(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern);

