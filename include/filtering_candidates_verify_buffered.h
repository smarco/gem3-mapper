/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_verify_buffered.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "archive_text.h"
#include "text_collection.h"
#include "pattern.h"
#include "gpu_buffer_align_bpm.h"

/*
 * BPM-Buffered Add (Candidates Verification)
 */
uint64_t filtering_candidates_verify_buffered_add(
    filtering_candidates_t* const filtering_candidates,pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm);

/*
 * BPM-Buffered Retrieve (Candidates Verification)
 */
uint64_t filtering_candidates_verify_buffered_retrieve(
    filtering_candidates_t* const filtering_candidates,archive_text_t* const archive_text,
    text_collection_t* const text_collection,pattern_t* const pattern,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,const uint64_t candidate_offset_begin,
    const uint64_t candidate_offset_end,mm_stack_t* const mm_stack);

/*
 * Display/Benchmark
 */
void filtering_candidates_verify_buffered_print_benchmark(
    FILE* const stream,filtering_candidates_t* const filtering_candidates,
    archive_text_t* const archive_text,pattern_t* const pattern,
    mm_stack_t* const mm_stack);

