/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_bpm_buffer.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates.h"
#include "archive_text.h"
#include "text_collection.h"
#include "pattern.h"
#include "align_bpm_gpu.h"

/*
 * BPM-Buffer API (Candidates Verification)
 */
uint64_t filtering_candidates_bpm_buffer_add(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,bpm_gpu_buffer_t* const bpm_gpu_buffer);
uint64_t filtering_candidates_bpm_buffer_retrieve(
    filtering_candidates_t* const filtering_candidates,archive_text_t* const archive_text,
    text_collection_t* const text_collection,pattern_t* const pattern,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t candidate_offset_begin,
    const uint64_t candidate_offset_end,mm_stack_t* const mm_stack);
