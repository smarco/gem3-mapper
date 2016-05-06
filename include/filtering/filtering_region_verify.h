/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_verify.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_REGION_VERIFY_H_
#define FILTERING_REGION_VERIFY_H_

#include "filtering/filtering_candidates.h"
#include "filtering/filtering_region.h"
#include "data_structures/pattern.h"

/*
 * Verify Levenshtein
 */
void filtering_region_verify_levenshtein_bpm(
    filtering_candidates_t* const filtering_candidates,
    region_alignment_t* const region_alignment,
    const uint64_t filtering_max_error,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    text_trace_t* const text_trace,
    mm_stack_t* const mm_stack);

/*
 * Region Verification
 */
bool filtering_region_verify(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const bool kmer_filter);
uint64_t filtering_region_verify_extension(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t text_trace_offset,
    const uint64_t index_position,
    const pattern_t* const pattern);

#endif /* FILTERING_REGION_VERIFY_H_ */
