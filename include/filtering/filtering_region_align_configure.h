/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align_configure.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */
#ifndef FILTERING_REGION_ALIGN_CONFIGURE_H_
#define FILTERING_REGION_ALIGN_CONFIGURE_H_

#include "filtering/filtering_region.h"
#include "data_structures/pattern.h"
#include "archive/archive_search_parameters.h"

/*
 * Configure Basic Alignment
 */
void filtering_region_align_configure_exact(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    const bool emulated_rc_search);
void filtering_region_align_configure_hamming(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool emulated_rc_search);
void filtering_region_align_configure_levenshtein(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool emulated_rc_search,
    const bool left_gap_alignment,
    mm_stack_t* const mm_stack);

/*
 * Configure SWG-based Alignment
 */
void filtering_region_align_configure_swg(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool emulated_rc_search,
    const bool left_gap_alignment,
    const bool local_alignment,
    mm_stack_t* const mm_stack);

#endif /* FILTERING_REGION_ALIGN_CONFIGURE_H_ */
