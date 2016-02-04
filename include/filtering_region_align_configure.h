/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align_configure.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */
#ifndef FILTERING_REGION_ALIGN_CONFIGURE_H_
#define FILTERING_REGION_ALIGN_CONFIGURE_H_

#include "filtering_region.h"
#include "pattern.h"
#include "archive_search_parameters.h"

/*
 * Configure Basic Alignment
 */
void filtering_region_align_configure_exact(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,const bool emulated_rc_search);
void filtering_region_align_configure_hamming(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,text_trace_t* const text_trace,const bool emulated_rc_search);
void filtering_region_align_configure_levenshtein(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,text_trace_t* const text_trace,
    const bool emulated_rc_search,const bool left_gap_alignment);

/*
 * Configure SWG-based Alignment
 */
void filtering_region_align_configure_scaffold(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,text_trace_t* const text_trace,
    const bool left_gap_alignment);
void filtering_region_align_configure_swg(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    filtering_region_t* const filtering_region,const as_parameters_t* const as_parameters,
    pattern_t* const pattern,text_trace_t* const text_trace,
    const bool emulated_rc_search,const bool left_gap_alignment);

#endif /* FILTERING_REGION_ALIGN_CONFIGURE_H_ */
