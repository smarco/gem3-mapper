/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_REGION_ALIGN_H_
#define FILTERING_REGION_ALIGN_H_

#include "filtering_region.h"
#include "archive_search_parameters.h"
#include "archive_text.h"
#include "pattern.h"
#include "matches.h"

/*
 * Region (Re)Align by clone previous
 */
void filtering_region_align_clone(
    match_trace_t* const match_trace_src,match_trace_t* const match_trace_dst,
    filtering_region_t* const filtering_region_dst,text_collection_t* const text_collection,
    const uint64_t run_length);

/*
 * Adjust distance bound by scaffolding
 */
void filtering_region_align_adjust_distance_by_scaffolding(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const as_parameters_t* const as_parameters,pattern_t* const pattern,
    matches_t* const matches,text_collection_t* const text_collection,
    mm_stack_t* const mm_stack);

/*
 * Region (Re)Align
 */
bool filtering_region_align(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const text_collection_t* const text_collection,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,pattern_t* const pattern,matches_t* const matches,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack);
bool filtering_region_align_unbounded(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    const text_collection_t* const text_collection,const as_parameters_t* const as_parameters,
    const bool emulated_rc_search,pattern_t* const pattern,matches_t* const matches,
    match_trace_t* const match_trace,mm_stack_t* const mm_stack);

#endif /* FILTERING_REGION_ALIGN_H_ */
