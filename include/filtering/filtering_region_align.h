/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_align.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_REGION_ALIGN_H_
#define FILTERING_REGION_ALIGN_H_

#include "filtering/filtering_candidates.h"
#include "filtering/filtering_region.h"
#include "data_structures/pattern.h"
#include "matches/matches.h"

/*
 * Region (Re)Align by clone previous
 */
void filtering_region_align_clone(
    match_trace_t* const match_trace_src,
    match_trace_t* const match_trace_dst,
    filtering_region_t* const filtering_region_dst,
    const uint64_t run_length);

/*
 * Region (Re)Align
 */
bool filtering_region_align(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    const bool local_alignment,
    matches_t* const matches,
    match_trace_t* const match_trace);

#endif /* FILTERING_REGION_ALIGN_H_ */
