/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_align.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates.h"
#include "data_structures/pattern.h"

/*
 * Filtering Candidates Align Region
 */
bool filtering_candidates_align_region(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const region,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    const bool local_alignment,
    const bool extended_match,
    matches_t* const matches);

/*
 * Filtering Candidates Align
 */
uint64_t filtering_candidates_align_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    const bool extended_match,
    const bool local_alignment,
    matches_t* const matches);
