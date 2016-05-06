/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_align_local.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering/filtering_candidates.h"
#include "data_structures/pattern.h"

/*
 * Align Candidates (Local)
 */
void filtering_candidates_align_local(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool emulated_rc_search,
    matches_t* const matches);
