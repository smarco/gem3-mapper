/*
 * PROJECT: GEMMapper
 * FILE: match_align_swg.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGN_SWG_H_
#define MATCH_ALIGN_SWG_H_

#include "utils/essentials.h"
#include "matches/match_align_dto.h"
#include "matches/match_alignment.h"
#include "matches/match_scaffold.h"
#include "matches/matches.h"
#include "matches/matches_cigar.h"

/*
 * Compute alignment type (local/global) wrt identity/score thresholds
 */
void match_align_swg_compute_alignment_type(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_parameters_t* const align_parameters);

/*
 * SWG Alignment
 */
void match_align_swg(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,
    mm_stack_t* const mm_stack);

#endif /* MATCH_ALIGN_SWG_H_ */
