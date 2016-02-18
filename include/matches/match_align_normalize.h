/*
 * PROJECT: GEMMapper
 * FILE: match_align_normalize.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGN_NORMALIZE_H_
#define MATCH_ALIGN_NORMALIZE_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/match_align_dto.h"

/*
 * SWG Normalize CIGAR & Adjust Position (Translate RL-CIGAR if required)
 */
void match_align_normalize(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    mm_stack_t* const mm_stack);

#endif /* MATCH_ALIGN_NORMALIZE_H_ */
