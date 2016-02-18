/*
 * PROJECT: GEMMapper
 * FILE: match_align_swg_local.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGN_SWG_LOCAL_H_
#define MATCH_ALIGN_SWG_LOCAL_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/match_align_dto.h"

/*
 * SWG-Local Alignment
 */
void match_align_swg_local_alignment(
    matches_t* const matches,
    match_trace_t* const match_trace,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters);

#endif /* MATCH_ALIGN_SWG_LOCAL_H_ */
