/*
 * PROJECT: GEMMapper
 * FILE: match_align_rl.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGN_RL_H_
#define MATCH_ALIGN_RL_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/match_align_dto.h"

/*
 * RL-Space CIGAR Translation
 */
void match_align_rl_translate_cigar(
    match_trace_t* const match_trace,
    vector_t* const cigar_vector,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters);

#endif /* MATCH_ALIGN_RL_H_ */
