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
 * RL-Space Translate CIGAR
 */
void match_align_rl_translate_region_cigar(
    region_matching_t* const region_matching,
    match_align_input_t* const align_input,
    const bool left_gap_alignment,
    vector_t* const cigar_vector);

#endif /* MATCH_ALIGN_RL_H_ */
