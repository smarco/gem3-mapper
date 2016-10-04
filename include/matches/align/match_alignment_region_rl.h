/*
 * PROJECT: GEMMapper
 * FILE: match_alignment_region_rl.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGNMENT_REGION_RL_H_
#define MATCH_ALIGNMENT_REGION_RL_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/align/match_align_dto.h"
#include "matches/align/match_alignment_region.h"

/*
 * RL-Space Translate CIGAR
 */
void match_alignment_region_rl_translate(
    match_alignment_region_t* const match_alignment_region,
    match_align_input_t* const align_input,
    const bool left_gap_alignment,
    vector_t* const cigar_vector);

#endif /* MATCH_ALIGNMENT_REGION_RL_H_ */
