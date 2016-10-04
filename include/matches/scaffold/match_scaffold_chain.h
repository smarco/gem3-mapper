/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_chain.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_SCAFFOLD_CHAIN_H_
#define MATCH_SCAFFOLD_CHAIN_H_

#include "utils/essentials.h"
#include "matches/align/match_align_dto.h"
#include "matches/scaffold/match_scaffold.h"
#include "matches/matches.h"

/*
 * Alignment-region chain scaffolding
 */
void match_scaffold_chain(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    const bool exact_extend,
    mm_stack_t* const mm_stack);

#endif /* MATCH_SCAFFOLD_CHAIN_H_ */
