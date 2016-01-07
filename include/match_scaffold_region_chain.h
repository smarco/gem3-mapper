/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_region_chain.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_SCAFFOLD_REGION_CHAIN_H_
#define MATCH_SCAFFOLD_REGION_CHAIN_H_

#include "essentials.h"
#include "match_align_dto.h"
#include "match_scaffold.h"
#include "matches.h"

/*
 * Region-Chain Scaffolding
 */
bool match_scaffold_region_chain(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack);

#endif /* MATCH_SCAFFOLD_REGION_CHAIN_H_ */
