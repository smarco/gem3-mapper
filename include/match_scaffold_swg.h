/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_swg.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_SCAFFOLD_SWG_H_
#define MATCH_SCAFFOLD_SWG_H_

#include "essentials.h"
#include "match_scaffold.h"
#include "match_align_dto.h"
#include "matches.h"

/*
 * SWG Scaffolding (based on SWG-distance)
 */
bool match_scaffold_smith_waterman_gotoh(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack);

#endif /* MATCH_SCAFFOLD_SWG_H_ */
