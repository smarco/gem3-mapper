/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_levenshtein.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_SCAFFOLD_LEVENSHTEIN_H_
#define MATCH_SCAFFOLD_LEVENSHTEIN_H_

#include "utils/essentials.h"
#include "matches/match_scaffold.h"
#include "matches/match_align_dto.h"
#include "matches/matches.h"

/*
 * Levenshtein Scaffold Tiled
 */
bool match_scaffold_levenshtein(
    match_scaffold_t* const restrict match_scaffold,
    match_align_input_t* const restrict align_input,
    match_align_parameters_t* const restrict align_parameters,
    matches_t* const restrict matches,
    mm_stack_t* const restrict mm_stack);

#endif /* MATCH_SCAFFOLD_LEVENSHTEIN_H_ */
