/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_SCAFFOLD_H_
#define MATCH_SCAFFOLD_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/match_align_dto.h"

/*
 * Match Scaffold Type
 */
typedef enum {
  scaffold_none,
  scaffold_region_chain,
  scaffold_levenshtein,
  scaffold_swg,
} match_scaffold_type;

/*
 * Match Scaffold
 */
typedef struct {
  /* Scaffold Properties */
  match_scaffold_type scaffold_type;
  bool scaffold_regions_rl;
  /* Scaffold Matching Regions */
  region_matching_t* scaffold_regions;
  uint64_t num_scaffold_regions;
  uint64_t scaffolding_coverage;
  /* Underlying Alignment */
  match_alignment_t match_alignment;
} match_scaffold_t;

/*
 * Setup
 */
void match_scaffold_init(match_scaffold_t* const match_scaffold);

/*
 * Accessors
 */
bool match_scaffold_is_null(match_scaffold_t* const match_scaffold);

/*
 * Adaptive Scaffolding of the alignment (Best effort)
 */
void match_scaffold_adaptive(
    match_scaffold_t* const match_scaffold,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    matches_t* const matches,
    mm_stack_t* const mm_stack);

/*
 * Sorting
 */
void match_scaffold_sort_regions_matching(match_scaffold_t* const match_scaffold);

/*
 * Display
 */
void match_scaffold_print(
    FILE* const stream,
    matches_t* const matches,
    match_scaffold_t* const match_scaffold);
void match_scaffold_print_pretty(
    FILE* const stream,
    matches_t* const matches,
    match_scaffold_t* const match_scaffold,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    mm_stack_t* const mm_stack);

#endif /* MATCH_SCAFFOLD_H_ */
