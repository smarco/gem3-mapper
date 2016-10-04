/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_compose.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_SCAFFOLD_COMPOSE_H_
#define MATCH_SCAFFOLD_COMPOSE_H_

#include "utils/essentials.h"
#include "matches/scaffold/match_scaffold.h"
#include "matches/align/match_alignment.h"

/*
 * Add alignment-region
 */
void match_scaffold_compose_add_approximate_match(
    match_scaffold_t* const match_scaffold,
    const uint64_t key_begin,
    const uint64_t key_end,
    const uint64_t text_begin,
    const uint64_t text_end,
    const uint64_t cigar_offset,
    const uint64_t cigar_length,
    const uint64_t score);
match_alignment_region_t* match_scaffold_compose_add_exact_match(
    match_scaffold_t* const match_scaffold,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t cigar_offset,
    const uint64_t match_length);
void match_scaffold_compose_add_mismatch(
    match_scaffold_t* const match_scaffold,
    match_alignment_region_t* const last_alignment_region,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t match_length);

#endif /* MATCH_SCAFFOLD_COMPOSE_H_ */
