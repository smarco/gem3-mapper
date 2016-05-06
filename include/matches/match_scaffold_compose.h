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
#include "matches/match_scaffold.h"
#include "matches/match_alignment.h"

/*
 * Add Scaffold Region
 */
void match_scaffold_compose_add_matching_approximate(
    match_scaffold_t* const match_scaffold,
    const uint64_t key_begin,
    const uint64_t key_end,
    const uint64_t text_begin,
    const uint64_t text_end,
    const uint64_t cigar_offset,
    const uint64_t cigar_length,
    const uint64_t score);
region_matching_t* match_scaffold_compose_add_matching_exact(
    match_scaffold_t* const match_scaffold,
    uint8_t* const text,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t cigar_offset,
    const uint64_t match_length);
void match_scaffold_compose_add_mismatch(
    match_scaffold_t* const match_scaffold,
    region_matching_t* const last_scaffold_region,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t match_length);

#endif /* MATCH_SCAFFOLD_COMPOSE_H_ */
