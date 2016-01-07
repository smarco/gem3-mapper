/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_compose.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_SCAFFOLD_COMPOSE_H_
#define MATCH_SCAFFOLD_COMPOSE_H_

#include "essentials.h"
#include "match_scaffold.h"
#include "match_alignment.h"

bool match_scaffold_compose_is_homopolymer(
    const uint8_t* const text,const uint64_t text_pos,
    const uint64_t homopolymer_min_context,const uint8_t* indel_text,
    const uint64_t indel_length);
void match_scaffold_compose_add_matching_approximate(
    const uint64_t key_begin,const uint64_t key_end,
    const uint64_t text_begin,const uint64_t text_end,
    const uint64_t cigar_offset,const uint64_t cigar_length,
    const uint64_t score,match_scaffold_t* const match_scaffold);
region_matching_t* match_scaffold_compose_add_matching_exact(
    uint8_t* const text,uint64_t* const key_pos,uint64_t* const text_pos,
    const uint64_t cigar_offset,const uint64_t match_length,
    match_scaffold_t* const match_scaffold);
void match_scaffold_compose_add_mismatch(
    uint64_t* const key_pos,uint64_t* const text_pos,const uint64_t match_length,
    match_scaffold_t* const match_scaffold,region_matching_t* const last_scaffold_region);
void match_scaffold_compose_add_homopolymer(
    uint8_t* const text,uint64_t* const key_pos,uint64_t* const text_pos,
    const uint64_t context_cigar_offset,const uint64_t context_match_length,
    const cigar_t indel_type,const uint64_t indel_length,
    match_scaffold_t* const match_scaffold,region_matching_t* const last_scaffold_region);

#endif /* MATCH_SCAFFOLD_COMPOSE_H_ */
