/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_compose.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/match_scaffold_compose.h"

bool match_scaffold_compose_is_homopolymer(
    const uint8_t* const text,
    const uint64_t text_pos,
    const uint64_t homopolymer_min_context,
    const uint8_t* indel_text,
    const uint64_t indel_length) {
  if (text_pos < homopolymer_min_context) return false;
  // Search for minimum context support for an homopolymer event
  const uint64_t context_init_position = text_pos - homopolymer_min_context;
  const uint64_t context_end_position = text_pos;
  const uint8_t monomer = text[context_init_position];
  // Check context
  uint64_t j;
  for (j=context_init_position+1;j<=context_end_position;j++) {
    if (monomer!=text[j]) return false;
  }
  // Check indel-text
  for (j=0;j<indel_length;++j) {
    if (monomer!=indel_text[j]) return false;
  }
  return true;
}
void match_scaffold_compose_add_matching_approximate(
    const uint64_t key_begin,
    const uint64_t key_end,
    const uint64_t text_begin,
    const uint64_t text_end,
    const uint64_t cigar_offset,
    const uint64_t cigar_length,
    const uint64_t score,
    match_scaffold_t* const match_scaffold) {
  // Add matching region
  region_matching_t* const region_matching = match_scaffold->scaffold_regions + match_scaffold->num_scaffold_regions;
  ++match_scaffold->num_scaffold_regions;
  // Setup matching-region
  region_matching->matching_type = region_matching_approximate;
  region_matching->error = score;
  region_matching->key_begin = key_begin;
  region_matching->key_end = key_end;
  region_matching->text_begin = text_begin;
  region_matching->text_end = text_end;
  region_matching->cigar_buffer_offset = cigar_offset;
  region_matching->cigar_length = cigar_length;
}
region_matching_t* match_scaffold_compose_add_matching_exact(
    uint8_t* const text,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t cigar_offset,
    const uint64_t match_length,
    match_scaffold_t* const match_scaffold) {
  // Add matching region
  region_matching_t* const region_matching = match_scaffold->scaffold_regions + match_scaffold->num_scaffold_regions;
  ++match_scaffold->num_scaffold_regions;
  // Set-up matching region
  region_matching->matching_type = region_matching_exact;
  region_matching->key_begin = *key_offset;
  region_matching->text_begin = *text_offset;
  region_matching->cigar_buffer_offset = cigar_offset;
  region_matching->cigar_length = 1;
  *key_offset += match_length;
  *text_offset += match_length;
  match_scaffold->scaffolding_coverage += match_length;
  region_matching->key_end = *key_offset;
  region_matching->text_end = *text_offset;
  // Return matching region
  return region_matching;
}
void match_scaffold_compose_add_mismatch(
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t match_length,
    match_scaffold_t* const match_scaffold,
    region_matching_t* const last_scaffold_region) {
  // Extend matching-region
  last_scaffold_region->matching_type = region_matching_approximate;
  last_scaffold_region->cigar_length += 2; // Add the mismatch + matching
  *key_offset += match_length;
  *text_offset += match_length;
  match_scaffold->scaffolding_coverage += match_length;
  last_scaffold_region->key_end = *key_offset;
  last_scaffold_region->text_end = *text_offset;
}
void match_scaffold_compose_add_homopolymer(
    uint8_t* const text,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t context_cigar_offset,
    const uint64_t context_match_length,
    const cigar_t indel_type,
    const uint64_t indel_length,
    match_scaffold_t* const match_scaffold,
    region_matching_t* const last_scaffold_region) {
  // Check previous matching-region
  if (last_scaffold_region!=NULL) {
    region_matching_t* const region_matching = last_scaffold_region;
    region_matching->matching_type = region_matching_approximate;
    ++region_matching->cigar_length;
    if (indel_type==cigar_ins) {
      *text_offset += indel_length;
      region_matching->text_end = *text_offset;
    } else { // cigar_del
      *key_offset += indel_length;
      region_matching->key_end = *key_offset;
      match_scaffold->scaffolding_coverage += indel_length;
    }
  } else {
    // Add matching region
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + match_scaffold->num_scaffold_regions;
    ++match_scaffold->num_scaffold_regions;
    // Set-up homopolymer matching region
    region_matching->matching_type = region_matching_approximate;
    region_matching->cigar_buffer_offset = context_cigar_offset;
    region_matching->cigar_length = 2;
    region_matching->key_begin = *key_offset - context_match_length;
    region_matching->text_begin = *text_offset - context_match_length;
    if (indel_type==cigar_ins) {
      *text_offset += indel_length;
      match_scaffold->scaffolding_coverage += context_match_length;
    } else { // cigar_del
      *key_offset += indel_length;
      match_scaffold->scaffolding_coverage += context_match_length + indel_length;
    }
    region_matching->key_end = *key_offset;
    region_matching->text_end = *text_offset;
  }
}
