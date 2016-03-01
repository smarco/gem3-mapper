/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_compose.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/match_scaffold_compose.h"

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
    const uint64_t score) {
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
    match_scaffold_t* const match_scaffold,
    uint8_t* const text,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t cigar_offset,
    const uint64_t match_length) {
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
    match_scaffold_t* const match_scaffold,
    region_matching_t* const last_scaffold_region,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t match_length) {
  // Extend matching-region
  last_scaffold_region->matching_type = region_matching_approximate;
  last_scaffold_region->cigar_length += 2; // Add the mismatch + matching
  *key_offset += match_length;
  *text_offset += match_length;
  match_scaffold->scaffolding_coverage += match_length;
  last_scaffold_region->key_end = *key_offset;
  last_scaffold_region->text_end = *text_offset;
}
