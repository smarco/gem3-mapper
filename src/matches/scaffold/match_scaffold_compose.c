/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_compose.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/scaffold/match_scaffold_compose.h"

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
    const uint64_t score) {
  // Add alignment-region
  match_alignment_region_t* const match_alignment_region =
      match_scaffold->alignment_regions + match_scaffold->num_alignment_regions;
  ++match_scaffold->num_alignment_regions;
  // Setup matching-region
  match_alignment_region->region_type = match_alignment_region_approximate;
  match_alignment_region->error = score;
  match_alignment_region->key_begin = key_begin;
  match_alignment_region->key_end = key_end;
  match_alignment_region->text_begin = text_begin;
  match_alignment_region->text_end = text_end;
  match_alignment_region->cigar_buffer_offset = cigar_offset;
  match_alignment_region->cigar_length = cigar_length;
}
match_alignment_region_t* match_scaffold_compose_add_exact_match(
    match_scaffold_t* const match_scaffold,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t cigar_offset,
    const uint64_t match_length) {
  // Add alignment-region
  match_alignment_region_t* const match_alignment_region =
      match_scaffold->alignment_regions + match_scaffold->num_alignment_regions;
  ++match_scaffold->num_alignment_regions;
  // Set-up alignment-region
  match_alignment_region->region_type = match_alignment_region_exact;
  match_alignment_region->key_begin = *key_offset;
  match_alignment_region->text_begin = *text_offset;
  match_alignment_region->cigar_buffer_offset = cigar_offset;
  match_alignment_region->cigar_length = 1;
  *key_offset += match_length;
  *text_offset += match_length;
  match_scaffold->scaffolding_coverage += match_length;
  match_alignment_region->key_end = *key_offset;
  match_alignment_region->text_end = *text_offset;
  // Return alignment-region
  return match_alignment_region;
}
void match_scaffold_compose_add_mismatch(
    match_scaffold_t* const match_scaffold,
    match_alignment_region_t* const last_alignment_region,
    uint64_t* const key_offset,
    uint64_t* const text_offset,
    const uint64_t match_length) {
  // Extend matching-region
  last_alignment_region->region_type = match_alignment_region_approximate;
  last_alignment_region->cigar_length += 2; // Add the mismatch + matching
  *key_offset += match_length;
  *text_offset += match_length;
  match_scaffold->scaffolding_coverage += match_length;
  last_alignment_region->key_end = *key_offset;
  last_alignment_region->text_end = *text_offset;
}
