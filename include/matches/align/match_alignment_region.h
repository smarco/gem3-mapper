/*
 * PROJECT: GEMMapper
 * FILE: match_alignment_region.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCH_ALIGNMENT_REGION_H_
#define MATCH_ALIGNMENT_REGION_H_

#include "utils/essentials.h"
#include "matches/matches.h"

/*
 * Match Alignment Region
 *   A chain of regions aligning determines the shape of the final
 *   alignment CIGAR. If allowed, the chain of alignment-regions will
 *   be followed as to rule out regions (key<->text) that match on the high level
 */
typedef enum {
  match_alignment_region_exact,
  match_alignment_region_approximate,
} match_alignment_region_type;

typedef struct {
  /* Alignment-region type */
  match_alignment_region_type region_type;
  /* Error of the alignment-region */
  uint64_t error;
  uint64_t cigar_buffer_offset;
  uint64_t cigar_length;
  /* Coordinates of the region */
  uint64_t key_begin;
  uint64_t key_end;
  uint64_t text_begin;
  uint64_t text_end;
} match_alignment_region_t;

/*
 * Key/Text Operators
 */
uint64_t match_alignment_region_text_coverage(
    match_alignment_region_t* const match_alignment_region);
uint64_t match_alignment_region_text_distance(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b);
bool match_alignment_region_text_overlap(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b);

/*
 * Compare
 */
int match_alignment_region_key_cmp(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b);
int match_alignment_region_cmp_text_position(
    const match_alignment_region_t* const a,
    const match_alignment_region_t* const b);

/*
 * Display
 */
void match_alignment_region_print(
    FILE* const stream,
    match_alignment_region_t* const match_alignment_region,
    const uint64_t match_alignment_region_id,
    matches_t* const matches);
void match_alignment_region_print_pretty(
    FILE* const stream,
    match_alignment_region_t* const match_alignment_region,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text);

#endif /* MATCH_ALIGNMENT_REGION_H_ */
