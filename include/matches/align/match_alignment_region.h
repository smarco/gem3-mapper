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
  match_alignment_region_exact       = 0,
  match_alignment_region_approximate = 1,
} match_alignment_region_type;

typedef struct {
  /* Alignment-region type */
  uint16_t _error;
  uint16_t _type;
  /* CIGAR alignment-region */
  uint32_t _cigar_buffer_offset;
  uint32_t _cigar_length;
  /* Coordinates of the region */
  uint32_t _key_begin;
  uint32_t _key_end;
  uint32_t _text_begin;
  uint32_t _text_end;
} match_alignment_region_t;

/*
 * Accessors
 */
void match_alignment_region_init(
    match_alignment_region_t* const match_alignment_region,
    const match_alignment_region_type type,
    const uint64_t error,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length,
    const uint64_t key_begin,
    const uint64_t key_end,
    const uint64_t text_begin,
    const uint64_t text_end);

match_alignment_region_type match_alignment_region_get_type(
    const match_alignment_region_t* const match_alignment_region);
void match_alignment_region_set_type(
    match_alignment_region_t* const match_alignment_region,
    const match_alignment_region_type type);
uint64_t match_alignment_region_get_error(
    const match_alignment_region_t* const match_alignment_region);
void match_alignment_region_set_error(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t error);

uint64_t match_alignment_region_get_cigar_buffer_offset(
    const match_alignment_region_t* const match_alignment_region);
uint64_t match_alignment_region_get_cigar_length(
    const match_alignment_region_t* const match_alignment_region);
uint64_t match_alignment_region_get_key_begin(
    const match_alignment_region_t* const match_alignment_region);
uint64_t match_alignment_region_get_key_end(
    const match_alignment_region_t* const match_alignment_region);
uint64_t match_alignment_region_get_text_begin(
    const match_alignment_region_t* const match_alignment_region);
uint64_t match_alignment_region_get_text_end(
    const match_alignment_region_t* const match_alignment_region);

void match_alignment_region_set_cigar_buffer_offset(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t cigar_buffer_offset);
void match_alignment_region_set_cigar_length(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t cigar_length);
void match_alignment_region_set_key_begin(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t key_begin);
void match_alignment_region_set_key_end(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t key_end);
void match_alignment_region_set_text_begin(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t text_begin);
void match_alignment_region_set_text_end(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t text_end);

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
