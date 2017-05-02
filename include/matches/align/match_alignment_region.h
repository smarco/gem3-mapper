/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
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
  uint16_t error;
  uint16_t type;
  /* CIGAR alignment-region */
  uint32_t cigar_buffer_offset;
  uint32_t cigar_length;
  /* Coordinates of the region */
  uint32_t key_begin;
  uint32_t key_end;
  uint32_t text_begin;
  uint32_t text_end;
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
int match_alignment_region_text_coverage(
    match_alignment_region_t* const match_alignment_region);
int match_alignment_region_text_distance(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b);
bool match_alignment_region_text_overlap(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b);

/*
 * Compare
 */
int64_t match_alignment_region_cmp_key_position(
    const match_alignment_region_t* const a,
    const match_alignment_region_t* const b);
int64_t match_alignment_region_cmp_text_position(
    const match_alignment_region_t* const a,
    const match_alignment_region_t* const b);
/*
 * Check
 */
void match_alignment_region_check(
    match_alignment_region_t* const match_alignment_region,
    uint8_t* const key,
    uint8_t* const text,
    vector_t* const cigar_vector);
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
