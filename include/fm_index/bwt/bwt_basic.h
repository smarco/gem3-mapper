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
 * DESCRIPTION:
 *   BWT data structure (basic implementation). Provides basic
 *   routines to encode a DNA text into a Burrows-Wheeler
 *   transform using a compact bit representation and counter
 *   buckets as to enhance Occ/rank queries
 *   MODEL:
 *     NAME:      bwt_basic
 *     ALPHABET:  3bits (8 chars)
 *     LEVELS:    2-levels
 *     STEP:      1s
 *     COUNTERS:  (8x64b,8x16b)
 *     BITMAP:    64bits
 *     SAMPLED:   No
 *     OTHERS:    -
 */

#ifndef BWT_BASIC_H_
#define BWT_BASIC_H_

#include "fm_index/bwt/bwt_commons.h"
#include "text/dna_text.h"

/*
 * BWT Structure
 */
typedef struct {
  /* Meta-Data */
  uint64_t length;                    // Length of the BWT
  uint64_t* c;                        // Occurrences of each character
  uint64_t* C;                        // The cumulative occurrences ("ranks") of the symbols of the string
  /* BWT Buckets-Structure */
  uint64_t num_minor_blocks;          // Total number of minor blocks
  uint64_t num_mayor_blocks;          // Total number of mayor blocks
  uint64_t* mayor_counters;           // Pointer to the Mayor Counters (Rank)
  uint64_t* bwt_mem;                  // Pointer to the BWT structure in memory
  /* MM */
  mm_t* mm_counters;
  mm_t* mm_bwt_mem;
} bwt_basic_t;
typedef struct {
  /* Meta-Data */
  bwt_basic_t bwt;
  /* Auxiliary variables */
  uint64_t mayor_counter;             // Current mayorCounter position
  uint64_t* minor_block_mem;          // Pointer to current minorBlock position
  uint64_t* mayor_occ;                // Mayor accumulated counts
  uint64_t* minor_occ;                // Minor accumulated counts
} bwt_basic_builder_t;

/*
 * BWT Builder
 */
bwt_basic_builder_t* bwt_basic_builder_new(
    dna_text_t* const bwt_text,
    const uint64_t* const character_occurrences,
    const bool verbose);
void bwt_basic_builder_write(
    fm_t* const file_manager,
    bwt_basic_builder_t* const bwt_builder);
void bwt_basic_builder_delete(bwt_basic_builder_t* const bwt_builder);

/*
 * BWT Loader
 */
bwt_basic_t* bwt_basic_read_mem(mm_t* const memory_manager,const bool check);
void bwt_basic_delete(bwt_basic_t* const bwt);

/*
 * BWT General Accessors
 */
uint64_t bwt_basic_builder_get_length(const bwt_basic_builder_t* const bwt_builder);
uint64_t bwt_basic_builder_get_size(bwt_basic_builder_t* const bwt_builder);

uint64_t bwt_basic_get_length(const bwt_basic_t* const bwt);
uint64_t bwt_basic_get_size(bwt_basic_t* const bwt);

bool bwt_basic_is_same_bucket(
    const uint64_t lo,
    const uint64_t hi);

/*
 * BWT Character Accessors
 */
uint8_t bwt_basic_char(
    const bwt_basic_t* const bwt,
    const uint64_t position);
char bwt_basic_char_character(
    const bwt_basic_t* const bwt,
    const uint64_t position);

/*
 * BWT ERank (Exclusive Rank Function)
 */
uint64_t bwt_basic_builder_erank(
    const bwt_basic_builder_t* const bwt_builder,
    const uint8_t char_enc,
    const uint64_t position);
uint64_t bwt_basic_erank(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t position);
uint64_t bwt_basic_erank_character(
    const bwt_basic_t* const bwt,
    const char character,
    const uint64_t position);
void bwt_basic_erank_interval(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out);

/*
 * BWT Prefetched ERank
 */
void bwt_basic_prefetch(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    bwt_block_locator_t* const block_loc);
uint64_t bwt_basic_prefetched_erank(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t position,
    const bwt_block_locator_t* const block_loc);
void bwt_basic_prefetched_erank_interval(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc);

/*
 *  BWT Precomputed ERank (Precomputation of the block's elements)
 */
void bwt_basic_precompute(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms);
void bwt_basic_precompute_interval(
    const bwt_basic_t* const bwt,
    const uint64_t lo,
    const uint64_t hi,
    bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms);
void bwt_basic_prefetched_precompute(
    const bwt_basic_t* const bwt,
    const bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms);
void bwt_basic_prefetched_precompute_interval(
    const bwt_basic_t* const bwt,
    const uint64_t lo,
    const bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms);

uint64_t bwt_basic_precomputed_erank(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const bwt_block_locator_t* const block_loc,
    const bwt_block_elms_t* const block_elms);
void bwt_basic_precomputed_erank_interval(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    uint64_t* const lo_out,
    uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc,
    const bwt_block_elms_t* const block_elms);

/*
 * BWT LF (Last to first)
 */
uint64_t bwt_basic_LF(
    const bwt_basic_t* const bwt,
    const uint64_t position);
uint64_t bwt_basic_prefetched_LF(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    const bwt_block_locator_t* const block_loc);

uint64_t bwt_basic_LF__enc(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    uint8_t* const char_enc);
uint64_t bwt_basic_LF__character(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    char* const character);
uint64_t bwt_basic_prefetched_LF__enc(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    uint8_t* const char_enc,
    const bwt_block_locator_t* const block_loc);

/*
 * Display
 */
void bwt_basic_builder_print(
    FILE* const stream,
    bwt_basic_builder_t* const bwt_builder);
void bwt_basic_print(
    FILE* const stream,
    bwt_basic_t* const bwt);

#endif /* BWT_BASIC_H_ */
