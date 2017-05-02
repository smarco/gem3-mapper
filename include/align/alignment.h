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
 *   Basic alignment module to check edit-distance based alignments/CIGARs
 *   and to generate the best correct CIGAR by means of a simple algorithm
 *   as to double-check
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "utils/essentials.h"
#include "align/pattern/pattern_tiled.h"

/*
 * Constants
 */
#define ALIGN_COLUMN_INF            UINT64_MAX
#define ALIGN_DISTANCE_INF         (UINT32_MAX)
#define ALIGN_DISTANCE_UNKNOWN     (UINT32_MAX-1)

/*
 * Alignment (tiled)
 */
typedef struct {
  uint64_t distance;                        // Distance
  uint64_t text_begin_offset;               // Text begin offset
  uint64_t text_end_offset;                 // Text end offset
} alignment_tile_t;
typedef struct {
  uint64_t num_tiles;                       // Total number of tiles
  uint64_t distance_min_bound;              // Distance min-bound (Sum all tile distances)
  uint64_t distance_rank;
  alignment_tile_t* alignment_tiles;        // Alignment of all tiles
} alignment_t;

/*
 * Setup
 */
void alignment_init(
    alignment_t* const alignment,
    const uint64_t key_length,
    const uint64_t text_begin_offset,
    const uint64_t text_end_offset,
    const uint64_t max_error,
    const uint64_t num_tiles,
    const uint64_t tile_length);

/*
 * Compute edit distance (Basic DP-Matrix Alignment)
 */
int64_t alignment_compute_edit_distance(
    const char* const key,
    const uint64_t key_length,
    const char* const text,
    const uint64_t text_length,
    const bool ends_free,
    uint64_t* const position);

/*
 * Verify edit distance using kmer-counting
 */
void alignment_verify_edit_kmer(
    alignment_t* const alignment,
    pattern_tiled_t* const pattern_tiled,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error,
    const uint64_t kmer_tiles,
    const uint64_t kmer_length);

/*
 * Verify edit distance using BPM
 */
void alignment_verify_edit_bpm(
    alignment_t* const alignment,
    pattern_tiled_t* const pattern_tiled,
    uint8_t* const key,
    uint8_t* const text,
    const uint64_t max_error);

/*
 * Check alignment (CIGAR string against text & pattern)
 */
bool alignment_check(
    FILE* const stream,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    uint64_t const cigar_offset,
    uint64_t const cigar_length,
    const bool local_alignment,
    const bool verbose);

/*
 * Display
 */
void alignment_print_pretty(
    FILE* const stream,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    uint64_t const cigar_offset,
    uint64_t const cigar_length,
    mm_allocator_t* const mm_allocator);

#endif /* ALIGNMENT_H_ */
