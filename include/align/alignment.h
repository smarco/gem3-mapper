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
#include "align/align_bpm_pattern.h"
#include "text/text_collection.h"

/*
 * Constants
 */
#define ALIGN_DISTANCE_INF     UINT32_MAX
#define ALIGN_COLUMN_INF       UINT64_MAX
#define ALIGN_DISABLED         (UINT32_MAX-1)

/*
 * Region Alignment
 */
typedef struct {
  uint64_t distance;                        // Distance
  uint64_t text_begin_offset;               // Text begin offset
  uint64_t text_end_offset;                 // Text end offset
} alignment_tile_t;
typedef struct {
  uint64_t num_tiles;                // Total number of tiles
  uint64_t distance_min_bound;       // Distance min-bound (Sum all tile distances)
  alignment_tile_t* alignment_tiles; // Alignment of all tiles
} alignment_t;

/*
 * Check matches (CIGAR string against text & pattern)
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
    const bool verbose);

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
 * Verify levenshtein using BPM
 */
void alignment_verify_levenshtein_bpm(
    alignment_t* const alignment,
    const uint64_t filtering_max_error,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    text_trace_t* const text_trace);

#endif /* ALIGNMENT_H_ */
