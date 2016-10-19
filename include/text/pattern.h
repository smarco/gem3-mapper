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

#ifndef PATTERN_H_
#define PATTERN_H_

#include "utils/essentials.h"
#include "align/align_bpm_pattern.h"
#include "archive/search/archive_search_se_parameters.h"
#include "align/alignment_filters.h"

/*
 * Approximate Search Pattern
 */
typedef struct {
  /* Processed Search Pattern */
  uint8_t* key;                      // Encoded Pattern
  uint8_t* quality_mask;             // Quality Mask
  uint64_t key_length;               // Total Length
  uint64_t clip_left;                // Input sequence left-clipped
  uint64_t clip_right;               // Input sequence right-clipped
  /* Run-Length Pattern */
  bool run_length;
  uint8_t* rl_key;                   // RL-Encoded Text
  uint64_t rl_key_length;            // RL-Encoded Text length
  uint32_t* rl_runs_acc;             // Length of each run (accumulated)
  /* Pattern Properties */
  uint64_t num_wildcards;
  uint64_t num_low_quality_bases;
  uint64_t num_non_canonical_bases;
  uint64_t max_effective_filtering_error;
  uint64_t max_effective_bandwidth;
  /* Alignment filters */
  alignment_filters_t alignment_filters;
} pattern_t;

/*
 * Tiled Pattern
 */
typedef struct {
  // Current tile dimensions
  uint64_t tile_offset;
  uint64_t tile_wide;
  uint64_t tile_tall;
  uint64_t tile_next_offset_inc;
  // Complete tile dimensions (defaults)
  uint64_t pattern_band_width;
  uint64_t pattern_tile_offset;
  uint64_t pattern_tile_wide;
  uint64_t pattern_tile_tall;
  uint64_t pattern_max_error;
  uint64_t pattern_remaining_length;
  uint64_t sequence_length;
  // Tile matching information
  uint64_t tile_distance;
  uint64_t tile_match_column;
  uint64_t prev_tile_match_position;
} pattern_tiled_t;

/*
 * Pattern Prepare
 */
void pattern_init(
    pattern_t* const pattern,
    sequence_t* const sequence,
    bool* const do_quality_search,
    const search_parameters_t* const parameters,
    const bool run_length_pattern,
    const bool kmer_filter_compile,
    mm_stack_t* const mm_stack);
void pattern_clear(pattern_t* const pattern);
bool pattern_is_null(pattern_t* const pattern);

/*
 * Pattern Tiling
 */
void pattern_tiled_init(
    pattern_tiled_t* const pattern_tiled,
    const uint64_t pattern_length,
    const uint64_t pattern_tile_length,
    const uint64_t sequence_length,
    const uint64_t max_error);
void pattern_tiled_calculate_next(pattern_tiled_t* const pattern_tiled);
uint64_t pattern_tiled_bound_matching_path(pattern_tiled_t* const pattern_tiled);

/*
 * Display
 */
void pattern_enc_print(
    FILE* const stream,
    const uint8_t* const key,
    const uint64_t key_length);

#endif /* PATTERN_H_ */
