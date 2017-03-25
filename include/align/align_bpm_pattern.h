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
 *   Alignment module providing BPM-pattern data structure used by BPM
 *   algorithms (full or tiled)
 */

#ifndef ALIGN_BPM_PATTERN_H_
#define ALIGN_BPM_PATTERN_H_

#include "utils/essentials.h"

/*
 * BPM Pattern
 */
typedef struct _bpm_pattern_t bpm_pattern_t;
struct _bpm_pattern_t {
  /* BMP Pattern */
  uint64_t* PEQ;                // Pattern equalities (Bit vector for Myers-DP)
  uint64_t pattern_length;      // Length
  uint64_t pattern_num_words64; // ceil(Length / |w|)
  uint64_t pattern_mod;         // Length % |w|
  /* BPM Auxiliary data */
  uint64_t* P;
  uint64_t* M;
  uint64_t* level_mask;
  int64_t* score;
  int64_t* init_score;
  uint64_t* pattern_left;
  /* BPM tiles (Pattern split in tiles) */
  uint64_t num_pattern_tiles;  // Total number of tiles
  uint64_t tile_length;        // Ideal tile length (the last can be shorter)
};

/*
 * Constants
 */
#define BPM_MIN_TILE_LENGTH 128
#define BPM_PREFERED_TILE_LENGTH 256
#define BPM_MAX_TILE_LENGTH 512

/*
 * Pattern Accessors
 */
#define BPM_PATTERN_PEQ_IDX(word_pos,encoded_character)   ((word_pos*DNA__N_RANGE)+(encoded_character))
#define BPM_PATTERN_BDP_IDX(position,num_words,word_pos)  ((position)*(num_words)+(word_pos))

/*
 * Compile Pattern
 */
void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,
    uint8_t* const pattern,
    const uint64_t pattern_length,
    mm_allocator_t* const mm_allocator);
void bpm_pattern_destroy(
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator);

/*
 * Compile Pattern Tiles
 */
void bpm_pattern_compile_tiles(
    bpm_pattern_t* const bpm_pattern,
    const uint64_t offset_words64,
    const uint64_t tile_length,
    bpm_pattern_t* const bpm_pattern_tile);

#endif /* ALIGN_BPM_PATTERN_H_ */
