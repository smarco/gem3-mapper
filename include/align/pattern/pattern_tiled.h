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
 */

#ifndef PATTERN_TILED_H_
#define PATTERN_TILED_H_

#include "utils/essentials.h"
#include "align/align_bpm_pattern.h"
#include "align/align_kmer_filter_nway.h"

/*
 * Pattern tiled
 */
typedef struct {
  /* Tile attributes */
  uint64_t tile_offset;
  uint64_t tile_length;
  uint64_t max_error;
  /* BPM Filter */
  bpm_pattern_t bpm_pattern_tile;
} pattern_tile_t;
typedef struct {
  /* Tile dimensions */
  uint64_t num_tiles;
  uint64_t tile_length;
  /* Global filters */
  bpm_pattern_t bpm_pattern;
  kmer_counting_nway_t kmer_filter_nway;
  /* Tiles */
  pattern_tile_t* tiles;
  /* MM */
  mm_allocator_t* mm_allocator;
} pattern_tiled_t;

/*
 * Compile alignment tiles (and filters)
 */
void pattern_tiled_compile(
    pattern_tiled_t* const pattern_tiled,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t max_error,
    const uint64_t kmer_tiles,
    const uint64_t kmer_length,
    const uint64_t kmer_enabled,
    mm_allocator_t* const mm_allocator);
void pattern_tiled_destroy(
    pattern_tiled_t* const pattern_tiled,
    mm_allocator_t* const mm_allocator);

#endif /* PATTERN_TILED_H_ */
