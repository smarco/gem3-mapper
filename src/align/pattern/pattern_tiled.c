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

#include "align/pattern/pattern_tiled.h"
#include "gpu/gpu_config.h"

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
    mm_allocator_t* const mm_allocator) {
  // Init
  pattern_tiled->mm_allocator = mm_allocator;
  // Compute tiling dimensions (adjust to min-128/max-512 tile lengths)
  const uint64_t min_tiles = DIV_CEIL(max_error,BPM_MIN_TILE_LENGTH);
  const uint64_t min_tiles_length = min_tiles * BPM_MIN_TILE_LENGTH;
  uint64_t tile_length = min_tiles_length;
  if (key_length >= tile_length) {
    // Read is larger than the tile (increase tile size)
    const uint64_t prefered_tile_length = BPM_PREFERED_TILE_LENGTH;
    if (key_length < prefered_tile_length) {
      tile_length = prefered_tile_length;
    } else {
      tile_length = BPM_MAX_TILE_LENGTH;
    }
  }
  if (tile_length > BPM_MAX_TILE_LENGTH) tile_length = BPM_MAX_TILE_LENGTH;
  pattern_tiled->tile_length = tile_length;
  pattern_tiled->num_tiles = DIV_CEIL(key_length,tile_length);
  // Compile BPM pattern (Global filters)
  bpm_pattern_compile(&pattern_tiled->bpm_pattern,key,key_length,mm_allocator);
#ifdef GPU_CHECK_KMER_FILTER
  kmer_counting_compile_nway(
      &pattern_tiled->kmer_filter_nway,key,key_length,
      kmer_tiles,kmer_length,true,mm_allocator);
#else
  kmer_counting_compile_nway(
      &pattern_tiled->kmer_filter_nway,key,key_length,
      kmer_tiles,kmer_length,!kmer_enabled,mm_allocator);
#endif
  // Configure tiles
  if (pattern_tiled->num_tiles==1) {
    pattern_tiled->tiles = mm_allocator_alloc(mm_allocator,pattern_tile_t);
    pattern_tiled->tiles->tile_offset = 0;
    pattern_tiled->tiles->tile_length = key_length;
    pattern_tiled->tiles->max_error = MIN(max_error,key_length);
    pattern_tiled->tiles->bpm_pattern_tile = pattern_tiled->bpm_pattern;
  } else {
    // Allocate
    pattern_tiled->tiles = mm_allocator_calloc(mm_allocator,pattern_tiled->num_tiles,pattern_tile_t,false);
    // Compile tiles
    const double max_error_rate = (double)max_error/(double)key_length;
    uint64_t key_length_left = key_length, key_offset = 0;
    uint64_t offset_words64 = 0;
    uint64_t i;
    for (i=0;i<pattern_tiled->num_tiles;++i) {
      // Fetch tile
      pattern_tile_t* const tile = pattern_tiled->tiles + i;
      // Compute dimensions
      const uint64_t actual_tile_length = MIN(tile_length,key_length_left);
      tile->tile_offset = key_offset;
      tile->tile_length = actual_tile_length;
      tile->max_error = (uint64_t)ceil(max_error_rate*(double)actual_tile_length);
      key_length_left -= actual_tile_length;
      key_offset += actual_tile_length;
      // Compile BPM-Tile
      bpm_pattern_compile_tiles(
          &pattern_tiled->bpm_pattern,offset_words64,
          tile->tile_length,&tile->bpm_pattern_tile);
      offset_words64 += tile->bpm_pattern_tile.pattern_num_words64;
    }
  }
}
void pattern_tiled_destroy(
    pattern_tiled_t* const pattern_tiled,
    mm_allocator_t* const mm_allocator) {
  // BPM
  bpm_pattern_destroy(&pattern_tiled->bpm_pattern,mm_allocator);
  // kmer-filter
  kmer_counting_destroy(&pattern_tiled->kmer_filter_nway,mm_allocator);
  mm_allocator_free(mm_allocator,pattern_tiled->tiles);
}
