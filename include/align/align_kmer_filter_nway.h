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
 *   Filter based on general k-mer counting as to quickly filter out
 *   candidates that cannot align against its region-text
 */

#ifndef ALIGN_KMER_FILTER_NWAY_H_
#define ALIGN_KMER_FILTER_NWAY_H_

#include "utils/essentials.h"

/*
 * Kmer counting filter
 */
typedef struct {
  uint64_t text_begin;                    // Begin position of the chunk
  uint64_t text_end;                      // End position of the chunk
  uint64_t num_key_kmers;                 // Total kmers in chunk
  uint64_t curr_text_kmers;               // Current number of kmers contained in text (wrt pattern profile)
  uint64_t max_text_kmers;                // Maximum number of kmers contained in text (wrt pattern profile)
} kmer_counting_text_tile_t;
typedef struct {
  uint64_t begin;                         // Begin position of the tile
  uint64_t end;                           // End position of the tile
} kmer_counting_key_tile_t;
typedef struct {
  // State
  bool enabled;                           // Enabled kmer-filtering
  // Filter parameters
  uint64_t kmer_length;                   // Kmer length
  uint64_t kmer_mask;                     // Kmer mask to extract kmer offset
  uint64_t num_kmers;                     // Total number of possible kmers in table
  // Key
  uint8_t* key;                           // Key
  uint64_t key_length;                    // Key length
  uint64_t key_tile_length;               // Key chunk length
  // Key tiled
  uint64_t num_tiles;                     // Total number of tiles
  kmer_counting_key_tile_t* key_tiles;    // Key tiles (N-way)
  // Text tiled
  uint64_t sliding_window_length;         // Kmer window length
  kmer_counting_text_tile_t* text_tiles;  // Text tiles (N-way)
  // Profile tables
  uint16_t* kmer_count_text;              // Text profile (kmers on text)
  uint16_t* kmer_count_pattern;           // Key chunks profile (kmers on each key chunk)
} kmer_counting_nway_t;

/*
 * Setup
 */
void kmer_counting_compile_nway(
    kmer_counting_nway_t* const kmer_counting,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t num_tiles,
    const uint64_t kmer_length,
    const bool count_pattern_kmers,
    mm_allocator_t* const mm_allocator);
void kmer_counting_destroy(
    kmer_counting_nway_t* const kmer_counting,
    mm_allocator_t* const mm_allocator);

/*
 * Tiling
 */
void kmer_counting_prepare_tiling(
    kmer_counting_nway_t* const kmer_counting,
    const uint64_t text_length,
    const uint64_t max_error);

/*
 * Min-error bound
 */
uint64_t kmer_counting_min_bound_nway(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error,
    mm_allocator_t* const mm_allocator);

#endif /* ALIGN_KMER_FILTER_NWAY_H_ */
