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
 *
 *   Verification tests can be heuristic, like discarding a candidate region if
 *   if doesn't have enough seed supporting it. Out it can be exact like k-mer counting.
 *   Both can give an upper/lower bound on the matching distance.
 *     Seed-counting => At least 1 error per seed not matching (it needs at least
 *                      MismsRegions to make a match)
 *     K-mer counting => ...
 */

#include "align/alignment.h"
#include "align/align_kmer_filter_nway.h"
#include "align/pattern/pattern.h"
#include "text/dna_text.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Flags
 */
// #define KMER_COUNTING_SLIDING_WINDOW

/*
 * Constants
 */
#define KMER_COUNTING_MASK_3   0x000000000000003Full /* 256 B */
#define KMER_COUNTING_MASK_4   0x00000000000000FFull /*  1 KB */
#define KMER_COUNTING_MASK_5   0x00000000000003FFull /*  4 KB */
#define KMER_COUNTING_MASK_6   0x0000000000000FFFull /* 16 KB */
#define KMER_COUNTING_MASK_7   0x0000000000003FFFull /* 64 KB */

#define KMER_COUNTING_MASK_INDEX(kmer_idx) ((kmer_idx) & kmer_counting->kmer_mask)
#define KMER_COUNTING_ADD_INDEX(kmer_idx,enc_char) kmer_idx = (kmer_idx<<2 | (enc_char))
#define KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char) kmer_idx = KMER_COUNTING_MASK_INDEX(kmer_idx<<2 | (enc_char))

/*
 * Uncalled bases handling
 */
#define KMER_COUNTING_FILTER_CHAR(enc_char) ((enc_char) % ENC_DNA_CHAR_N)

/*
 * Compile Pattern
 */
void kmer_counting_compile_nway_allocate_tables(
    kmer_counting_nway_t* const kmer_counting,
    const bool count_pattern_kmers,
    mm_allocator_t* const mm_allocator) {
  // Profile tables
  const uint64_t num_tiles = kmer_counting->num_tiles;
  if (count_pattern_kmers) {
    const uint64_t key_tiles_size = num_tiles * sizeof(kmer_counting_key_tile_t);
    const uint64_t text_tiles_size = num_tiles * sizeof(kmer_counting_text_tile_t);
    const uint64_t kmer_table_size = kmer_counting->num_kmers*num_tiles*sizeof(uint16_t);
    void* const memory = mm_allocator_calloc(mm_allocator,key_tiles_size+text_tiles_size+2*kmer_table_size,uint8_t,true);
    kmer_counting->key_tiles = memory;
    kmer_counting->text_tiles = memory + key_tiles_size;
    kmer_counting->kmer_count_text = memory + key_tiles_size + text_tiles_size;
    kmer_counting->kmer_count_pattern = memory + key_tiles_size + text_tiles_size + kmer_table_size;
  } else {
    const uint64_t key_tiles_size = num_tiles * sizeof(kmer_counting_key_tile_t);
    const uint64_t text_tiles_size = num_tiles * sizeof(kmer_counting_text_tile_t);
    void* const memory = mm_allocator_calloc(mm_allocator,key_tiles_size+text_tiles_size,uint8_t,true);
    kmer_counting->key_tiles = memory;
    kmer_counting->text_tiles = memory + key_tiles_size;
    kmer_counting->kmer_count_text = NULL;
    kmer_counting->kmer_count_pattern = NULL;
  }
}
void kmer_counting_compile_nway(
    kmer_counting_nway_t* const kmer_counting,
    uint8_t* const key,
    const uint64_t key_length,
    const uint64_t num_tiles,
    const uint64_t kmer_length,
    const bool count_pattern_kmers,
    mm_allocator_t* const mm_allocator) {
  // Filter parameters
  kmer_counting->kmer_length = kmer_length;
  switch (kmer_length) { // Check kmer length
    case 3: kmer_counting->kmer_mask = KMER_COUNTING_MASK_3; break;
    case 4: kmer_counting->kmer_mask = KMER_COUNTING_MASK_4; break;
    case 5: kmer_counting->kmer_mask = KMER_COUNTING_MASK_5; break;
    case 6: kmer_counting->kmer_mask = KMER_COUNTING_MASK_6; break;
    default:
      gem_warn_msg("K-mer counting. Invalid proposed k-mer length. Set to k=7");
      kmer_counting->kmer_length = 7;
      // no break
    case 7: kmer_counting->kmer_mask = KMER_COUNTING_MASK_7; break;
  }
  kmer_counting->num_kmers = POW4(kmer_counting->kmer_length);
  // Pattern parameters
  kmer_counting->key = key;
  kmer_counting->key_length = key_length;
  kmer_counting->num_tiles = num_tiles;
  kmer_counting->key_tile_length = DIV_CEIL(key_length,num_tiles);
  // Check min-tile length
  if (kmer_counting->key_tile_length < kmer_counting->kmer_length) {
    kmer_counting->enabled = false;
    kmer_counting->key_tiles = NULL;
    return;
  }
  kmer_counting->enabled = true;
  // Profile tables
  kmer_counting_compile_nway_allocate_tables(kmer_counting,count_pattern_kmers,mm_allocator);
  // Count kmers in pattern (chunks)
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  uint64_t chunk_idx = 0, chunk_offset = 0;
  kmer_counting->sliding_window_length = 0;
  while (chunk_idx < num_tiles) {
    // Compute chunk length
    const uint64_t chunk_end = (chunk_idx < num_tiles-1) ?
        chunk_offset + kmer_counting->key_tile_length : key_length;
    kmer_counting->key_tiles[chunk_idx].begin = chunk_offset;
    kmer_counting->key_tiles[chunk_idx].end = chunk_end;
    const uint64_t actual_key_tile_length = chunk_end-chunk_offset;
    kmer_counting->sliding_window_length = MAX(kmer_counting->sliding_window_length,actual_key_tile_length);
    kmer_counting->text_tiles[chunk_idx].num_key_kmers = actual_key_tile_length - (kmer_counting->kmer_length-1);
    // Count until chunk end
    if (count_pattern_kmers) {
      uint64_t pos, kmer_idx=0, acc=0;
      for (pos=chunk_offset;pos<chunk_end;++pos) {
        const uint8_t enc_char = key[pos];
        if (gem_expect_false(enc_char==ENC_DNA_CHAR_N)) {
          acc = 0;
        } else {
          KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char); // Update kmer-index
          if (acc < kmer_counting->kmer_length-1) {
            ++acc; // Inc accumulator
          } else {
            ++(kmer_count_pattern[kmer_idx*num_tiles+chunk_idx]); // Increment kmer-count
          }
        }
      }
    }
    // Next chunk
    chunk_offset = chunk_end;
    ++chunk_idx;
  }
}
void kmer_counting_destroy(
    kmer_counting_nway_t* const kmer_counting,
    mm_allocator_t* const mm_allocator) {
  if (kmer_counting->key_tiles!=NULL) {
    mm_allocator_free(mm_allocator,kmer_counting->key_tiles);
    kmer_counting->key_tiles = NULL;
  }
}
/*
 * K-mer counting tiling
 */
void kmer_counting_prepare_tiling(
    kmer_counting_nway_t* const kmer_counting,
    const uint64_t text_length,
    const uint64_t max_error) {
  const uint64_t num_tiles = kmer_counting->num_tiles;
  if (num_tiles==1) {
    kmer_counting->text_tiles[0].max_text_kmers = 0;
    kmer_counting->text_tiles[0].curr_text_kmers = 0;
  } else {
    const uint64_t adjusted_max_error = MIN(max_error,kmer_counting->key_tile_length);
    pattern_tiling_t pattern_tiling;
    pattern_tiling_init(
        &pattern_tiling,kmer_counting->key_length,
        kmer_counting->key_tile_length,text_length,adjusted_max_error);
    kmer_counting_text_tile_t* const text_tiles = kmer_counting->text_tiles;
    uint64_t chunk_idx;
    for (chunk_idx=0;chunk_idx<num_tiles;++chunk_idx) {
      // Setup chunk
      text_tiles[chunk_idx].text_begin = pattern_tiling.tile_offset;
      text_tiles[chunk_idx].text_end = pattern_tiling.tile_offset + pattern_tiling.tile_wide;
      text_tiles[chunk_idx].max_text_kmers = 0;
      text_tiles[chunk_idx].curr_text_kmers = 0;
      // Calculate next tile
      pattern_tiling_next(&pattern_tiling);
    }
  }
}
/*
 * K-mer counting operators
 */
/*
 * Kmer counting
 */
#define KMER_COUNTING_INC_COUNT(count_text,count_pattern,kmers_in_text,max_kmers_in_text) \
  if (gem_expect_false(count_pattern > 0 && count_text < count_pattern)) { \
    ++kmers_in_text; \
    max_kmers_in_text = MAX(max_kmers_in_text,kmers_in_text);\
  }
#define KMER_COUNTING_DEC_COUNT(count_text,count_pattern,kmers_in_text) \
  if (gem_expect_false(count_pattern > 0 && count_text <= count_pattern)) { \
    --kmers_in_text; \
  }
void kmer_counting_decrement_kmer_count(
    kmer_counting_text_tile_t* const chunk,
    const uint64_t kmer_begin,
    const uint64_t kmer_end,
    const uint64_t sliding_window_length,
    uint16_t* const text_count_ptr,
    uint16_t* const pattern_count_ptr) {
  const uint16_t text_count = *text_count_ptr;
  const uint16_t pattern_count = *pattern_count_ptr;
  // KMER_COUNTING_DEC_COUNT(text_count,pattern_count,chunk->curr_text_kmers);
  if (gem_expect_false(pattern_count > 0 && text_count <= pattern_count)) {
    --(chunk->curr_text_kmers);
  }
  --(*text_count_ptr);
}
#define KMER_COUNTING_DECREMENT(chunk_idx) \
  /* Check Tiling boundaries & Update kmer counters */ \
  if (chunks[chunk_idx].text_begin+sliding_window_length <= kmer_begin && \
      kmer_end < text_tiles[chunk_idx].text_end) { \
    kmer_counting_decrement_kmer_count( \
        text_tiles+chunk_idx,kmer_begin,kmer_end,sliding_window_length, \
        text_count_ptr+chunk_idx,pattern_count_ptr+chunk_idx); \
  }
void kmer_counting_increment_kmer_count(
    kmer_counting_text_tile_t* const chunk,
    const uint64_t kmer_begin,
    const uint64_t kmer_end,
    uint16_t* const text_count_ptr,
    uint16_t* const pattern_count_ptr) {
  const uint16_t text_count = *text_count_ptr;
  const uint16_t pattern_count = *pattern_count_ptr;
  // KMER_COUNTING_INC_COUNT(text_count,pattern_count,chunk->curr_text_kmers,chunk->max_text_kmers);
  if (gem_expect_false(pattern_count > 0 && text_count < pattern_count)) {
    ++(chunk->curr_text_kmers);
    chunk->max_text_kmers = MAX(chunk->max_text_kmers,chunk->curr_text_kmers);
  }
  ++(*text_count_ptr);
}
#define KMER_COUNTING_INCREMENT(chunk_idx) \
  /* Check Tiling boundaries & Update kmer counters */ \
  if (text_tiles[chunk_idx].text_begin<=kmer_begin && kmer_end<text_tiles[chunk_idx].text_end) { \
    kmer_counting_increment_kmer_count( \
        text_tiles+chunk_idx,kmer_begin,kmer_end, \
        text_count_ptr+chunk_idx,pattern_count_ptr+chunk_idx); \
  }
/* Compute k-mer offsets & fetch counters*/
#define KMER_COUNTING_IN_FETCH_COUNTERS(text_count_ptr,pattern_count_ptr) \
  const uint8_t char_end = KMER_COUNTING_FILTER_CHAR(text[kmer_end]); \
  KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,char_end); \
  const uint64_t kmer_offset = kmer_idx*num_tiles; \
  uint16_t* const text_count_ptr = kmer_count_text + kmer_offset; \
  uint16_t* const pattern_count_ptr = kmer_count_pattern + kmer_offset
/* Load/Store window counter pointers */
#ifdef KMER_COUNTING_SLIDING_WINDOW
  #define KMER_COUNTING_WINDOW_IN(text_count_ptr,pattern_count_ptr) \
    kmer_window[kmer_window_in++] = text_count_ptr; \
    kmer_window[kmer_window_in++] = pattern_count_ptr; \
    if (kmer_window_in >= kmer_window_length) kmer_window_in = 0
  #define KMER_COUNTING_WINDOW_OUT(text_count_ptr,pattern_count_ptr) \
    uint16_t* const text_count_ptr = kmer_window[kmer_window_out++]; \
    uint16_t* const pattern_count_ptr = kmer_window[kmer_window_out++]; \
    if (kmer_window_out >= kmer_window_length) kmer_window_out = 0
#else
  #define KMER_COUNTING_WINDOW_IN(text_count_ptr,pattern_count_ptr)
  #define KMER_COUNTING_WINDOW_OUT(text_count_ptr,pattern_count_ptr)
#endif
/*
 * Kernel specializations
 */
uint64_t __attribute__ ((noinline)) kmer_counting_min_bound_nway_ks1(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t sliding_window_length,
    uint16_t** const kmer_window,
    const uint64_t kmer_window_length) {
  // Parameters
  const uint64_t kmer_length = kmer_counting->kmer_length;
  kmer_counting_text_tile_t* const text_tiles = kmer_counting->text_tiles;
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
#ifdef KMER_COUNTING_SLIDING_WINDOW
  uint64_t kmer_window_in = 0, kmer_window_out = 0;
#endif
  uint64_t kmer_idx = 0, kmer_end, kmer_begin;
  // Initial fill (kmer)
  for (kmer_end=0;kmer_end<kmer_length-1;++kmer_end) {
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,KMER_COUNTING_FILTER_CHAR(text[kmer_end]));
  }
  // Sliding window
  for (kmer_begin=0;kmer_end<text_length;++kmer_begin,++kmer_end) {
#ifdef KMER_COUNTING_SLIDING_WINDOW
    // OUT sliding window
    if (kmer_begin >= sliding_window_length) {
      // Fetch counts
      uint16_t* const text_count_ptr = kmer_window[kmer_window_out++];
      uint16_t* const pattern_count_ptr = kmer_window[kmer_window_out++];
      if (kmer_window_out >= kmer_window_length) kmer_window_out = 0;
      // Decrement kmer-count
      kmer_counting_decrement_kmer_count(
          text_tiles,kmer_begin,kmer_end,sliding_window_length,
          text_count_ptr,pattern_count_ptr);
    }
#endif
    // IN sliding window
    {
      // Fetch counters & store them in window
      const uint8_t char_end = KMER_COUNTING_FILTER_CHAR(text[kmer_end]);
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,char_end);
      const uint64_t kmer_offset = kmer_idx;
      uint16_t* const text_count_ptr = kmer_count_text + kmer_offset;
      uint16_t* const pattern_count_ptr = kmer_count_pattern + kmer_offset;
      KMER_COUNTING_WINDOW_IN(text_count_ptr,pattern_count_ptr);
      // Increment kmer counts
      kmer_counting_increment_kmer_count(
          text_tiles,kmer_begin,kmer_end,
          text_count_ptr,pattern_count_ptr);
    }
  }
  // Compute min-error bound
  return (text_tiles->num_key_kmers-text_tiles->max_text_kmers)/kmer_length;;
}
uint64_t __attribute__ ((noinline)) kmer_counting_min_bound_nway_ks2(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t sliding_window_length,
    uint16_t** const kmer_window,
    const uint64_t kmer_window_length) {
  // Parameters
  const uint64_t kmer_length = kmer_counting->kmer_length;
  const uint64_t num_tiles = kmer_counting->num_tiles;
  kmer_counting_text_tile_t* const text_tiles = kmer_counting->text_tiles;
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
#ifdef KMER_COUNTING_SLIDING_WINDOW
  uint64_t kmer_window_in = 0, kmer_window_out = 0;
#endif
  uint64_t kmer_idx = 0, kmer_end, kmer_begin;
  // Initial fill (kmer)
  for (kmer_end=0;kmer_end<kmer_length-1;++kmer_end) {
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,KMER_COUNTING_FILTER_CHAR(text[kmer_end]));
  }
  // Sliding window
  for (kmer_begin=0;kmer_end<text_length;++kmer_begin,++kmer_end) {
#ifdef KMER_COUNTING_SLIDING_WINDOW
    // OUT sliding window
    if (kmer_begin >= sliding_window_length) {
      /* Load counters from window */
      KMER_COUNTING_WINDOW_OUT(text_count_ptr,pattern_count_ptr);
      // Check Tiling boundaries & Update kmer counters
      KMER_COUNTING_DECREMENT(0);
      KMER_COUNTING_DECREMENT(1);
    }
#endif
    // IN sliding window
    {
      // Fetch counters & store them in window
      KMER_COUNTING_IN_FETCH_COUNTERS(text_count_ptr,pattern_count_ptr);
      KMER_COUNTING_WINDOW_IN(text_count_ptr,pattern_count_ptr);
      // Check Tiling boundaries & Update kmer counters
      KMER_COUNTING_INCREMENT(0);
      KMER_COUNTING_INCREMENT(1);
    }
  }
  // Compute min-error bound
  uint64_t min_error = 0;
  min_error += (text_tiles[0].num_key_kmers-text_tiles[0].max_text_kmers)/kmer_length;
  min_error += (text_tiles[1].num_key_kmers-text_tiles[1].max_text_kmers)/kmer_length;
  return min_error;
}
uint64_t __attribute__ ((noinline)) kmer_counting_min_bound_nway_ks3(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t sliding_window_length,
    uint16_t** const kmer_window,
    const uint64_t kmer_window_length) {
  // Parameters
  const uint64_t kmer_length = kmer_counting->kmer_length;
  const uint64_t num_tiles = kmer_counting->num_tiles;
  kmer_counting_text_tile_t* const text_tiles = kmer_counting->text_tiles;
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
#ifdef KMER_COUNTING_SLIDING_WINDOW
  uint64_t kmer_window_in = 0, kmer_window_out = 0;
#endif
  uint64_t kmer_idx = 0, kmer_end, kmer_begin;
  // Initial fill (kmer)
  for (kmer_end=0;kmer_end<kmer_length-1;++kmer_end) {
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,KMER_COUNTING_FILTER_CHAR(text[kmer_end]));
  }
  // Sliding window
  for (kmer_begin=0;kmer_end<text_length;++kmer_begin,++kmer_end) {
#ifdef KMER_COUNTING_SLIDING_WINDOW
    // OUT sliding window
    if (kmer_begin >= sliding_window_length) {
      /* Load counters from window */
      KMER_COUNTING_WINDOW_OUT(text_count_ptr,pattern_count_ptr);
      // Check Tiling boundaries & Update kmer counters
      KMER_COUNTING_DECREMENT(0);
      KMER_COUNTING_DECREMENT(1);
      KMER_COUNTING_DECREMENT(2);
    }
#endif
    // IN sliding window
    {
      // Fetch counters & store them in window
      KMER_COUNTING_IN_FETCH_COUNTERS(text_count_ptr,pattern_count_ptr);
      KMER_COUNTING_WINDOW_IN(text_count_ptr,pattern_count_ptr);
      // Check Tiling boundaries & Update kmer counters
      KMER_COUNTING_INCREMENT(0);
      KMER_COUNTING_INCREMENT(1);
      KMER_COUNTING_INCREMENT(2);
    }
  }
  // Compute min-error bound
  uint64_t min_error = 0;
  min_error += (text_tiles[0].num_key_kmers-text_tiles[0].max_text_kmers)/kmer_length;
  min_error += (text_tiles[1].num_key_kmers-text_tiles[1].max_text_kmers)/kmer_length;
  min_error += (text_tiles[2].num_key_kmers-text_tiles[2].max_text_kmers)/kmer_length;
  return min_error;
}
uint64_t __attribute__ ((noinline)) kmer_counting_min_bound_nway_ks4(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t sliding_window_length,
    uint16_t** const kmer_window,
    const uint64_t kmer_window_length) {
  // Parameters
  const uint64_t kmer_length = kmer_counting->kmer_length;
  const uint64_t num_tiles = kmer_counting->num_tiles;
  kmer_counting_text_tile_t* const text_tiles = kmer_counting->text_tiles;
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
#ifdef KMER_COUNTING_SLIDING_WINDOW
  uint64_t kmer_window_in = 0, kmer_window_out = 0;
#endif
  uint64_t kmer_idx = 0, kmer_end, kmer_begin;
  // Initial fill (kmer)
  for (kmer_end=0;kmer_end<kmer_length-1;++kmer_end) {
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,KMER_COUNTING_FILTER_CHAR(text[kmer_end]));
  }
  // Sliding window
  for (kmer_begin=0;kmer_end<text_length;++kmer_begin,++kmer_end) {
#ifdef KMER_COUNTING_SLIDING_WINDOW
    // OUT sliding window
    if (kmer_begin >= sliding_window_length) {
      /* Load counters from window */
      KMER_COUNTING_WINDOW_OUT(text_count_ptr,pattern_count_ptr);
      // Check Tiling boundaries & Update kmer counters
      KMER_COUNTING_DECREMENT(0);
      KMER_COUNTING_DECREMENT(1);
      KMER_COUNTING_DECREMENT(2);
      KMER_COUNTING_DECREMENT(3);
    }
#endif
    // IN sliding window
    {
      // Fetch counters & store them in window
      KMER_COUNTING_IN_FETCH_COUNTERS(text_count_ptr,pattern_count_ptr);
      KMER_COUNTING_WINDOW_IN(text_count_ptr,pattern_count_ptr);
      // Check Tiling boundaries & Update kmer counters
      KMER_COUNTING_INCREMENT(0);
      KMER_COUNTING_INCREMENT(1);
      KMER_COUNTING_INCREMENT(2);
      KMER_COUNTING_INCREMENT(3);
    }
  }
  // Compute min-error bound
  uint64_t min_error = 0;
  min_error += (text_tiles[0].num_key_kmers-text_tiles[0].max_text_kmers)/kmer_length;
  min_error += (text_tiles[1].num_key_kmers-text_tiles[1].max_text_kmers)/kmer_length;
  min_error += (text_tiles[2].num_key_kmers-text_tiles[2].max_text_kmers)/kmer_length;
  min_error += (text_tiles[3].num_key_kmers-text_tiles[3].max_text_kmers)/kmer_length;
  return min_error;
}
/*
 * Kernel specialization N-WAY
 */
uint64_t __attribute__ ((noinline)) kmer_counting_min_bound_nway_ksn(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t sliding_window_length,
    uint16_t** const kmer_window,
    const uint64_t kmer_window_length) {
  // Parameters
  const uint64_t kmer_length = kmer_counting->kmer_length;
  const uint64_t num_tiles = kmer_counting->num_tiles;
  kmer_counting_text_tile_t* const text_tiles = kmer_counting->text_tiles;
  uint16_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  uint16_t* const kmer_count_text = kmer_counting->kmer_count_text;
#ifdef KMER_COUNTING_SLIDING_WINDOW
  uint64_t kmer_window_in = 0, kmer_window_out = 0;
#endif
  uint64_t kmer_idx = 0, kmer_end, kmer_begin, chunk_idx;
  // Initial fill (kmer)
  for (kmer_end=0;kmer_end<kmer_length-1;++kmer_end) {
    const uint8_t enc_char = KMER_COUNTING_FILTER_CHAR(text[kmer_end]);
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char); // Update kmer-index
  }
  // Sliding window
  for (kmer_begin=0;kmer_end<text_length;++kmer_begin,++kmer_end) {
#ifdef KMER_COUNTING_SLIDING_WINDOW
    // OUT sliding window
    if (kmer_begin >= sliding_window_length) {
      // Fetch counts
      uint16_t* const text_count_ptr = kmer_window[kmer_window_out++];
      uint16_t* const pattern_count_ptr = kmer_window[kmer_window_out++];
      if (kmer_window_out >= kmer_window_length) kmer_window_out = 0;
      // Decrement kmer-count
      for (chunk_idx=0;chunk_idx<num_tiles;++chunk_idx) {
        // Check Tiling boundaries & Update kmer counters
        if (text_tiles[chunk_idx].text_begin+sliding_window_length > kmer_begin) break;
        if (kmer_end >= text_tiles[chunk_idx].text_end) continue;
        kmer_counting_decrement_kmer_count(
            text_tiles+chunk_idx,kmer_begin,kmer_end,sliding_window_length,
            text_count_ptr+chunk_idx,pattern_count_ptr+chunk_idx);
      }
    }
#endif
    // IN sliding window
    {
      // Fetch counters & store them in window
      const uint8_t char_end = KMER_COUNTING_FILTER_CHAR(text[kmer_end]);
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,char_end);
      const uint64_t kmer_offset = kmer_idx*num_tiles;
      uint16_t* const text_count_ptr = kmer_count_text + kmer_offset;
      uint16_t* const pattern_count_ptr = kmer_count_pattern + kmer_offset;
      KMER_COUNTING_WINDOW_IN(text_count_ptr,pattern_count_ptr);
      // Increment kmer counts
      for (chunk_idx=0;chunk_idx<num_tiles;++chunk_idx) {
        // Check Tiling boundaries & Update kmer counters
        if (text_tiles[chunk_idx].text_begin > kmer_begin) break;
        if (kmer_end >= text_tiles[chunk_idx].text_end) continue;
        kmer_counting_increment_kmer_count(
            text_tiles+chunk_idx,kmer_begin,kmer_end,
            text_count_ptr+chunk_idx,pattern_count_ptr+chunk_idx);
      }
    }
  }
  // Compute min-error bound
  uint64_t min_error = 0;
  for (chunk_idx=0;chunk_idx<num_tiles;++chunk_idx) {
    min_error += (text_tiles[chunk_idx].num_key_kmers-text_tiles[chunk_idx].max_text_kmers)/kmer_length;
  }
  return min_error;
}
/*
 * K-mer counting compute min-error bound
 *   Init + Switch
 */
uint64_t kmer_counting_min_bound_nway(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error,
    mm_allocator_t* const mm_allocator) {
  PROFILE_START(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
  // Check enabled
  if (!kmer_counting->enabled) {
    PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_NA);
    return 0; // Don't filter
  }
  // Prepare filter
  if (kmer_counting->kmer_count_text==NULL) return 0; // Don't filter
  const uint64_t num_tiles = kmer_counting->num_tiles;
  memset(kmer_counting->kmer_count_text,0,kmer_counting->num_kmers*num_tiles*sizeof(uint16_t));
  // Prepare tiling
  kmer_counting_prepare_tiling(kmer_counting,text_length,max_error);
  // MM Push
  mm_allocator_push_state(mm_allocator);
  // Prepare sliding window (pointers to the profile counters going out the sliding window)
  const uint64_t sliding_window_length = MIN(text_length,kmer_counting->sliding_window_length);
  const uint64_t kmer_window_length = 2*sliding_window_length;
  uint16_t** const kmer_window = mm_allocator_calloc(mm_allocator,kmer_window_length,uint16_t*,false);
  // Switch kernel specialization
  uint64_t min_error;
  switch (num_tiles) {
    case 1:
      min_error = kmer_counting_min_bound_nway_ks1(kmer_counting,
          text,text_length,sliding_window_length,kmer_window,kmer_window_length);
      break;
    case 2:
      min_error = kmer_counting_min_bound_nway_ks2(kmer_counting,
          text,text_length,sliding_window_length,kmer_window,kmer_window_length);
      break;
    case 3:
      min_error = kmer_counting_min_bound_nway_ks3(kmer_counting,
          text,text_length,sliding_window_length,kmer_window,kmer_window_length);
      break;
    case 4:
      min_error = kmer_counting_min_bound_nway_ks4(kmer_counting,
          text,text_length,sliding_window_length,kmer_window,kmer_window_length);
      break;
    default:
      min_error = kmer_counting_min_bound_nway_ksn(kmer_counting,
          text,text_length,sliding_window_length,kmer_window,kmer_window_length);
      break;
  }
  // Free
  mm_allocator_pop_state(mm_allocator);
  // Return min-error
  PROFILE_STOP(GP_FC_KMER_COUNTER_FILTER,PROFILE_LEVEL);
  return min_error;
}
