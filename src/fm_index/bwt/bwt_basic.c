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
 *   BWT data structure (basic implementation not-sampled).
 *   Provides basic routines to encode a DNA text into a
 *   Burrows-Wheeler transform using a compact bit representation
 *   and counter buckets as to enhance Occ/rank queries
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

#include "text/dna_text.h"
#include "fm_index/bwt/bwt_basic.h"

/*
 * BWT Model & Version
 */
#define BWT_BASIC_MODEL_NO  2000ull

/*
 * BWT Dimensions
 *   Alphabet = 3bits => 8 characters dense
 *   MinorBlock => 16B (minorCounters) + 24B (BM) => 40B
 *   MayorBlock => 40KB (1024 MinorBlocks)
 *   |MinorCounters| = 8*UINT16
 *   |MayorCounters| = 8*UINT64  (Stored separately)
 *   MinorBlock.BM = 3*UINT64
 *     - Bitmap-LayOut Bitwise (3 x 64word)
 *     - [64bits-Bit0-{LSB}][64bits-Bit1][64bits-Bit2-{MSB}]
 * Remarks::
 *   Number of MayorBlocks (3GB genome) => 3*(1024^3)/(1024×64b) = 49152 MBlocks
 *   Assuming SysPage = 4KB =>
 *       1 SysPage => 102 minorBlocks (1 unaligned, will cause 2 Pagefaults)
 *       5 SysPage => 512 minorBlocks (4 unaligned) => 1,0078125 PageFaults per access [OK]
 *       40KB/4KB = 10 SysPage per MayorBlock
 */
#define BWT_MINOR_COUNTER_RANGE         8
#define BWT_MINOR_BLOCK_COUNTER_LENGTH 16
#define BWT_MINOR_BLOCK_LENGTH         64
#define BWT_MINOR_BLOCK_BITS            3

#define BWT_MINOR_BLOCK_SIZE \
  (BWT_MINOR_BLOCK_COUNTER_LENGTH*BWT_MINOR_COUNTER_RANGE + \
   BWT_MINOR_BLOCK_LENGTH*BWT_MINOR_BLOCK_BITS)/8                                        /* 40B */
#define BWT_MINOR_BLOCK_64WORDS (BWT_MINOR_BLOCK_SIZE/8)                                 /* Offset 5 */

#define BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS (1<<10)                                         /* (((2^16)−1)÷64)+1 = 1024 */
#define BWT_MAYOR_BLOCK_LENGTH (BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS*BWT_MINOR_BLOCK_LENGTH) /* 1024*64 = 65536 */
#define BWT_MAYOR_BLOCK_SIZE   (BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS*BWT_MINOR_BLOCK_SIZE)   /* 1024*40B = 40KB */

#define BWT_MAYOR_COUNTER_SIZE  UINT64_SIZE
#define BWT_MAYOR_COUNTER_RANGE 8

/*
 * BWT Builder Auxiliary functions
 */
void bwt_basic_builder_initialize(
    bwt_basic_builder_t* const bwt_builder,
    const uint64_t bwt_text_length,
    const uint64_t* const character_occurrences) {
  bwt_basic_t* const bwt = &bwt_builder->bwt;
  /*
   * Meta-Data
   */
  bwt->length = bwt_text_length;
  // Set characters count
  uint8_t i;
  bwt->c = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  for (i=0;i<BWT_MAYOR_COUNTER_RANGE;++i) bwt->c[i] = character_occurrences[i];
  // Set the cumulative occurrences (the relative ordering of the symbols is encoded)
  bwt->C = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  bwt->C[ENC_DNA_CHAR_A]    = 0;
  bwt->C[ENC_DNA_CHAR_C]    = bwt->C[ENC_DNA_CHAR_A] + bwt->c[ENC_DNA_CHAR_A];
  bwt->C[ENC_DNA_CHAR_G]    = bwt->C[ENC_DNA_CHAR_C] + bwt->c[ENC_DNA_CHAR_C];
  bwt->C[ENC_DNA_CHAR_T]    = bwt->C[ENC_DNA_CHAR_G] + bwt->c[ENC_DNA_CHAR_G];
  bwt->C[ENC_DNA_CHAR_N]    = bwt->C[ENC_DNA_CHAR_T] + bwt->c[ENC_DNA_CHAR_T];
  bwt->C[ENC_DNA_CHAR_SEP]  = bwt->C[ENC_DNA_CHAR_N] + bwt->c[ENC_DNA_CHAR_N];
  bwt->C[ENC_DNA_CHAR_JUMP] = bwt->C[ENC_DNA_CHAR_SEP] + bwt->c[ENC_DNA_CHAR_SEP];
  bwt->num_minor_blocks = DIV_CEIL(bwt->length,BWT_MINOR_BLOCK_LENGTH);
  bwt->num_mayor_blocks = DIV_CEIL(bwt->num_minor_blocks,BWT_MAYOR_BLOCK_NUM_MINOR_BLOCKS);
  /*
   * Allocate BWT Structures
   */
  // BWT Buckets-Structure
  bwt->mayor_counters = mm_calloc(bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  bwt->bwt_mem = (uint64_t*) mm_calloc(bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE,uint8_t,true);
  // Auxiliary variables
  bwt_builder->mayor_counter = 0;
  bwt_builder->minor_block_mem = bwt->bwt_mem;
  bwt_builder->mayor_occ = mm_calloc(BWT_MAYOR_COUNTER_RANGE,uint64_t,true);
  bwt_builder->minor_occ = mm_calloc(BWT_MINOR_COUNTER_RANGE,uint64_t,true);
}
#define bwt_basic_builder_write_enc(enc,bit_mask,layer_0,layer_1,layer_2) \
  switch (enc) { \
    case ENC_DNA_CHAR_A:   /* 000 */ \
      /* layer_0 |= bit_mask; */ \
      /* layer_1 |= bit_mask; */ \
      /* layer_2 |= bit_mask; */ \
      break; \
    case ENC_DNA_CHAR_C:   /* 001 */ \
      layer_0 |= bit_mask; \
      /* layer_1 |= bit_mask; */ \
      /* layer_2 |= bit_mask; */ \
      break; \
    case ENC_DNA_CHAR_G:   /* 010 */ \
      /* layer_0 |= bit_mask; */ \
      layer_1 |= bit_mask; \
      /* layer_2 |= bit_mask; */ \
      break; \
    case ENC_DNA_CHAR_T:   /* 011 */ \
      layer_0 |= bit_mask; \
      layer_1 |= bit_mask; \
      /* layer_2 |= bit_mask; */ \
      break; \
    case ENC_DNA_CHAR_N:   /* 100 */ \
      /* layer_0 |= bit_mask; */ \
      /* layer_1 |= bit_mask; */ \
      layer_2 |= bit_mask; \
      break; \
    case ENC_DNA_CHAR_SEP: /* 101 */ \
      layer_0 |= bit_mask; \
      /* layer_1 |= bit_mask; */ \
      layer_2 |= bit_mask; \
      break; \
    case ENC_DNA_CHAR_JUMP: /* 110 */ \
      /* layer_0 |= bit_mask; */ \
      layer_1 |= bit_mask; \
      layer_2 |= bit_mask; \
      break; \
    default: \
      GEM_INVALID_CASE(); \
      break; \
  } \
  bit_mask <<= 1

void bwt_basic_builder_write_mayor_counters(bwt_basic_builder_t* const bwt_builder) {
  uint64_t* const mayor_counters = bwt_builder->bwt.mayor_counters + bwt_builder->mayor_counter;
  uint64_t i;
  for (i=0;i<BWT_MAYOR_COUNTER_RANGE;++i) {
    // Update sentinel counters
    bwt_builder->mayor_occ[i] += bwt_builder->minor_occ[i];
    bwt_builder->minor_occ[i] = 0; // Reset
    // Dump counters
    mayor_counters[i] = bwt_builder->mayor_occ[i] + bwt_builder->bwt.C[i];
  }
  // Update MayorCounter position
  bwt_builder->mayor_counter += BWT_MAYOR_COUNTER_RANGE;
}
void bwt_basic_builder_write_minor_counters(bwt_basic_builder_t* const bwt_builder) {
  uint16_t* const minor_counters_mem = (uint16_t*) bwt_builder->minor_block_mem;
  uint64_t i;
  // Dump counters
  for (i=0;i<BWT_MINOR_COUNTER_RANGE;++i) {
    minor_counters_mem[i] = bwt_builder->minor_occ[i]; // Write minorCounter
  }
  // Update minorBlock pointer (Skip counters)
  bwt_builder->minor_block_mem += 2; // 2*64bits == 8*16bits
}
void bwt_basic_builder_write_minor_block(
    bwt_basic_builder_t* const bwt_builder,
    const uint64_t layer_0,
    const uint64_t layer_1,
    const uint64_t layer_2) {
  *bwt_builder->minor_block_mem = layer_0; ++(bwt_builder->minor_block_mem);
  *bwt_builder->minor_block_mem = layer_1; ++(bwt_builder->minor_block_mem);
  *bwt_builder->minor_block_mem = layer_2; ++(bwt_builder->minor_block_mem);
}
/*
 * BWT Builder
 */
bwt_basic_builder_t* bwt_basic_builder_new(
    dna_text_t* const bwt_text,
    const uint64_t* const character_occurrences,
    const bool verbose) {
  /*
   * Allocate & initialize builder
   */
  bwt_basic_builder_t* const bwt_builder = mm_alloc(bwt_basic_builder_t);
  const uint64_t bwt_length = dna_text_get_length(bwt_text);
  bwt_basic_builder_initialize(bwt_builder,bwt_length,character_occurrences);
  /*
   * Compute BWT & write
   */
  // Ticker
  ticker_t ticker;
  ticker_percentage_reset(&ticker,verbose,"Building-BWT::Generating BWT-Bitmap",bwt_length,10,true);
  // Cursors/Locators
  uint64_t mayorBlock_pos = BWT_MAYOR_BLOCK_LENGTH;
  uint64_t minorBlock_pos = BWT_MINOR_BLOCK_LENGTH;
  // BM Layers
  uint64_t layer_0=0, layer_1=0, layer_2=0, bit_mask=UINT64_ONE_MASK;
  // Iterate over the BWT
  const uint8_t* const bwt = dna_text_get_text(bwt_text);
  uint64_t bwt_pos = 0;
  while (bwt_pos < bwt_length) {
    // Get BWT character
    const uint8_t enc = bwt[bwt_pos++];
    // Print MayorCounters
    if (mayorBlock_pos==BWT_MAYOR_BLOCK_LENGTH) {
      ticker_update(&ticker,mayorBlock_pos); // Update ticker
      bwt_basic_builder_write_mayor_counters(bwt_builder);
      mayorBlock_pos = 0;
    }
    // Print MinorCounters
    if (minorBlock_pos==BWT_MINOR_BLOCK_LENGTH) {
      bwt_basic_builder_write_minor_counters(bwt_builder);
      minorBlock_pos = 0;
    }
    // Write character to BitMap
    bwt_basic_builder_write_enc(enc,bit_mask,layer_0,layer_1,layer_2);
    // Account current char & Next
    ++(bwt_builder->minor_occ[enc]);
    ++minorBlock_pos; ++mayorBlock_pos;
    // Close MinorBlock (dump)
    if (minorBlock_pos==BWT_MINOR_BLOCK_LENGTH) {
      bwt_basic_builder_write_minor_block(bwt_builder,layer_0,layer_1,layer_2);
      // Reset
      layer_0=0; layer_1=0; layer_2=0;
      bit_mask = UINT64_ONE_MASK;
    }
  }
  // Close last MinorBlock (dump)
  if (minorBlock_pos!=BWT_MINOR_BLOCK_LENGTH) {
    bwt_basic_builder_write_minor_block(bwt_builder,layer_0,layer_1,layer_2);
  }
  // Ticker
  ticker_finish(&ticker);
  return bwt_builder;
}
void bwt_basic_builder_write(
    fm_t* const file_manager,
    bwt_basic_builder_t* const bwt_builder) {
  BWT_BUILDER_CHECK(bwt_builder);
  /* Meta-Data */
  fm_write_uint64(file_manager,BWT_BASIC_MODEL_NO); // Model Number
  fm_write_uint64(file_manager,bwt_builder->bwt.length); // Length of the BWT
  fm_write_mem(file_manager,bwt_builder->bwt.c,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // Occurrences of each character
  fm_write_mem(file_manager,bwt_builder->bwt.C,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // The cumulative occurrences
  /* BWT Buckets-Structure */
  fm_write_uint64(file_manager,bwt_builder->bwt.num_minor_blocks); // Total number of minor blocks
  fm_write_uint64(file_manager,bwt_builder->bwt.num_mayor_blocks); // Total number of mayor blocks
  /* Mayor Counters [Aligned to 4KB] */
  fm_skip_align_4KB(file_manager);
  fm_write_mem(file_manager,bwt_builder->bwt.mayor_counters,
      bwt_builder->bwt.num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*BWT_MAYOR_COUNTER_SIZE);
  /* BWT Structure [Aligned to 4KB] */
  fm_skip_align_4KB(file_manager);
  fm_write_mem(file_manager,bwt_builder->bwt.bwt_mem,bwt_builder->bwt.num_minor_blocks*BWT_MINOR_BLOCK_SIZE);
}
void bwt_basic_builder_delete(bwt_basic_builder_t* const bwt_builder) {
  BWT_BUILDER_CHECK(bwt_builder);
  mm_free(bwt_builder->bwt.c);
  mm_free(bwt_builder->bwt.C);
  mm_free(bwt_builder->bwt.mayor_counters);
  mm_free(bwt_builder->bwt.bwt_mem);
  mm_free(bwt_builder->mayor_occ);
  mm_free(bwt_builder->minor_occ);
  mm_free(bwt_builder);
}
/*
 * BWT Loader
 */
bwt_basic_t* bwt_basic_read_mem(mm_t* const memory_manager,const bool check) {
  // Allocate handler
  bwt_basic_t* const bwt = mm_alloc(bwt_basic_t);
  /* Meta-Data */
  const uint64_t bwt_serial_no = mm_read_uint64(memory_manager); // Serial Number
  gem_cond_error(bwt_serial_no!=BWT_BASIC_MODEL_NO,BWT_WRONG_MODEL_NO,bwt_serial_no,(uint64_t)BWT_BASIC_MODEL_NO);
  bwt->length = mm_read_uint64(memory_manager); // Length of the BWT
  bwt->c = mm_read_mem(memory_manager,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // Occurrences of each character
  bwt->C = mm_read_mem(memory_manager,BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE); // The cumulative occurrences
  /* BWT Buckets-Structure */
  bwt->num_minor_blocks = mm_read_uint64(memory_manager); // Total number of minor blocks
  bwt->num_mayor_blocks = mm_read_uint64(memory_manager); // Total number of mayor blocks
  /* Mayor Counters [Aligned to 4KB] */
  mm_skip_align_4KB(memory_manager);
  bwt->mm_counters = NULL;
  bwt->mayor_counters = mm_read_mem(memory_manager,bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*BWT_MAYOR_COUNTER_SIZE);
  /* BWT Structure [Aligned to 4KB] */
  mm_skip_align_4KB(memory_manager);
  bwt->mm_bwt_mem = NULL;
  bwt->bwt_mem = mm_read_mem(memory_manager,bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE);
  // Return
  return bwt;
}
void bwt_basic_delete(bwt_basic_t* const bwt) {
  BWT_CHECK(bwt);
  // Free counters memory
  if (bwt->mm_counters) {
    mm_free(bwt->c);
    mm_free(bwt->c);
    mm_bulk_free(bwt->mm_counters);
  }
  // Free bwt memory
  if (bwt->mm_bwt_mem) mm_bulk_free(bwt->mm_bwt_mem);
  // Free handler
  mm_free(bwt);
}
/*
 * BWT General Accessors
 */
uint64_t bwt_basic_get_size(bwt_basic_t* const bwt) {
  // Compute sizes
  const uint64_t minor_blocks_size = bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE;
  const uint64_t mayor_counters_size = bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE;
  return (BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE) +   /* bwt->c */
         (BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE) +   /* bwt->C */
         mayor_counters_size +                     /* bwt->mayor_counters */
         minor_blocks_size;                        /* bwt->bwt_mem */
}
uint64_t bwt_basic_builder_get_length(const bwt_basic_builder_t* const bwt_builder) {
  BWT_BUILDER_CHECK(bwt_builder);
  return bwt_basic_get_length(&bwt_builder->bwt);
}
uint64_t bwt_basic_get_length(const bwt_basic_t* const bwt) {
  BWT_CHECK(bwt);
  return bwt->length;
}
uint64_t bwt_basic_builder_get_size(bwt_basic_builder_t* const bwt_builder) {
  BWT_BUILDER_CHECK(bwt_builder);
  return bwt_basic_get_size(&bwt_builder->bwt);
}
bool bwt_basic_is_same_bucket(
    const uint64_t lo,
    const uint64_t hi) {
  return (lo/BWT_MINOR_BLOCK_LENGTH) == (hi/BWT_MINOR_BLOCK_LENGTH);
}
/* Gets the mayor_counters & minor_block corresponding to position @i */
#define BWT_LOCATE_BLOCK(bwt,position,block_pos,block_mod,mayor_counters,block_mem) \
  const uint64_t block_pos = position / BWT_MINOR_BLOCK_LENGTH; \
  const uint64_t block_mod = position % BWT_MINOR_BLOCK_LENGTH; \
  const uint64_t* const mayor_counters = bwt->mayor_counters + (position/BWT_MAYOR_BLOCK_LENGTH)*BWT_MAYOR_COUNTER_RANGE; \
  const uint64_t* const block_mem = bwt->bwt_mem + block_pos*BWT_MINOR_BLOCK_64WORDS
/* Gets the mayor_counters & minor_block corresponding to position @i */
void bwt_basic_get_block_location(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    bwt_block_locator_t* const block_loc) {
  block_loc->block_pos = position / BWT_MINOR_BLOCK_LENGTH;
  block_loc->block_mod = position % BWT_MINOR_BLOCK_LENGTH;
  block_loc->mayor_counters = bwt->mayor_counters + (position/BWT_MAYOR_BLOCK_LENGTH)*BWT_MAYOR_COUNTER_RANGE;
  block_loc->block_mem = bwt->bwt_mem + block_loc->block_pos*BWT_MINOR_BLOCK_64WORDS;
}
/* Computes and returns the encoded letter */
uint8_t bwt_basic_char_(
    const uint64_t block_mod,
    const uint64_t* const block_mem) {
  const uint64_t letter_mask = 1ull << block_mod;
  const uint8_t bit_1 = ((block_mem[2] & letter_mask) != 0);
  const uint8_t bit_2 = ((block_mem[3] & letter_mask) != 0);
  const uint8_t bit_3 = ((block_mem[4] & letter_mask) != 0);
  return (bit_3 << 2) | (bit_2 << 1) | bit_1;
}
/* Computes and returns the rank */
uint64_t bwt_basic_erank_(
    const uint8_t char_enc,
    const uint64_t block_mod,
    const uint64_t* const mayor_counters,
    const uint64_t* const block_mem) {
  // Calculate the exclusive rank for the given DNA character
  const uint64_t sum_counters = mayor_counters[char_enc] + ((uint16_t*)block_mem)[char_enc];
  const uint64_t bitmap = (block_mem[2]^xor_table_3[char_enc]) &
                          (block_mem[3]^xor_table_2[char_enc]) &
                          (block_mem[4]^xor_table_1[char_enc]);
  // Return rank
  return sum_counters + POPCOUNT_64(bitmap & uint64_erank_mask(block_mod));
}
/* Computes and returns the rank of an interval located in the same block */
void bwt_basic_erank_interval_(
    const uint8_t char_enc,
    const uint64_t lo_value,
    const uint64_t block_mod,
    const uint64_t* const mayor_counters,
    const uint64_t* const block_mem,
    uint64_t* const lo,
    uint64_t* const hi) {
  // Fetching Regular DNA Characters
  const uint64_t sum_counters = mayor_counters[char_enc] + ((uint16_t*)block_mem)[char_enc];
  const uint64_t bitmap = (block_mem[2]^xor_table_3[char_enc]) &
                          (block_mem[3]^xor_table_2[char_enc]) &
                          (block_mem[4]^xor_table_1[char_enc]);
  *hi = sum_counters +
        ( POPCOUNT_64(bitmap & uint64_erank_mask(block_mod)) );
  *lo = sum_counters +
        ( POPCOUNT_64(bitmap & uint64_erank_mask(lo_value % BWT_MINOR_BLOCK_LENGTH)) );
}
/* Pre-computes the block elements (faster computation of all possible ranks) */
void bwt_basic_precompute_(
    const uint64_t block_mod,
    const uint64_t* const block_mem,
    bwt_block_elms_t* const bwt_block_elms) {
  const uint64_t erase_mask = uint64_erank_mask(block_mod);
  const uint64_t bitmap_1 = block_mem[2] & erase_mask;
  const uint64_t bitmap_N1 = (block_mem[2]^(-1ull)) & erase_mask;
  const uint64_t bitmap_2 = block_mem[3];
  const uint64_t bitmap_N2 = bitmap_2^(-1ull);
  const uint64_t bitmap_3 = block_mem[4];
  bwt_block_elms->bitmap_1__2[0] = bitmap_N2 & bitmap_N1;
  bwt_block_elms->bitmap_1__2[1] = bitmap_N2 & bitmap_1;
  bwt_block_elms->bitmap_1__2[2] = bitmap_2  & bitmap_N1;
  bwt_block_elms->bitmap_1__2[3] = bitmap_2  & bitmap_1;
  bwt_block_elms->bitmap_3[0] = (bitmap_3^(-1ull));
  bwt_block_elms->bitmap_3[1] = bitmap_3;
}
/*
 * BWT Character Accessors
 */
uint8_t bwt_basic_char(
    const bwt_basic_t* const bwt,
    const uint64_t position) {
  BWT_CHAR_TICK();
  /* Locate Block */
  const uint64_t block_pos = position / BWT_MINOR_BLOCK_LENGTH;
  const uint64_t block_mod = position % BWT_MINOR_BLOCK_LENGTH;
  const uint64_t* const block_mem = bwt->bwt_mem + block_pos*BWT_MINOR_BLOCK_64WORDS;
  return bwt_basic_char_(block_mod,block_mem);
}
char bwt_basic_char_character(
    const bwt_basic_t* const bwt,
    const uint64_t position) {
  return dna_decode(bwt_basic_char(bwt,position));
}
/*
 * BWT ERank (Exclusive Rank Function)
 */
uint64_t bwt_basic_builder_erank(
    const bwt_basic_builder_t* const bwt_builder,
    const uint8_t char_enc,
    const uint64_t position) {
  BWT_ERANK_TICK();
  BWT_LOCATE_BLOCK((&bwt_builder->bwt),position,block_pos,block_mod,mayor_counters,block_mem);
  return bwt_basic_erank_(char_enc,block_mod,mayor_counters,block_mem);
}
uint64_t bwt_basic_erank(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t position) {
  BWT_ERANK_TICK();
  BWT_LOCATE_BLOCK(bwt,position,block_pos,block_mod,mayor_counters,block_mem);
  return bwt_basic_erank_(char_enc,block_mod,mayor_counters,block_mem);
}
uint64_t bwt_basic_erank_character(
    const bwt_basic_t* const bwt,
    const char character,
    const uint64_t position) {
  return bwt_basic_erank(bwt,dna_encode(character),position);
}
void bwt_basic_erank_interval(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out) {
  BWT_ERANK_INTERVAL_TICK();
  BWT_LOCATE_BLOCK(bwt,hi_in,block_pos,block_mod,mayor_counters,block_mem);
  bwt_basic_erank_interval_(char_enc,lo_in,block_mod,mayor_counters,block_mem,lo_out,hi_out);
}
/*
 * BWT Prefetched ERank
 */
void bwt_basic_prefetch(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    bwt_block_locator_t* const block_loc) {
  BWT_PREFETCH_TICK();
  bwt_basic_get_block_location(bwt,position,block_loc);
  BWT_PREFETCH_BLOCK(block_loc);
}
uint64_t bwt_basic_prefetched_erank(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t position,
    const bwt_block_locator_t* const block_loc) {
  BWT_ERANK_TICK();
  return bwt_basic_erank_(char_enc,block_loc->block_mod,block_loc->mayor_counters,block_loc->block_mem);
}
void bwt_basic_prefetched_erank_interval(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const uint64_t lo_in,
    const uint64_t hi_in,
    uint64_t* const lo_out,
    uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc) {
  BWT_ERANK_INTERVAL_TICK();
  bwt_basic_erank_interval_(char_enc,lo_in,block_loc->block_mod,
      block_loc->mayor_counters,block_loc->block_mem,lo_out,hi_out);
}
/*
 *  BWT Precomputed ERank (Precomputation of the block's elements)
 */
void bwt_basic_precompute(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms) {
  BWT_PRECOMPUTE_TICK();
  bwt_basic_get_block_location(bwt,position,block_loc);
  bwt_basic_precompute_(block_loc->block_mod,block_loc->block_mem,block_elms);
}
void bwt_basic_precompute_interval(
    const bwt_basic_t* const bwt,
    const uint64_t lo,
    const uint64_t hi,
    bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms) {
  bwt_basic_precompute(bwt,hi,block_loc,block_elms);
  block_elms->gap_mask = uint64_erank_inv_mask(lo % BWT_MINOR_BLOCK_LENGTH);
}
void bwt_basic_prefetched_precompute(
    const bwt_basic_t* const bwt,
    const bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms) {
  BWT_PRECOMPUTE_TICK();
  bwt_basic_precompute_(block_loc->block_mod,block_loc->block_mem,block_elms);
}
void bwt_basic_prefetched_precompute_interval(
    const bwt_basic_t* const bwt,
    const uint64_t lo,
    const bwt_block_locator_t* const block_loc,
    bwt_block_elms_t* const block_elms) {
  bwt_basic_prefetched_precompute(bwt,block_loc,block_elms);
  block_elms->gap_mask = uint64_erank_inv_mask(lo % BWT_MINOR_BLOCK_LENGTH);
}
uint64_t bwt_basic_precomputed_erank(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    const bwt_block_locator_t* const block_loc,
    const bwt_block_elms_t* const block_elms) {
  BWT_ERANK_TICK();
  // Return precomputed-erank
  const uint64_t sum_counters = block_loc->mayor_counters[char_enc] + ((uint16_t*)block_loc->block_mem)[char_enc];
  const uint64_t bitmap = block_elms->bitmap_1__2[char_enc & 3] & block_elms->bitmap_3[char_enc>>2];
  return sum_counters + POPCOUNT_64(bitmap);
}
void bwt_basic_precomputed_erank_interval(
    const bwt_basic_t* const bwt,
    const uint8_t char_enc,
    uint64_t* const lo_out,
    uint64_t* const hi_out,
    const bwt_block_locator_t* const block_loc,
    const bwt_block_elms_t* const block_elms) {
  BWT_ERANK_INTERVAL_TICK();
  const uint64_t bitmap = block_elms->bitmap_1__2[char_enc & 3] & block_elms->bitmap_3[char_enc>>2];
  const uint64_t bitmap_gap = bitmap & block_elms->gap_mask;
  if (gem_expect_false(bitmap_gap==0)) {
    *hi_out = 0;
    *lo_out = 0;
  } else {
    *hi_out = block_loc->mayor_counters[char_enc] + ((uint16_t*)block_loc->block_mem)[char_enc] + POPCOUNT_64(bitmap);
    *lo_out = *hi_out - POPCOUNT_64(bitmap_gap);
  }
}
/*
 * BWT LF (Last to first)
 */
uint64_t bwt_basic_LF(
    const bwt_basic_t* const bwt,
    const uint64_t position) {
  uint8_t char_enc;
  return bwt_basic_LF__enc(bwt,position,&char_enc);
}
uint64_t bwt_basic_prefetched_LF(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    const bwt_block_locator_t* const block_loc) {
  uint8_t char_enc;
  return bwt_basic_prefetched_LF__enc(bwt,position,&char_enc,block_loc);
}
uint64_t bwt_basic_LF__enc(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    uint8_t* const char_enc) {
  BWT_LF_TICK();
  BWT_LOCATE_BLOCK(bwt,position,block_pos,block_mod,mayor_counters,block_mem);
  *char_enc = bwt_basic_char_(block_mod,block_mem);  // Retrieve char_enc
  return bwt_basic_erank_(*char_enc,block_mod,mayor_counters,block_mem);
}
uint64_t bwt_basic_LF__character(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    char* const character) {
  uint8_t char_enc = 0;
  const uint64_t rank_LF = bwt_basic_LF__enc(bwt,position,&char_enc);
  *character = dna_decode(char_enc);
  return rank_LF;
}
uint64_t bwt_basic_prefetched_LF__enc(
    const bwt_basic_t* const bwt,
    const uint64_t position,
    uint8_t* const char_enc,
    const bwt_block_locator_t* const block_loc) {
  BWT_LF_TICK();
  *char_enc = bwt_basic_char_(block_loc->block_mod,block_loc->block_mem); // Retrieve char_enc
  return bwt_basic_erank_(*char_enc,block_loc->block_mod,block_loc->mayor_counters,block_loc->block_mem);
}
/*
 * Display
 */
void bwt_basic_print_(
    FILE* const stream,
    bwt_basic_t* const bwt) {
  // Compute sizes
  const uint64_t mayor_counters_size = bwt->num_mayor_blocks*BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE;
  const uint64_t minor_blocks_size = bwt->num_minor_blocks*BWT_MINOR_BLOCK_SIZE;
  const uint64_t bwt_basic_total_size =
      (BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE) /* bwt->c */ +
      (BWT_MAYOR_COUNTER_RANGE*UINT64_SIZE) /* bwt->C */ +
      mayor_counters_size + minor_blocks_size;
  const uint64_t minor_counters_size = bwt->num_minor_blocks*BWT_MINOR_COUNTER_RANGE*BWT_MINOR_BLOCK_COUNTER_LENGTH/8;
  const uint64_t minor_bitmap_size = bwt->num_minor_blocks*BWT_MINOR_BLOCK_LENGTH*BWT_MINOR_BLOCK_BITS/8;
  const uint64_t minor_sampling_size = bwt->num_minor_blocks*UINT64_SIZE;
  // Display BWT
  tab_fprintf(stream,"[GEM]>BWT\n");
  tab_fprintf(stream,"  => Architecture\tGBWT.l2.s1.8xC64.8xc16.3xbm64\n");
  tab_fprintf(stream,"    => Buckets 2-levels\n");
  tab_fprintf(stream,"    => Step    1-step\n");
  tab_fprintf(stream,"    => MayorCounters.length 8x64bits\n");
  tab_fprintf(stream,"    => MinorCounters.length 8x16bits\n");
  tab_fprintf(stream,"    => Bitmap.Bitwise 3x64bits\n");
  tab_fprintf(stream,"    => Sampled.Positions\n");
  tab_fprintf(stream,"      => Sampling.MayorCounters 64bits\n");
  tab_fprintf(stream,"      => Sampling.MinorCounters 16bits\n");
  tab_fprintf(stream,"      => Sampling.Bitmap        64bits\n");
  tab_fprintf(stream,"  => Total.length %"PRIu64"\n",bwt->length);
  tab_fprintf(stream,"  => Total.Size %"PRIu64" MB\n",CONVERT_B_TO_MB(bwt_basic_total_size));
  tab_fprintf(stream,"    => Mayor.counters %"PRIu64" (%"PRIu64" MB) [%2.3f%%]\n",
      bwt->num_mayor_blocks,CONVERT_B_TO_MB(mayor_counters_size),
      PERCENTAGE(mayor_counters_size,bwt_basic_total_size));
  tab_fprintf(stream,"    => Minor.Blocks %"PRIu64" (%"PRIu64" MB) [%2.3f%%]\n",
      bwt->num_minor_blocks,CONVERT_B_TO_MB(minor_blocks_size),
      PERCENTAGE(minor_blocks_size,bwt_basic_total_size));
  tab_fprintf(stream,"      => Minor.Counters (%"PRIu64" MB) [%2.3f%%]\n",
      CONVERT_B_TO_MB(minor_counters_size),
      PERCENTAGE(minor_counters_size,bwt_basic_total_size));
  tab_fprintf(stream,"      => Minor.Bitmap (%"PRIu64" MB) [%2.3f%%]\n",
      CONVERT_B_TO_MB(minor_bitmap_size),
      PERCENTAGE(minor_bitmap_size,bwt_basic_total_size));
  tab_fprintf(stream,"      => Minor.Sampling (%"PRIu64" MB) [%2.3f%%]\n",
      CONVERT_B_TO_MB(minor_sampling_size),
      PERCENTAGE(minor_sampling_size,bwt_basic_total_size));
  tab_fprintf(stream,"  => Occurrences\tA=[%"PRIu64"] C=[%"PRIu64"] G=[%"PRIu64"] T=[%"PRIu64"] N=[%"PRIu64"] |=[%"PRIu64"]\n",
      bwt->c[ENC_DNA_CHAR_A],bwt->c[ENC_DNA_CHAR_C],bwt->c[ENC_DNA_CHAR_G],
      bwt->c[ENC_DNA_CHAR_T],bwt->c[ENC_DNA_CHAR_N],bwt->c[ENC_DNA_CHAR_SEP]);
  tab_fprintf(stream,"  => Cumulative.Occ\tA=[%"PRIu64"] C=[%"PRIu64"] G=[%"PRIu64"] T=[%"PRIu64"] N=[%"PRIu64"] |=[%"PRIu64"]\n",
      bwt->C[ENC_DNA_CHAR_A],bwt->C[ENC_DNA_CHAR_C],bwt->C[ENC_DNA_CHAR_G],
      bwt->C[ENC_DNA_CHAR_T],bwt->C[ENC_DNA_CHAR_N],bwt->C[ENC_DNA_CHAR_SEP]);
  /*
   * Stats // todo
   */
  // Flush
  fflush(stream);
}
void bwt_basic_builder_print(
    FILE* const stream,
    bwt_basic_builder_t* const bwt_builder) {
  bwt_basic_print_(stream,&bwt_builder->bwt);
}
void bwt_basic_print(
    FILE* const stream,
    bwt_basic_t* const bwt) {
  bwt_basic_print_(stream,bwt);
}

