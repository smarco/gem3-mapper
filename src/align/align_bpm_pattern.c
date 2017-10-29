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

#include "gpu/gpu_buffer_bpm_distance.h"
#include "text/dna_text.h"
#include "align/align_bpm_pattern.h"

/*
 * Checks
 */
#define BPM_PATTERN_CHECK(bpm_pattern) GEM_CHECK_NULL(bpm_pattern->PEQ)

/*
 * Constants
 */
#define BPM_W64_LENGTH UINT64_LENGTH
#define BPM_W64_SIZE   UINT64_SIZE
#define BPM_W64_ONES   UINT64_MAX
#define BPM_W64_MASK   (1ull<<63)

/*
 * Compile Pattern
 */
void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,
    uint8_t* const pattern,
    const uint64_t pattern_length,
    mm_allocator_t* const mm_allocator) {
  // Calculate dimensions
  const uint64_t pattern_num_words64 = DIV_CEIL(pattern_length,BPM_W64_LENGTH);
  const uint64_t PEQ_length = pattern_num_words64*BPM_W64_LENGTH;
  const uint64_t pattern_mod = pattern_length%BPM_W64_LENGTH;
  // Init fields
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words64 = pattern_num_words64;
  bpm_pattern->pattern_mod = pattern_mod;
  // Allocate memory
  const uint64_t aux_vector_size = pattern_num_words64*BPM_W64_SIZE;
  const uint64_t PEQ_size = DNA__N_RANGE*aux_vector_size;
  const uint64_t score_size = pattern_num_words64*UINT64_SIZE;
  const uint64_t total_memory = PEQ_size + 3*aux_vector_size + 2*score_size + (pattern_num_words64+1)*UINT64_SIZE;
  void* memory = mm_allocator_malloc(mm_allocator,total_memory);
  bpm_pattern->PEQ = memory; memory += PEQ_size;
  bpm_pattern->P = memory; memory += aux_vector_size;
  bpm_pattern->M = memory; memory += aux_vector_size;
  bpm_pattern->level_mask = memory; memory += aux_vector_size;
  bpm_pattern->score = memory; memory += score_size;
  bpm_pattern->init_score = memory; memory += score_size;
  bpm_pattern->pattern_left = memory;
  // Init PEQ
  memset(bpm_pattern->PEQ,0,PEQ_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = pattern[i];
    if (enc_char==ENC_DNA_CHAR_N) continue; // N's Inequality
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,enc_char)] |= mask;
  }
  for (;i<PEQ_length;++i) { // Padding
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    uint64_t j;
    for (j=0;j<DNA__N_RANGE;++j) {
      bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,j)] |= mask;
    }
  }
  // Init auxiliary data
  uint64_t pattern_left = pattern_length;
  const uint64_t top = pattern_num_words64-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i] = BPM_W64_MASK;
    bpm_pattern->init_score[i] = BPM_W64_LENGTH;
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left-BPM_W64_LENGTH : 0;
  }
  for (;i<=pattern_num_words64;++i) {
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left-BPM_W64_LENGTH : 0;
  }
  if (pattern_mod > 0) {
    const uint64_t mask_shift = pattern_mod-1;
    bpm_pattern->level_mask[top] = 1ull<<(mask_shift);
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top] = BPM_W64_MASK;
    bpm_pattern->init_score[top] = BPM_W64_LENGTH;
  }
}
void bpm_pattern_destroy(
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,bpm_pattern->PEQ);
}
/*
 * Compile Pattern Tiles
 */
void bpm_pattern_compile_tiles(
    bpm_pattern_t* const bpm_pattern,
    const uint64_t offset_words64,
    const uint64_t tile_length,
    bpm_pattern_t* const bpm_pattern_tile) {
  // Initialize pattern-chunk variables (BPM chunks)
  bpm_pattern_tile->pattern_length = tile_length;
  bpm_pattern_tile->pattern_num_words64 = DIV_CEIL(tile_length,BPM_W64_LENGTH);
  bpm_pattern_tile->pattern_mod = tile_length % BPM_W64_LENGTH;
  bpm_pattern_tile->PEQ = bpm_pattern->PEQ + offset_words64*DNA__N_RANGE;
  bpm_pattern_tile->P = bpm_pattern->P;
  bpm_pattern_tile->M = bpm_pattern->M;
  bpm_pattern_tile->level_mask = bpm_pattern->level_mask + offset_words64;
  bpm_pattern_tile->score = bpm_pattern->score + offset_words64;
  bpm_pattern_tile->init_score = bpm_pattern->init_score + offset_words64;
  bpm_pattern_tile->pattern_left = NULL;
}
