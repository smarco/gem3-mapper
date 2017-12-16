/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 * DESCRIPTION:
 */

#include "gpu/gpu_buffer_bpm_pattern.h"
#include "align/align_bpm_distance.h"
#include "text/dna_text.h"

/*
 * Compile/Decompile Query Entries [DTO]
 */
void gpu_buffer_bpm_pattern_compile(
    gpu_bpm_peq_entry_t* const pattern_entry,
    const bpm_pattern_t* const bpm_pattern) {
  // Copy PEQ pattern
  const uint64_t bpm_pattern_num_words = bpm_pattern->pattern_num_words64;
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  const uint64_t gpu_pattern_num_words = DIV_CEIL(pattern_length,GPU_BPM_ENTRY_LENGTH)*2;
  const uint32_t* PEQ = (uint32_t*) bpm_pattern->PEQ;
  uint64_t i, entry=0, subentry=0;
  for (i=0;i<bpm_pattern_num_words;++i) {
    // Update location
    if (subentry==GPU_BPM_NUM_SUB_ENTRIES) {
      subentry = 0; ++entry;
    }
    // Copy pattern
    uint8_t enc_char;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      pattern_entry[entry].bitmap[enc_char][subentry] = *PEQ; ++PEQ;
      pattern_entry[entry].bitmap[enc_char][subentry+1] = *PEQ; ++PEQ;
    }
    subentry += 2;
  }
  for (;i<gpu_pattern_num_words;++i) {
    // Update location
    if (subentry==GPU_BPM_NUM_SUB_ENTRIES) {
      subentry = 0; ++entry;
    }
    // Copy pattern
    uint8_t enc_char;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      pattern_entry[entry].bitmap[enc_char][subentry] = UINT32_ONES;
      pattern_entry[entry].bitmap[enc_char][subentry+1] = UINT32_ONES;
    }
    subentry += 2;
  }
}
void gpu_buffer_bpm_pattern_decompile(
    gpu_bpm_peq_entry_t* const pattern_entry,
    const uint32_t pattern_length,
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator) {
  // Calculate dimensions
  const uint64_t pattern_num_words = DIV_CEIL(pattern_length,UINT64_LENGTH);
  const uint64_t pattern_mod = pattern_length%UINT64_LENGTH;
  // Init fields
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words64 = pattern_num_words;
  bpm_pattern->pattern_mod = pattern_mod;
  // Allocate memory
  const uint64_t gpu_pattern_num_words = DIV_CEIL(pattern_length,GPU_BPM_ENTRY_LENGTH);
  const uint64_t aux_vector_size = gpu_pattern_num_words*GPU_BPM_ENTRY_SIZE;
  const uint64_t PEQ_size = DNA__N_RANGE*aux_vector_size;
  const uint64_t score_size = pattern_num_words*UINT64_SIZE;
  uint64_t* PEQ = mm_allocator_malloc(mm_allocator,PEQ_size);
  bpm_pattern->PEQ = PEQ;
  bpm_pattern->P = mm_allocator_malloc(mm_allocator,aux_vector_size);
  bpm_pattern->M = mm_allocator_malloc(mm_allocator,aux_vector_size);
  bpm_pattern->level_mask = mm_allocator_malloc(mm_allocator,aux_vector_size);
  bpm_pattern->score = mm_allocator_malloc(mm_allocator,score_size);
  bpm_pattern->init_score = mm_allocator_malloc(mm_allocator,score_size);
  bpm_pattern->pattern_left = mm_allocator_malloc(mm_allocator,(pattern_num_words+1)*UINT64_SIZE);
  // Init PEQ
  uint8_t enc_char;
  uint64_t i, entry=0, subentry=0;
  for (i=0;i<gpu_pattern_num_words;++i) {
    subentry = 0;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      *PEQ = ((uint64_t)pattern_entry[entry].bitmap[enc_char][subentry]) |
             ((uint64_t)pattern_entry[entry].bitmap[enc_char][subentry+1] << 32);
      ++PEQ;
    }
    subentry += 2;
    for (enc_char=0;enc_char<DNA__N_RANGE;++enc_char) {
      *PEQ = ((uint64_t)pattern_entry[entry].bitmap[enc_char][subentry]) |
             ((uint64_t)pattern_entry[entry].bitmap[enc_char][subentry+1] << 32);
      ++PEQ;
    }
    ++entry;
  }
  // Init auxiliary data
  uint64_t pattern_left = pattern_length;
  const uint64_t top = pattern_num_words-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i] = BMP_W64_MASK;
    bpm_pattern->init_score[i] = UINT64_LENGTH;
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > UINT64_LENGTH) ? pattern_left-UINT64_LENGTH : 0;
  }
  for (;i<=pattern_num_words;++i) {
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > UINT64_LENGTH) ? pattern_left-UINT64_LENGTH : 0;
  }
  if (pattern_mod>0) {
    const uint64_t mask_shift = pattern_mod-1;
    bpm_pattern->level_mask[top] = 1ull<<(mask_shift);
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top] = BMP_W64_MASK;
    bpm_pattern->init_score[top] = UINT64_LENGTH;
  }
}
