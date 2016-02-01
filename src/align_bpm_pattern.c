/*
 * PROJECT: GEMMapper
 * FILE: align_bpm_pattern.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "align_bpm_pattern.h"
#include "gpu_buffer_align_bpm.h"

/*
 * Checks
 */
#define BPM_PATTERN_CHECK(bpm_pattern) GEM_CHECK_NULL(bpm_pattern->PEQ)

/*
 * Constants
 */
#define BMP_W64_LENGTH UINT64_LENGTH
#define BMP_W64_ONES   UINT64_MAX
#define BMP_W64_MASK   (1ull<<63)

/*
 * Compile Pattern
 */
bpm_pattern_t* bpm_pattern_compile(
    uint8_t* const pattern,const uint64_t pattern_length,
    const uint64_t max_error,mm_stack_t* const mm_stack) {
  // Alloc
  bpm_pattern_t* const bpm_pattern = mm_stack_alloc(mm_stack,bpm_pattern_t);
  // Calculate dimensions
  const uint64_t pattern_num_words = DIV_CEIL(pattern_length,BPM_ALIGN_WORD_LENGTH);
  const uint64_t PEQ_length = pattern_num_words*BPM_ALIGN_WORD_LENGTH;
  const uint64_t pattern_mod = pattern_length%BPM_ALIGN_WORD_LENGTH;
  // Init fields
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words = pattern_num_words;
  bpm_pattern->pattern_mod = pattern_mod;
  // Allocate memory
  const uint64_t aux_vector_size = pattern_num_words*BPM_ALIGN_WORD_SIZE;
  const uint64_t PEQ_size = DNA__N_RANGE*aux_vector_size;
  const uint64_t score_size = pattern_num_words*UINT64_SIZE;
  bpm_pattern->PEQ = mm_stack_malloc(mm_stack,PEQ_size);
  bpm_pattern->P = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->M = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->level_mask = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->score = mm_stack_malloc(mm_stack,score_size);
  bpm_pattern->init_score = mm_stack_malloc(mm_stack,score_size);
  bpm_pattern->pattern_left = mm_stack_malloc(mm_stack,(pattern_num_words+1)*UINT64_SIZE);
  // Init PEQ
  memset(bpm_pattern->PEQ,0,PEQ_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = pattern[i];
    const uint64_t block = i/BPM_ALIGN_WORD_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_ALIGN_WORD_LENGTH);
    bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,enc_char)] |= mask;
  }
  for (;i<PEQ_length;++i) { // Padding
    const uint64_t block = i/BPM_ALIGN_WORD_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_ALIGN_WORD_LENGTH);
    uint64_t j;
    for (j=0;j<DNA__N_RANGE;++j) {
      bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,j)] |= mask;
    }
  }
  // Init auxiliary data
  uint64_t pattern_left = pattern_length;
  const uint64_t top = pattern_num_words-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i] = BMP_W64_MASK;
    bpm_pattern->init_score[i] = BPM_ALIGN_WORD_LENGTH;
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_ALIGN_WORD_LENGTH) ? pattern_left-BPM_ALIGN_WORD_LENGTH : 0;
  }
  for (;i<=pattern_num_words;++i) {
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_ALIGN_WORD_LENGTH) ? pattern_left-BPM_ALIGN_WORD_LENGTH : 0;
  }
  if (pattern_mod>0) {
    const uint64_t mask_shift = pattern_mod-1;
    bpm_pattern->level_mask[top] = 1ull<<(mask_shift);
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top] = BMP_W64_MASK;
    bpm_pattern->init_score[top] = BPM_ALIGN_WORD_LENGTH;
  }
  // GPU Compile
  gpu_bpm_pattern_compile(bpm_pattern,GPU_WORDS128_PER_TILE,max_error);
  // Return
  return bpm_pattern;
}
/*
 * Compile Pattern Tiles
 */
bpm_pattern_t* bpm_pattern_compile_tiles(
    bpm_pattern_t* const bpm_pattern,uint64_t words64_per_tile,
    const uint64_t max_error,mm_stack_t* const mm_stack) {
  // Set tile tall (64-bits words per tile)
  const uint64_t min_words64_per_tile = DIV_CEIL(max_error,UINT64_LENGTH);
  if (min_words64_per_tile > words64_per_tile) words64_per_tile = min_words64_per_tile;
  const uint64_t num_tiles = DIV_CEIL(bpm_pattern->pattern_num_words,words64_per_tile);
  // Alloc
  bpm_pattern_t* const bpm_pattern_tiles = mm_stack_calloc(mm_stack,num_tiles,bpm_pattern_t,false);
  // Init Tiles (Store meta-data in the first tile)
  bpm_pattern_tiles->words_per_tile = words64_per_tile;
  bpm_pattern_tiles->num_pattern_tiles = num_tiles;
  // Init BPM chunks
  const uint64_t tile_tall = words64_per_tile * UINT64_LENGTH;
  uint64_t pattern_length_left = bpm_pattern->pattern_length;
  uint64_t offset_words = 0, i;
  for (i=0;i<num_tiles;++i) {
    // Initialize pattern-chunk variables
    bpm_pattern_t* const bpm_pattern_tile = bpm_pattern_tiles + i;
    bpm_pattern_tile->pattern_length =  MIN(tile_tall,pattern_length_left);
    bpm_pattern_tile->pattern_num_words = words64_per_tile;
    bpm_pattern_tile->pattern_mod = bpm_pattern_tile->pattern_length % BPM_ALIGN_WORD_LENGTH;
    bpm_pattern_tile->PEQ = bpm_pattern->PEQ + offset_words*DNA__N_RANGE;
    bpm_pattern_tile->P = bpm_pattern->P;
    bpm_pattern_tile->M = bpm_pattern->M;
    bpm_pattern_tile->level_mask = bpm_pattern->level_mask + offset_words;
    bpm_pattern_tile->score = bpm_pattern->score + offset_words;
    bpm_pattern_tile->init_score = bpm_pattern->init_score + offset_words;
    bpm_pattern_tile->pattern_left = NULL;
    offset_words += words64_per_tile;
    pattern_length_left -= tile_tall;
  }
  // Return
  return bpm_pattern_tiles;
}
