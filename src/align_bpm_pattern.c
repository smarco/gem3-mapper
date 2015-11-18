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
GEM_INLINE void bpm_pattern_compile_chunks(
    bpm_pattern_t* const bpm_pattern,uint8_t* const pattern,
    const uint64_t pattern_length,const uint64_t max_error,mm_stack_t* const mm_stack) {
  // Init BPM chunks
  const uint64_t words_per_chunk = DIV_CEIL(max_error,BPM_ALIGN_WORD_LENGTH);
  const uint64_t num_chunks = DIV_CEIL(bpm_pattern->pattern_num_words,words_per_chunk);
  bpm_pattern->words_per_chunk = words_per_chunk;
  bpm_pattern->num_pattern_chunks = num_chunks;
  bpm_pattern->bpm_pattern_chunks = mm_stack_calloc(mm_stack,num_chunks,bpm_pattern_t,false);
  uint64_t offset_words = 0, i;
  for (i=0;i<num_chunks;++i) {
    // Initialize pattern-chunk variables
    bpm_pattern_t* const bpm_pattern_chunk = bpm_pattern->bpm_pattern_chunks + i;
    bpm_pattern_chunk->P = bpm_pattern->P;
    bpm_pattern_chunk->M = bpm_pattern->M;
    bpm_pattern_chunk->PEQ = bpm_pattern->PEQ + offset_words*DNA__N_RANGE;
    bpm_pattern_chunk->pattern_num_words = words_per_chunk;
    bpm_pattern_chunk->level_mask = bpm_pattern->level_mask + offset_words;
    bpm_pattern_chunk->score = bpm_pattern->score + offset_words;
    bpm_pattern_chunk->init_score = bpm_pattern->init_score + offset_words;
    bpm_pattern_chunk->pattern_left = NULL;
    offset_words += words_per_chunk;
  }
}
GEM_INLINE void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,uint8_t* const pattern,
    const uint64_t pattern_length,const uint64_t max_error,mm_stack_t* const mm_stack) {
  GEM_CHECK_NULL(bpm_pattern);
  GEM_CHECK_NULL(pattern);
  GEM_CHECK_ZERO(pattern_length);
  // Calculate dimensions
  const uint64_t word_length = BPM_ALIGN_WORD_LENGTH;
  const uint64_t word_size = BPM_ALIGN_WORD_SIZE;
  const uint64_t pattern_num_words = DIV_CEIL(pattern_length,BPM_ALIGN_WORD_LENGTH);
  const uint64_t PEQ_length = pattern_num_words*word_length;
  const uint64_t pattern_mod = pattern_length%word_length;
  // Init fields
  bpm_pattern->pattern_word_size = word_size;
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words = pattern_num_words;
  bpm_pattern->pattern_mod = pattern_mod;
  bpm_pattern->PEQ_length = PEQ_length;
  // Allocate memory
  const uint64_t aux_vector_size = pattern_num_words*word_size;
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
  for (;i<PEQ_length;++i) {
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
    bpm_pattern->init_score[i] = word_length;
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > word_length) ? pattern_left-word_length : 0;
  }
  for (;i<=pattern_num_words;++i) {
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > word_length) ? pattern_left-word_length : 0;
  }
  if (pattern_mod>0) {
    const uint64_t mask_shift = pattern_mod-1;
    bpm_pattern->level_mask[top] = 1ull<<(mask_shift);
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top] = BMP_W64_MASK;
    bpm_pattern->init_score[top] = word_length;
  }
#ifdef BPM_TILED
  // Init BPM chunks
  bpm_pattern_compile_chunks(bpm_pattern,pattern,pattern_length,max_error,mm_stack);
#endif
}
