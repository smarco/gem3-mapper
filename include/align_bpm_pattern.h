/*
 * PROJECT: GEMMapper
 * FILE: align_bpm_pattern.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ALIGN_BPM_PATTERN_H_
#define ALIGN_BPM_PATTERN_H_

#include "essentials.h"

/*
 * Constants
 */
#define BPM_ALIGN_WORD_LENGTH  UINT64_LENGTH
#define BPM_ALIGN_WORD_SIZE    UINT64_SIZE

/*
 * BPM Pattern
 */
typedef struct _bpm_pattern_t bpm_pattern_t;
struct _bpm_pattern_t {
  /* BMP Pattern */
  uint64_t* PEQ;              // Pattern equalities (Bit vector for Myers-DP)
  uint64_t pattern_word_size; // Word size (in bytes)
  uint64_t pattern_length;    // Length
  uint64_t pattern_num_words; // ceil(Length / |w|)
  uint64_t pattern_mod;       // Length % |w|
  uint64_t PEQ_length;        // ceil(Length / |w|) * |w|
  /* BPM Auxiliary data */
  uint64_t* P;
  uint64_t* M;
  uint64_t* level_mask;
  int64_t* score;
  int64_t* init_score;
  uint64_t* pattern_left;
  /* BPM tiles (Pattern split in tiles) */
  uint64_t words_per_tile;
  uint64_t num_pattern_tiles;
  /* BPM-GPU Dimensions */
  uint64_t gpu_num_entries;
  uint64_t gpu_num_tiles;
  uint64_t gpu_entries_per_tile;
  bpm_pattern_t* bpm_pattern_tiles;
};

/*
 * Pattern Accessors
 */
#define BPM_PATTERN_PEQ_IDX(word_pos,encoded_character)   ((word_pos*DNA__N_RANGE)+(encoded_character))
#define BPM_PATTERN_BDP_IDX(position,num_words,word_pos)  ((position)*(num_words)+(word_pos))

/*
 * Compile Pattern
 */
void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,uint8_t* const pattern,
    const uint64_t pattern_length,const uint64_t max_error,
    mm_stack_t* const mm_stack);

#endif /* ALIGN_BPM_PATTERN_H_ */
