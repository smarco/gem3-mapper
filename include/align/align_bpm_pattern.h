/*
 * PROJECT: GEMMapper
 * FILE: align_bpm_pattern.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ALIGN_BPM_PATTERN_H_
#define ALIGN_BPM_PATTERN_H_

#include "utils/essentials.h"

/*
 * BPM Pattern
 */
typedef struct _bpm_pattern_t bpm_pattern_t;
struct _bpm_pattern_t {
  /* BMP Pattern */
  uint64_t* PEQ;                // Pattern equalities (Bit vector for Myers-DP)
  uint64_t pattern_length;      // Length
  uint64_t pattern_num_words64; // ceil(Length / |w|)
  uint64_t pattern_mod;         // Length % |w|
  /* BPM Auxiliary data */
  uint64_t* P;
  uint64_t* M;
  uint64_t* level_mask;
  int64_t* score;
  int64_t* init_score;
  uint64_t* pattern_left;
  /* BPM tiles (Pattern split in tiles) */
  uint64_t num_pattern_tiles;  // Total number of tiles
  uint64_t tile_length;        // Ideal tile length (the last can be shorter)
};

/*
 * Pattern Accessors
 */
#define BPM_PATTERN_PEQ_IDX(word_pos,encoded_character)   ((word_pos*DNA__N_RANGE)+(encoded_character))
#define BPM_PATTERN_BDP_IDX(position,num_words,word_pos)  ((position)*(num_words)+(word_pos))

/*
 * Compile Pattern
 */
bpm_pattern_t* bpm_pattern_compile(
    uint8_t* const pattern,
    const uint64_t pattern_length,
    const uint64_t max_error,
    mm_stack_t* const mm_stack);

/*
 * Compile Pattern Tiles
 */
bpm_pattern_t* bpm_pattern_compile_tiles(
    bpm_pattern_t* const bpm_pattern,
    const uint64_t prefered_words64_per_tile,
    const uint64_t max_error,
    mm_stack_t* const mm_stack);

#endif /* ALIGN_BPM_PATTERN_H_ */
