/*
 * PROJECT: GEMMapper
 * FILE: bpm_align.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef BPM_ALIGN_H_
#define BPM_ALIGN_H_

#include "essentials.h"

typedef struct {
  /* BMP Pattern */
  uint8_t* peq;               // Pattern equalities (Bit vector for Myers-DP)
  uint64_t pattern_length;    // Length
  uint64_t pattern_num_words; // ceil(Length / |w|)
  uint64_t pattern_mod;       // Length % |w|
  uint64_t peq_length;        // ceil(Length / |w|) * |w|
  /* BPM Auxiliary data */
  uint8_t* P;
  uint8_t* M;
  uint8_t* level_mask;
  int64_t* score;
  int64_t* init_score;
} bpm_pattern_t;

/*
 * Compile Pattern
 */

GEM_INLINE void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,const uint64_t word_size,
    uint8_t* const pattern,const uint64_t pattern_length,mm_stack_t* const mm_stack);

/*
 * Bit-compressed Alignment
 *   BMP[BitParalellMyers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
GEM_INLINE bool bpm_get_distance(
    bpm_pattern_t* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance);
GEM_INLINE bool bpm_get_distance__cutoff(
    bpm_pattern_t* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance);

#endif /* BPM_ALIGN_H_ */
