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

/*
 * BPM Pattern
 */
typedef struct {
  /* BMP Pattern */
  uint8_t* PEQ;               // Pattern equalities (Bit vector for Myers-DP)
  uint64_t pattern_word_size; // Word size (in bytes)
  uint64_t pattern_length;    // Length
  uint64_t pattern_num_words; // ceil(Length / |w|)
  uint64_t pattern_mod;       // Length % |w|
  uint64_t PEQ_length;        // ceil(Length / |w|) * |w|
  /* BPM Auxiliary data */
  uint8_t* P;
  uint8_t* M;
  uint8_t* level_mask;
  int64_t* score;
  int64_t* init_score;
} bpm_pattern_t;

/*
 * (Re)alignment Basic: Dynamic Programming - LEVENSHTEIN/EDIT
 */
GEM_INLINE int64_t align_levenshtein_get_distance(
    const char* const pattern,const uint64_t pattern_length,
    const char* const sequence,const uint64_t sequence_length,
    const bool ends_free,uint64_t* const position);

/*
 * Pattern Accessors
 */
#define BPM_PATTERN_PEQ_IDX(encoded_character,word_pos,num_words) ((encoded_character)*(num_words)+(word_pos))
#define BPM_PATTERN_BDP_IDX(position,word_pos,num_blocks) ((position)*(num_blocks)+(word_pos))

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
    bpm_pattern_t* const bpm_pattern,const char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance);
GEM_INLINE bool bpm_get_distance__cutoff(
    const bpm_pattern_t* const bpm_pattern,const uint8_t* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance);
GEM_INLINE void bpm_align_match(
    const uint8_t* const key,const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const sequence,uint64_t* const match_position,
    const uint64_t matching_distance,const uint64_t matching_column,
    vector_t* const cigar_buffer,uint64_t* const cigar_buffer_offset,uint64_t* const cigar_length,
    mm_stack_t* const mm_stack);

#endif /* BPM_ALIGN_H_ */
