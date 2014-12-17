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
#include "../resources/myers_gpu/myers-interface.h"

/*
 * Constants
 */
#define ALIGN_DISTANCE_INF     UINT32_MAX
#define ALIGN_COLUMN_INF       UINT64_MAX
#define BPM_ALIGN_WORD_LENGTH  UINT64_LENGTH
#define BPM_ALIGN_WORD_SIZE    UINT64_SIZE

/*
 * GPU Constants
 */
#define BPM_GPU_PATTERN_NUM_SUB_ENTRIES  BMP_GPU_PEQ_SUBENTRIES
#define BPM_GPU_PATTERN_ENTRY_LENGTH     BMP_GPU_PEQ_ENTRY_LENGTH
#define BPM_GPU_PATTERN_SUBENTRY_LENGTH  BMP_GPU_PEQ_SUBENTRY_LENGTH
#define BPM_GPU_PATTERN_ENTRY_SIZE       (BPM_GPU_PATTERN_ENTRY_LENGTH/UINT8_SIZE)
#define BPM_GPU_PATTERN_ALPHABET_LENGTH  BMP_GPU_PEQ_ALPHABET_SIZE

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
  /* BPM chunks (Pattern split in chunks) */
  uint64_t words_per_chunk;
  uint64_t num_pattern_chunks;
  /* BPM-GPU Dimensions */
  uint64_t gpu_num_entries;
  uint64_t gpu_entries_per_chunk;
  uint64_t gpu_num_chunks;
  bpm_pattern_t* bpm_pattern_chunks;
};

/*
 * Compile Pattern
 */
GEM_INLINE void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,uint8_t* const pattern,
    const uint64_t pattern_length,const uint64_t max_error,mm_stack_t* const mm_stack);

/*
 * Bit-compressed Alignment BMP (BitParalellMyers) - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
GEM_INLINE bool bpm_get_distance(
    bpm_pattern_t* const bpm_pattern,const uint8_t* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance);
GEM_INLINE bool bpm_get_distance__cutoff(
    const bpm_pattern_t* const bpm_pattern,const uint8_t* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance,const bool quick_abandon);
GEM_INLINE void bpm_align_match(
    const uint8_t* const key,const bpm_pattern_t* const bpm_pattern,
    uint64_t* const match_position,const uint8_t* const sequence,
    const uint64_t sequence_length,const uint64_t max_distance,
    vector_t* const cigar_vector,uint64_t* const cigar_vector_offset,uint64_t* const cigar_length,
    uint64_t* const distance,int64_t* const effective_length,mm_stack_t* const mm_stack);

/*
 * BMP Tiled (bound)
 */
GEM_INLINE void bpm_bound_distance_tiled(
    bpm_pattern_t* const bpm_pattern,const uint8_t* const sequence,const uint64_t sequence_length,
    uint64_t* const levenshtein_distance,uint64_t* const levenshtein_match_pos,const uint64_t max_error);

#endif /* BPM_ALIGN_H_ */
