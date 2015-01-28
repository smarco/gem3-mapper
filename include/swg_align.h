/*
 * PROJECT: GEMMapper
 * FILE: swg_align.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef SWG_ALIGN_H_
#define SWG_ALIGN_H_

#include "essentials.h"

/*
 * SWG Query Profile
 */
typedef struct {
  int32_t M; // Alignment matching/mismatching
  int32_t I; // Alignment ends with a gap in the reference (insertion)
  int32_t D; // Alignment ends with a gap in the read (deletion)
} swg_cell_t;
typedef int32_t matching_score_t[DNA__N_RANGE][DNA__N_RANGE];
typedef struct {
  int32_t gap_open_penalty;
  int32_t gap_extension_penalty;
  int32_t generic_match_score;
  int32_t generic_mismatch_score;
  matching_score_t matching_score;
} swg_penalties_t;
typedef struct {
  uint8_t* query_profile_uint8[DNA__N_RANGE];
  int16_t* query_profile_int16[DNA__N_RANGE];
  uint8_t weight_bias;
} swg_query_profile_t;

/*
 * Check CIGAR string
 */
GEM_INLINE bool align_check_match(
    FILE* const stream,const uint8_t* const key,const uint64_t key_length,const uint8_t* const text,
    const uint64_t text_length,vector_t* const cigar_buffer,uint64_t const cigar_offset,
    uint64_t const cigar_length,const bool verbose);

/*
 * Levenshtein Alignment
 */
GEM_INLINE int64_t align_levenshtein_get_distance(
    const char* const key,const uint64_t key_length,
    const char* const text,const uint64_t text_length,
    const bool ends_free,uint64_t* const position);

/*
 * Compile SWG Query Profile
 */
GEM_INLINE void swg_compile_query_profile(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,const uint64_t key_length,mm_stack_t* const mm_stack);

/*
 * Smith-waterman-gotoh Alignment
 */
GEM_INLINE void swg_align_match(
    const uint8_t* const key,const uint64_t key_length,const bool* const allowed_enc,
    const swg_penalties_t* swg_penalties,uint64_t* const match_position,uint8_t* const text,
    const uint64_t text_length,uint64_t max_bandwidth,const bool begin_free,const bool end_free,
    vector_t* const cigar_buffer,uint64_t* const cigar_length,int64_t* const effective_length,
    int32_t* const alignment_score,mm_stack_t* const mm_stack);

/*
 * TMP
 */
GEM_INLINE void swg_align_match_full_32b(
    const uint8_t* const key,const uint64_t key_length,const swg_penalties_t* swg_penalties,
    uint64_t* const match_position,uint8_t* const text,const uint64_t text_length,
    vector_t* const cigar_buffer,uint64_t* const cigar_length,int64_t* const effective_length,
    int32_t* const alignment_score,mm_stack_t* const mm_stack);
GEM_INLINE void swg_align_match_banded_32b(
    const uint8_t* const key,const uint64_t key_length,const swg_penalties_t* swg_penalties,
    uint64_t* const match_position,uint8_t* const text,const uint64_t text_length,
    uint64_t max_bandwidth,const bool begin_free,const bool end_free,
    vector_t* const cigar_buffer,uint64_t* const cigar_length,int64_t* const effective_length,
    int32_t* const alignment_score,mm_stack_t* const mm_stack);

#endif /* SWG_ALIGN_H_ */
