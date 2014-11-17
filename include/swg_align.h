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
 * BPM Pattern
 */
#define SWG_INT_MIN INT16_MIN
typedef int32_t swg_int;
typedef struct {
  swg_int M; // Alignment matching/mismatching
  swg_int I; // Alignment ends with a gap in the reference (insertion)
  swg_int D; // Alignment ends with a gap in the read (deletion)
} swg_cell_t;
typedef struct {
  uint32_t matching_score;
  uint32_t mismatch_penalty;
  uint32_t gap_open_penalty;
  uint32_t gap_extension_penalty; // TODO adjust wide type
} swg_penalties_t;



/*
 * Utils
 */
#define ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_element_type,element_length) \
  if (cigar_buffer->type == cigar_element_type) { \
    cigar_buffer->length += element_length; \
  } else { \
    if (cigar_buffer->type!=cigar_null) ++(cigar_buffer); \
    cigar_buffer->type = cigar_element_type; \
    cigar_buffer->length = element_length; \
  }

/*
 * Levenshtein Distance
 */
GEM_INLINE int64_t align_levenshtein_get_distance(
    const char* const pattern,const uint64_t pattern_length,
    const char* const sequence,const uint64_t sequence_length,
    const bool ends_free,uint64_t* const position);

/*
 * Smith-waterman-gotoh
 */
GEM_INLINE void swg_align_match(
    const uint8_t* const key,const uint64_t key_length,const swg_penalties_t* swgp,
    uint64_t* const match_position,const uint8_t* const sequence,const uint64_t sequence_length,
    const uint64_t max_distance,vector_t* const cigar_vector,uint64_t* const cigar_vector_offset,
    uint64_t* const cigar_length,uint64_t* const alignment_score,uint64_t* const effective_length,
    mm_stack_t* const mm_stack);

#endif /* SWG_ALIGN_H_ */
