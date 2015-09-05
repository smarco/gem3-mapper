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
#include "match_elements.h"

/*
 * Constants
 */
#define SWG_SCORE_MIN (INT16_MIN)

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
  int32_t gap_open_score;
  int32_t gap_extension_score;
  int32_t generic_match_score;
  int32_t generic_mismatch_score;
  matching_score_t matching_score;
} swg_penalties_t;
typedef struct {
  // 8-bits cells
  uint8_t match_bias_uint8;  // Query Profile Bias
  uint8_t matrix_bias_uint8; // Full matrix Bias
  bool overflow_uint8;
  uint64_t segment_length_uint8;
  uint64_t key_effective_length_uint8;
  uint8_t* query_profile_uint8[DNA__N_RANGE];
  // 16-bits cells
  uint64_t segment_length_int16;
  uint64_t key_effective_length_int16;
  int16_t* query_profile_int16[DNA__N_RANGE];
} swg_query_profile_t;

/*
 * Check CIGAR string
 */
bool align_check_match(
    FILE* const stream,const uint8_t* const key,const uint64_t key_length,const uint8_t* const text,
    const uint64_t text_length,vector_t* const cigar_vector,uint64_t const cigar_offset,
    uint64_t const cigar_length,const bool verbose);

/*
 * Levenshtein Alignment
 */
int64_t align_levenshtein_get_distance(
    const char* const key,const uint64_t key_length,
    const char* const text,const uint64_t text_length,
    const bool ends_free,uint64_t* const position);

/*
 * Init SWG Query Profile
 */
void swg_init_query_profile(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint64_t max_expected_key_length,mm_stack_t* const mm_stack);
bool swg_compile_query_profile_uint8(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,const uint64_t key_length,mm_stack_t* const mm_stack);
bool swg_compile_query_profile_int16(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,const uint64_t key_length,mm_stack_t* const mm_stack);

/*
 * SWG Score
 */
int32_t swg_score_deletion(const swg_penalties_t* const swg_penalties,const int32_t length);
int32_t swg_score_insertion(const swg_penalties_t* const swg_penalties,const int32_t length);
int32_t swg_score_mismatch(const swg_penalties_t* const swg_penalties);
int32_t swg_score_match(const swg_penalties_t* const swg_penalties,const int32_t match_length);
int32_t swg_score_cigar_element(
    const swg_penalties_t* const swg_penalties,const cigar_element_t* const cigar_element);
int32_t swg_score_cigar(
    const swg_penalties_t* const swg_penalties,vector_t* const cigar_vector,
    const uint64_t cigar_offset,const uint64_t cigar_length);
int32_t swg_score_cigar__excluding_clipping(
    const swg_penalties_t* const swg_penalties,vector_t* const cigar_vector,
    const uint64_t cigar_offset,const uint64_t cigar_length);

/*
 * Smith-waterman-gotoh Alignment
 */
//typedef struct match_align_input_t match_align_input_t;
//typedef struct match_align_parameters_t match_align_parameters_t;
#include "match_align_dto.h"

/*
 * Smith-waterman-gotoh Base (ie. no-optimizations)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @match_alignment->match_position (Adjusted)
 */
void swg_align_match_base(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,vector_t* const cigar_vector,mm_stack_t* const mm_stack);
/*
 * Smith-Waterman-Gotoh - Main procedure (Dispatcher)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @match_alignment->match_position (Adjusted)
 *   @match_alignment->cigar_length (Cumulative)
 */
void swg_align_match(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    const bool begin_free,const bool end_free,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,mm_stack_t* const mm_stack);

#endif /* SWG_ALIGN_H_ */
