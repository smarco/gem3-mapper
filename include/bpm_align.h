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
#include "match_elements.h"

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
typedef struct {
  uint64_t* Pv;
  uint64_t* Mv;
  uint64_t min_score;
  uint64_t min_score_column;
} bpm_align_matrix_t;

/*
 * Compile Pattern
 */
void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,uint8_t* const pattern,
    const uint64_t pattern_length,const uint64_t max_error,mm_stack_t* const mm_stack);

/*
 * BPM-distance (BitParalellMyers, Bit-compressed Alignment)
 *   Myers' Fast Bit-Vector algorithm to compute levenshtein distance
 */
// Raw
bool bpm_get_distance_raw(
    bpm_pattern_t* const bpm_pattern,const uint8_t* const text,const uint64_t text_length,
    uint64_t* const position,uint64_t* const distance);
// Cut-off
bool bpm_get_distance_cutoff(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,const uint64_t text_length,
    uint64_t* const match_end_column,uint64_t* const distance,
    const uint64_t max_distance,const bool quick_abandon);
// BPM Tiled (bound)
void bpm_get_distance_cutoff_tiled(
    bpm_pattern_t* const bpm_pattern,const uint8_t* const text,const uint64_t text_length,
    uint64_t* const levenshtein_distance,uint64_t* const levenshtein_match_end_column,
    const uint64_t max_error);
// Find all local minimums
uint64_t bpm_search_all(
    const bpm_pattern_t* const bpm_pattern,vector_t* const filtering_regions,
    const uint64_t text_trace_offset,const uint64_t index_position,
    const uint8_t* const text,const uint64_t text_length,const uint64_t max_distance);

/*
 * BPM-alignment (BitParalellMyers, Bit-compressed Alignment)
 *   Myers' Fast Bit-Vector algorithm to compute levenshtein alignment (CIGAR)
 */
#include "match_align_dto.h"

/*
 * BPM. Compute BPM-DP-Matrix
 *   @align_input->key
 *   @align_input->bpm_pattern
 *   @align_input->text
 *   @align_input->text_length
 */
void bpm_align_compute_matrix(
    match_align_input_t* const align_input,const uint64_t max_distance,
    bpm_align_matrix_t* const bpm_align_matrix,mm_stack_t* const mm_stack);
/*
 * BPM. Recover CIGAR from a matching string
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text
 *   @bpm_align_matrix->Pv
 *   @bpm_align_matrix->Mv
 *   @bpm_align_matrix->min_score
 *   @match_alignment->match_position (Adjusted)
 */
void bpm_align_backtrack_matrix(
    match_align_input_t* const align_input,const bool left_gap_alignment,
    bpm_align_matrix_t* const bpm_align_matrix,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector);
/*
 * BPM Align match
 *   @align_input->key
 *   @align_input->bpm_pattern
 *   @align_input->text
 *   @align_input->text_length
 *   @match_alignment->match_position (Adjusted)
 */
void bpm_align_match(
    match_align_input_t* const align_input,const uint64_t max_distance,
    const bool left_gap_alignment,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,mm_stack_t* const mm_stack);

#endif /* BPM_ALIGN_H_ */
