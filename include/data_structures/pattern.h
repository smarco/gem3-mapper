/*
 * PROJECT: GEMMapper
 * FILE: pattern.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PATTERN_H_
#define PATTERN_H_

#include "utils/essentials.h"
#include "align/align_bpm_pattern.h"
#include "archive/archive_search_parameters.h"
#include "filtering/kmer_counting.h"

/*
 * Approximate Search Pattern
 */
typedef struct {
  /* Processed Search Pattern */
  uint8_t* key;                  // Encoded Pattern
  uint8_t* quality_mask;         // Quality Mask
  uint64_t key_length;           // Total Length
  /* Run-Length Pattern */
  bool run_length;
  uint8_t* rl_key;               // RL-Encoded Text
  uint64_t rl_key_length;        // RL-Encoded Text length
  uint32_t* rl_runs_acc;         // Length of each run (accumulated)
  /* Pattern Properties */
  uint64_t num_wildcards;
  uint64_t num_low_quality_bases;
  uint64_t max_effective_filtering_error;
  uint64_t max_effective_bandwidth;
  /* K-mers counting */
  kmer_counting_t kmer_counting;
  /* Pattern BitVector-Encoded (Myers-DP) */
  bpm_pattern_t* bpm_pattern;
  bpm_pattern_t* bpm_pattern_tiles;
} pattern_t;

/*
 * Tiled Pattern
 */
typedef struct {
  // Current tile dimensions
  uint64_t tile_offset;
  uint64_t tile_wide;
  uint64_t tile_tall;
  uint64_t tile_next_offset_inc;
  // Complete tile dimensions (defaults)
  uint64_t pattern_band_width;
  uint64_t pattern_tile_offset;
  uint64_t pattern_tile_wide;
  uint64_t pattern_tile_tall;
  uint64_t pattern_max_error;
  uint64_t pattern_remaining_length;
  uint64_t sequence_length;
  // Tile matching information
  uint64_t tile_distance;
  uint64_t tile_match_column;
  uint64_t prev_tile_match_position;
} pattern_tiled_t;

/*
 * Constants
 */
#define PATTERN_BPM_WORDS64_PER_TILE 4

/*
 * Pattern Prepare
 */
void pattern_init(
    pattern_t* const pattern,
    sequence_t* const sequence,
    bool* const do_quality_search,
    const search_parameters_t* const parameters,
    const bool run_length_pattern,
    const bool kmer_filter_compile,
    mm_stack_t* const mm_stack);
void pattern_clear(pattern_t* const pattern);
bool pattern_is_null(pattern_t* const pattern);

/*
 * Pattern Tiling
 */
void pattern_tiled_init(
    pattern_tiled_t* const pattern_tiled,
    const uint64_t pattern_length,
    const uint64_t pattern_tile_length,
    const uint64_t sequence_length,
    const uint64_t max_error);
void pattern_tiled_calculate_next(pattern_tiled_t* const pattern_tiled);
uint64_t pattern_tiled_bound_matching_path(pattern_tiled_t* const pattern_tiled);

/*
 * Pattern Trimmed
 */
void pattern_trimmed_init(
    pattern_t* const pattern,
    bpm_pattern_t** const bpm_pattern_trimmed,
    bpm_pattern_t** const bpm_pattern_trimmed_tiles,
    const uint64_t key_trimmed_length,
    const uint64_t key_trim_left,
    const uint64_t key_trim_right,
    mm_stack_t* const mm_stack);

/*
 * Display
 */
void pattern_enc_print(
    FILE* const stream,
    const uint8_t* const key,
    const uint64_t key_length);

#endif /* PATTERN_H_ */
