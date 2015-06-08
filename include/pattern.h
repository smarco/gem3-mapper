/*
 * PROJECT: GEMMapper
 * FILE: pattern.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PATTERN_H_
#define PATTERN_H_

#include "essentials.h"
#include "approximate_search_parameters.h"
#include "region_profile.h"
#include "bpm_align.h"
#include "kmer_counting.h"
#include "swg_align.h"

/*
 * Approximate Search Pattern
 */
typedef struct {
  /* Processed Search Pattern */
  uint8_t* key;            // Encoded Pattern
  uint8_t* quality_mask;   // Quality Mask
  uint64_t key_length;     // Total Length
  /* RunLength Processed Search Pattern */
  uint8_t* rl_key;         // Encoded RL-Pattern
  uint64_t rl_key_length;  // Total RL-Length
  /* Pattern Properties */
  uint64_t num_wildcards;
  uint64_t num_low_quality_bases;
  uint64_t max_effective_filtering_error;
  uint64_t max_effective_bandwidth;
  /* K-mers counting */
  kmer_counting_t kmer_counting;
  /* Pattern BitVector-Encoded (Myers-DP) */
  bpm_pattern_t bpm_pattern;
  /* SWG query profile */
  swg_query_profile_t swg_query_profile; // FIXME Remove me
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
 * Pattern Prepare
 */
GEM_INLINE void pattern_prepare(
    sequence_t* const sequence,pattern_t* const pattern,region_profile_t* const region_profile,
    const as_parameters_t* const actual_parameters,const bool prepare_rl_pattern,
    bool* const do_quality_search,mm_stack_t* const mm_stack);
GEM_INLINE void pattern_clear(pattern_t* const pattern);
GEM_INLINE bool pattern_is_null(pattern_t* const pattern);

/*
 * Pattern Tiling
 */
GEM_INLINE bool pattern_tiled_init(
    pattern_tiled_t* const pattern_tiled,
    const uint64_t pattern_length,const uint64_t pattern_tile_tall,
    const uint64_t sequence_length,const uint64_t max_error);
GEM_INLINE void pattern_tiled_calculate_next(pattern_tiled_t* const pattern_tiled);
GEM_INLINE uint64_t pattern_tiled_bound_matching_path(pattern_tiled_t* const pattern_tiled);

/*
 * Display
 */
GEM_INLINE void pattern_enc_print(FILE* const stream,const uint8_t* const key,const uint64_t key_length);

#endif /* PATTERN_H_ */
