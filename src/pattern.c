/*
 * PROJECT: GEMMapper
 * FILE: pattern.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "pattern.h"
#include "sampled_rl.h"

/*
 * Debug
 */
#define DEBUG_PATTERN_TILE_POSITION true

/*
 * Pattern Prepare
 */
GEM_INLINE void pattern_prepare(
    sequence_t* const sequence,pattern_t* const pattern,region_profile_t* const region_profile,
    const as_parameters_t* const actual_parameters,const bool prepare_rl_pattern,
    bool* const do_quality_search,mm_stack_t* const mm_stack) {
  // Parameters
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  // Set quality search
  *do_quality_search = (parameters->quality_format!=qualities_ignore) && sequence_has_qualities(sequence);
  // Allocate pattern memory
  const uint64_t read_length = sequence_get_length(sequence);
  pattern->key_length = read_length;
  pattern->key = mm_stack_calloc(mm_stack,read_length,uint8_t,false);
  if (prepare_rl_pattern) {
    pattern->rl_key = mm_stack_calloc(mm_stack,read_length,uint8_t,false);
  }
  // Build quality model & mask
  if (*do_quality_search) {
    pattern->quality_mask =  mm_stack_calloc(mm_stack,read_length,uint8_t,false);
    quality_model(sequence,parameters->quality_model,
        parameters->quality_format,parameters->quality_threshold,pattern->quality_mask);
  } else {
    pattern->quality_mask = NULL;
  }
  /*
   * Check all characters in the key & encode key
   * Counts the number of wildcards(characters not allowed as replacements) & low_quality_bases
   */
  uint64_t i, num_wildcards=0, num_low_quality_bases=0;
  const char* const read = sequence_get_read(sequence);
  if (pattern->quality_mask == NULL) {
    for (i=0;i<read_length;++i) {
      const char character = read[i];
      if (!parameters->allowed_chars[(uint8_t)character]) ++num_wildcards;
      pattern->key[i] = dna_encode(character);
    }
  } else {
    for (i=0;i<read_length;++i) {
      const char character = read[i];
      if (!parameters->allowed_chars[(uint8_t)character]) {
        ++num_low_quality_bases; ++num_wildcards;
      } else if (pattern->quality_mask[i]!=qm_real) {
        ++num_low_quality_bases;
      }
      pattern->key[i] = dna_encode(character);
    }
  }
  pattern->num_wildcards = num_wildcards;
  pattern->num_low_quality_bases = num_low_quality_bases;
  // Compute the RL-pattern
  if (prepare_rl_pattern) {
    uint64_t rl_key_length = 1, run_length = 1;
    pattern->rl_key[0] = pattern->key[0];
    for (i=1;i<read_length;++i) {
      if (pattern->key[i] == pattern->key[i-1] && run_length < SAMPLED_RL_MAX_RUN_LENGTH) {
        ++run_length;
      } else {
        pattern->rl_key[rl_key_length++] = pattern->key[i];
        run_length = 1;
      }
    }
    pattern->rl_key_length = rl_key_length;
  }
  // Compute the effective number of differences
  const int64_t max_allowed_error = (int64_t)read_length - (int64_t)actual_parameters->alignment_min_identity_nominal;
  uint64_t effective_filtering_max_error;
  if (gem_expect_false(max_allowed_error<=0)) { // Constrained by min_matching_length_nominal
    effective_filtering_max_error = 0;
  } else {
    // Constrained by num_low_quality_bases
    effective_filtering_max_error = actual_parameters->alignment_max_error_nominal + pattern->num_low_quality_bases;
    if (effective_filtering_max_error > max_allowed_error) {
      effective_filtering_max_error = max_allowed_error;
    }
  }
  pattern->max_effective_filtering_error = effective_filtering_max_error;
  pattern->max_effective_bandwidth = actual_parameters->max_bandwidth_nominal + pattern->num_low_quality_bases;
  if (effective_filtering_max_error > 0) {
//    // Prepare kmer-counting filter
//    kmer_counting_compile(&pattern->kmer_counting,(bool* const)parameters->allowed_enc,pattern->key,read_length,mm_stack);
    // Prepare BPM pattern
    bpm_pattern_compile(&pattern->bpm_pattern,pattern->key,read_length,effective_filtering_max_error,mm_stack);
//    // Prepare SWG query-profile
//    if (parameters->alignment_model == alignment_model_gap_affine) {
//      swg_init_query_profile(&pattern->swg_query_profile,&parameters->swg_penalties,read_length,mm_stack);
//    }
  }
}
GEM_INLINE void pattern_clear(pattern_t* const pattern) {
  pattern->key_length = 0; // Clear the pattern
}
GEM_INLINE bool pattern_is_null(pattern_t* const pattern) {
  return (pattern->key_length == 0);
}
/*
 * Pattern Tiling
 */
GEM_INLINE bool pattern_tiled_init(
    pattern_tiled_t* const pattern_tiled,
    const uint64_t pattern_length,const uint64_t pattern_tile_tall,
    const uint64_t sequence_length,const uint64_t max_error) {
  // Calculate default tile dimensions & position
  pattern_tiled->pattern_max_error = max_error;
  pattern_tiled->pattern_tile_tall = pattern_tile_tall;
  pattern_tiled->pattern_remaining_length = pattern_length;
  const int64_t pattern_band_width = sequence_length - pattern_length + 2*max_error;
  if (pattern_band_width < 0) {
    return false; // Cannot align
  }
  pattern_tiled->pattern_band_width = pattern_band_width;
  pattern_tiled->sequence_length = sequence_length;
  pattern_tiled->tile_offset = 0;
  // Calculate current tile dimensions (adjusted to initial conditions)
  pattern_tiled->tile_next_offset_inc = pattern_tile_tall - max_error;
  pattern_tiled->tile_tall = (pattern_length > pattern_tile_tall) ? pattern_tile_tall : pattern_length;
  pattern_tiled->tile_wide = pattern_tiled->pattern_band_width + pattern_tiled->tile_next_offset_inc;
  if (pattern_tiled->tile_offset+pattern_tiled->tile_wide > sequence_length) {
    pattern_tiled->tile_wide = sequence_length - pattern_tiled->tile_offset;
  }
  // Init last tile-matching column
  pattern_tiled->prev_tile_match_position = UINT64_MAX;
  // Return Ok
  return true;
}
GEM_INLINE void pattern_tiled_calculate_next(pattern_tiled_t* const pattern_tiled) {
  // DEBUG
//  gem_cond_debug_block(DEBUG_PATTERN_TILE_POSITION) {
//    fprintf(stderr,">Tile (pos=%"PRIu64",len=%"PRIu64",tall=%"PRIu64") [distance=%"PRIu64",match_col=%"PRIu64"]\n",
//        pattern_tiled->tile_offset,pattern_tiled->tile_wide,pattern_tiled->tile_tall,
//        pattern_tiled->tile_distance,pattern_tiled->tile_match_column+pattern_tiled->tile_offset);
//  }
  // Update tile dimensions
  pattern_tiled->tile_offset += pattern_tiled->tile_next_offset_inc;
  pattern_tiled->pattern_remaining_length -= pattern_tiled->tile_tall;
  // Calculate current tile dimensions
  pattern_tiled->tile_next_offset_inc = pattern_tiled->tile_tall;
  pattern_tiled->tile_tall = (pattern_tiled->pattern_remaining_length > pattern_tiled->tile_tall) ?
      pattern_tiled->tile_tall : pattern_tiled->pattern_remaining_length;
  pattern_tiled->tile_wide = pattern_tiled->pattern_band_width + pattern_tiled->tile_next_offset_inc;
  if (pattern_tiled->tile_offset+pattern_tiled->tile_wide > pattern_tiled->sequence_length) {
    pattern_tiled->tile_wide = pattern_tiled->sequence_length-pattern_tiled->tile_offset;
  }
}
GEM_INLINE uint64_t pattern_tiled_bound_matching_path(pattern_tiled_t* const pattern_tiled) {
  if (pattern_tiled->prev_tile_match_position!=UINT64_MAX) {
    const int64_t prev_tile_match_position = pattern_tiled->prev_tile_match_position;
    const int64_t tile_match_position = pattern_tiled->tile_match_column + pattern_tiled->tile_offset;
    const int64_t tile_match_position_proyection = tile_match_position - pattern_tiled->tile_tall;
    // Calculate plausible match-band limits
    const int64_t tile_match_band_begin =
        BOUNDED_SUBTRACTION(tile_match_position_proyection,pattern_tiled->tile_distance,0);
    const int64_t tile_match_band_end =
        BOUNDED_ADDITION(tile_match_position_proyection,pattern_tiled->tile_distance,pattern_tiled->sequence_length-1);
    // Calculate differences
    const int64_t band_begin_difference = ABS(prev_tile_match_position - tile_match_band_begin);
    const int64_t band_end_difference = ABS(prev_tile_match_position - tile_match_band_end);
    // Keep matching column
    pattern_tiled->prev_tile_match_position = tile_match_position;
    // Return bound
    return MAX(band_begin_difference,band_end_difference);
  } else {
    // Keep matching column
    pattern_tiled->prev_tile_match_position = pattern_tiled->tile_match_column;
    return 0;
  }
}
/*
 * Display
 */
GEM_INLINE void pattern_enc_print(FILE* const stream,const uint8_t* const key,const uint64_t key_length) {
  uint64_t i;
  for (i=0;i<key_length;++i) {
    fprintf(stream,"%c",dna_decode(key[i]));
  }
}
