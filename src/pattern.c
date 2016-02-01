/*
 * PROJECT: GEMMapper
 * FILE: pattern.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "pattern.h"
#include "sampled_rl.h"
#include "gpu_buffer_align_bpm.h"
#include "gpu_config.h"

/*
 * Debug
 */
#define DEBUG_PATTERN_TILE_POSITION true

/*
 * Pattern Prepare
 */
void pattern_init(
    pattern_t* const pattern,sequence_t* const sequence,
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
  uint64_t i, num_wildcards=0, num_low_quality_bases=0, num_non_canonical_bases = 0;
  const char* const read = sequence_get_read(sequence);
  if (pattern->quality_mask == NULL) {
    for (i=0;i<read_length;++i) {
      const char character = read[i];
      if (!parameters->allowed_chars[(uint8_t)character]) ++num_wildcards;
      if (!is_dna_canonical(character)) ++num_non_canonical_bases;
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
      if (!is_dna_canonical(character)) ++num_non_canonical_bases;
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
  // Constrained by num_low_quality_bases
  const uint64_t max_effective_filtering_error =
      actual_parameters->alignment_max_error_nominal + pattern->num_low_quality_bases;
  pattern->max_effective_filtering_error = max_effective_filtering_error;
  pattern->max_effective_bandwidth = actual_parameters->alignment_max_bandwidth_nominal + pattern->num_low_quality_bases;
  if (max_effective_filtering_error > 0) {
    // Prepare kmer-counting filter
    kmer_counting_compile(&pattern->kmer_counting,pattern->key,read_length,
        num_non_canonical_bases,max_effective_filtering_error,mm_stack);
    // Prepare BPM pattern
    pattern->bpm_pattern = bpm_pattern_compile(
        pattern->key,read_length,max_effective_filtering_error,mm_stack);
    pattern->bpm_pattern_tiles = bpm_pattern_compile_tiles(pattern->bpm_pattern,
        PATTERN_BPM_WORDS64_PER_TILE,max_effective_filtering_error,mm_stack);
  }
}
void pattern_clear(pattern_t* const pattern) {
  pattern->key_length = 0; // Clear the pattern
}
bool pattern_is_null(pattern_t* const pattern) {
  return (pattern->key_length == 0);
}
/*
 * Pattern Tiling
 */
bool pattern_tiled_init(
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
void pattern_tiled_calculate_next(pattern_tiled_t* const pattern_tiled) {
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
uint64_t pattern_tiled_bound_matching_path(pattern_tiled_t* const pattern_tiled) {
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
 * Pattern Trimmed
 */
void pattern_trimmed_init(
    pattern_t* const pattern,bpm_pattern_t** const bpm_pattern_trimmed,
    bpm_pattern_t** const bpm_pattern_trimmed_tiles,const uint64_t key_trim_left,
    const uint64_t key_trim_right,mm_stack_t* const mm_stack) {
  const uint64_t key_length_trimmed = pattern->key_length - key_trim_left - key_trim_right;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Compile BPM-Pattern Trimmed
  *bpm_pattern_trimmed = bpm_pattern_compile(pattern->key+key_trim_left,
      key_length_trimmed,MIN(max_error,key_length_trimmed),mm_stack);
  *bpm_pattern_trimmed_tiles = bpm_pattern_compile_tiles(
      *bpm_pattern_trimmed,PATTERN_BPM_WORDS64_PER_TILE,max_error,mm_stack);
}
/*
 * Display
 */
void pattern_enc_print(FILE* const stream,const uint8_t* const key,const uint64_t key_length) {
  uint64_t i;
  for (i=0;i<key_length;++i) {
    fprintf(stream,"%c",dna_decode(key[i]));
  }
}
