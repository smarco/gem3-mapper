/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "text/pattern.h"
#include "archive/archive_text_rl.h"
#include "gpu/gpu_buffer_align_bpm.h"
#include "gpu/gpu_config.h"

/*
 * Debug
 */
#define DEBUG_PATTERN_TILE_POSITION true

/*
 * Pattern Prepare
 */
void pattern_init_encode(
    pattern_t* const pattern,
    sequence_t* const sequence,
    const search_parameters_t* const parameters) {
  // Parameters
  const uint64_t sequence_length = sequence_get_length(sequence);
  const char* const read = sequence_get_read(sequence);
  // Aux
  uint64_t num_wildcards = 0;           // Number of bases not-allowed
  uint64_t num_low_quality_bases = 0;   // Number of bases with low quality value
  uint64_t num_non_canonical_bases = 0; // Number of bases not compliant with the k-mer filter
  int64_t clip_left, clip_right;
  int64_t i, j;
  // Set left-clipping (if any)
  switch (parameters->clipping) {
    case clipping_disabled:
      i = 0;
      clip_left = 0;
      break;
    case clipping_uncalled:
      i = 0;
      while (i<sequence_length && !is_dna_canonical(read[i])) ++i;
      clip_left = i;
      break;
    case clipping_masked:
      i = 0;
      while (i<sequence_length && !is_unmasked_dna(read[i])) ++i;
      clip_left = i;
      break;
    case clipping_fixed:
      i = parameters->clip_left;
      clip_left = MIN(parameters->clip_left,sequence_length);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Encode read
  if (pattern->quality_mask == NULL) {
    for (j=0;i<sequence_length;++i,++j) {
      const char character = read[i];
      if (!is_dna_canonical(character)) ++num_non_canonical_bases;
      if (!parameters->allowed_chars[(uint8_t)character]) ++num_wildcards;
      pattern->key[j] = dna_encode(character);
    }
  } else {
    // Encode read
    for (j=0;i<sequence_length;++i,++j) {
      const char character = read[i];
      if (!is_dna_canonical(character)) ++num_non_canonical_bases;
      if (!parameters->allowed_chars[(uint8_t)character]) {
        ++num_wildcards;
      } else if (pattern->quality_mask[i]!=qm_real) {
        ++num_low_quality_bases;
      }
      pattern->key[j] = dna_encode(character);
    }
  }
  // Skip right-clipping (if any)
  switch (parameters->clipping) {
    case clipping_disabled:
      clip_right = 0;
      break;
    case clipping_uncalled:
      i = sequence_length-1;
      while (i>=clip_left && !is_dna_canonical(read[i])) --i;
      clip_right = sequence_length - (i+1);
      break;
    case clipping_masked:
      i = sequence_length-1;
      while (i>=clip_left && !is_unmasked_dna(read[i])) --i;
      clip_right = sequence_length - (i+1);
      break;
    case clipping_fixed:
      clip_right = MIN(parameters->clip_right,sequence_length - clip_left);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Set effective length
  pattern->key_length = sequence_length - clip_right - clip_left;
  pattern->clip_left = clip_left;
  pattern->clip_right = clip_right;
  // Set wildcards & low-quality bases
  pattern->num_wildcards = num_wildcards;
  pattern->num_low_quality_bases = num_low_quality_bases;
  pattern->num_non_canonical_bases = num_non_canonical_bases;
}
void pattern_init_encode_rl(
    pattern_t* const pattern,
    const bool run_length_pattern,
    mm_stack_t* const mm_stack) {
  if (run_length_pattern) {
    pattern->run_length = true;
    // Allocate
    pattern->rl_key = mm_stack_calloc(mm_stack,pattern->key_length,uint8_t,false);
    pattern->rl_runs_acc = mm_stack_calloc(mm_stack,pattern->key_length,uint32_t,false);
    // RL encode
    archive_text_rl_encode(pattern->key,
        pattern->key_length,pattern->rl_key,
        &pattern->rl_key_length,pattern->rl_runs_acc);
    //    // DEBUG
    //    fprintf(stderr,">key\n");
    //    dna_buffer_print(stderr,pattern->regular_key,pattern->regular_key_length,false);
    //    fprintf(stderr,"\n");
    //    fprintf(stderr,">RL.key\n");
    //    dna_buffer_print(stderr,pattern->rl_key,pattern->rl_key_length,false);
    //    fprintf(stderr,"\n");
  } else {
    pattern->run_length = false;
    pattern->rl_key = NULL;
    pattern->rl_runs_acc = NULL;
  }
}
void pattern_init(
    pattern_t* const pattern,
    sequence_t* const sequence,
    bool* const do_quality_search,
    const search_parameters_t* const parameters,
    const bool run_length_pattern,
    const bool kmer_filter_compile,
    mm_stack_t* const mm_stack) {
  // Allocate pattern memory
  const uint64_t sequence_length = sequence_get_length(sequence);
  pattern->key = mm_stack_calloc(mm_stack,sequence_length,uint8_t,false);
  // Set quality search & Build quality model & mask
  *do_quality_search = (parameters->qualities_format!=sequence_qualities_ignore) && sequence_has_qualities(sequence);
  if (*do_quality_search) {
    pattern->quality_mask = mm_stack_calloc(mm_stack,sequence_length,uint8_t,false);
    sequence_qualities_model_process(sequence,
        parameters->qualities_model,parameters->qualities_format,
        parameters->quality_threshold,pattern->quality_mask);
  } else {
    pattern->quality_mask = NULL;
  }
  /*
   * Compute the encoded pattern
   *   Check all characters in the key & encode key
   *   Counts the number of wildcards (characters not allowed as replacements) and low_quality_bases
   */
  // Compute then pattern
  pattern_init_encode(pattern,sequence,parameters);
  // Compute the RL-pattern
  pattern_init_encode_rl(pattern,run_length_pattern,mm_stack);
  /*
   * Compute the effective number of differences
   * and compile the BMP-Pattern & k-mer filter
   */
  if (pattern->key_length == 0) { // Check length
    pattern->max_effective_filtering_error = 0;
    pattern->max_effective_bandwidth = 0;
    pattern->bpm_pattern = NULL;
    pattern->bpm_pattern_tiles = NULL;
  } else {
    const uint64_t max_effective_filtering_error =
        parameters->alignment_max_error_nominal +
        pattern->num_low_quality_bases +
        pattern->num_non_canonical_bases +
        pattern->num_wildcards;
    const uint64_t max_effective_bandwidth =
        parameters->alignment_max_bandwidth_nominal +
        pattern->num_low_quality_bases +
        pattern->num_non_canonical_bases +
        pattern->num_wildcards;
    pattern->max_effective_filtering_error = MIN(max_effective_filtering_error,pattern->key_length);
    pattern->max_effective_bandwidth = MIN(max_effective_bandwidth,pattern->key_length);
    // if (pattern->max_effective_filtering_error > 0) return; // Old-style
    // Prepare kmer-counting filter
    if (kmer_filter_compile) {
      kmer_counting_compile(&pattern->kmer_counting,pattern->key,
          pattern->key_length,pattern->max_effective_filtering_error,mm_stack);
    } else {
      pattern->kmer_counting.enabled = false;
    }
    // Prepare BPM pattern
    pattern->bpm_pattern = bpm_pattern_compile(pattern->key,pattern->key_length,mm_stack);
    pattern->bpm_pattern_tiles = bpm_pattern_compile_tiles(pattern->bpm_pattern,
        PATTERN_BPM_WORDS64_PER_TILE,pattern->max_effective_filtering_error,mm_stack);
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
void pattern_tiled_init(
    pattern_tiled_t* const pattern_tiled,
    const uint64_t pattern_length,
    const uint64_t pattern_tile_length,
    const uint64_t sequence_length,
    const uint64_t max_error) {
  // Calculate default tile dimensions & position
  pattern_tiled->pattern_max_error = max_error;
  pattern_tiled->pattern_tile_tall = pattern_tile_length;
  pattern_tiled->pattern_remaining_length = pattern_length;
  const int64_t pattern_band_width = sequence_length - pattern_length + 2*max_error;
  if (pattern_band_width < 0) {
    // It's supposed to be an invalid case because the filtering-process should have
    // detected (sequence_length < pattern_length) and trim the pattern accordingly
    gem_fatal_error_msg("Pattern.Tiled. Wrong Dimensions");
  }
  pattern_tiled->pattern_band_width = pattern_band_width;
  pattern_tiled->sequence_length = sequence_length;
  pattern_tiled->tile_offset = 0;
  // Calculate current tile dimensions (adjusted to initial conditions)
  pattern_tiled->tile_next_offset_inc = BOUNDED_SUBTRACTION(pattern_tile_length,max_error,0);
  pattern_tiled->tile_tall = MIN(pattern_tile_length,pattern_length);
  pattern_tiled->tile_wide = pattern_tiled->pattern_band_width + pattern_tiled->tile_next_offset_inc;
  if (pattern_tiled->tile_offset+pattern_tiled->tile_wide > pattern_tiled->sequence_length) {
    pattern_tiled->tile_wide = pattern_tiled->sequence_length - pattern_tiled->tile_offset;
  }
  pattern_tiled->prev_tile_match_position = UINT64_MAX; // Init last tile-matching column
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
  pattern_tiled->tile_tall = MIN(pattern_tiled->tile_tall,pattern_tiled->pattern_remaining_length);
  pattern_tiled->tile_wide = pattern_tiled->pattern_band_width + pattern_tiled->tile_next_offset_inc;
  if (pattern_tiled->tile_offset+pattern_tiled->tile_wide > pattern_tiled->sequence_length) {
    pattern_tiled->tile_wide = pattern_tiled->sequence_length - pattern_tiled->tile_offset;
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
    pattern_t* const pattern,
    bpm_pattern_t** const bpm_pattern_trimmed,
    bpm_pattern_t** const bpm_pattern_trimmed_tiles,
    const uint64_t key_trimmed_length,
    const uint64_t key_trim_left,
    mm_stack_t* const mm_stack) {
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Compile BPM-Pattern Trimmed
  *bpm_pattern_trimmed = bpm_pattern_compile(pattern->key+key_trim_left,key_trimmed_length,mm_stack);
  *bpm_pattern_trimmed_tiles = bpm_pattern_compile_tiles(
      *bpm_pattern_trimmed,PATTERN_BPM_WORDS64_PER_TILE,max_error,mm_stack);
}
/*
 * Display
 */
void pattern_enc_print(
    FILE* const stream,
    const uint8_t* const key,
    const uint64_t key_length) {
  uint64_t i;
  for (i=0;i<key_length;++i) {
    fprintf(stream,"%c",dna_decode(key[i]));
  }
}
