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

#include "align/pattern/pattern.h"
#include "archive/archive_text_rl.h"
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
  // Set left-clipping (if any)
  uint64_t num_wildcards = 0;           // Number of bases not-allowed
  int64_t clip_left, clip_right;
  int64_t i, j;
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
  for (j=0;i<sequence_length;++i,++j) {
    const char character = read[i];
    if (!is_dna_canonical(character)) ++num_wildcards;
    pattern->key[j] = dna_encode(character);
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
  // Encode qualities
  if (pattern->quality_mask!=NULL) {
    char* const qualities = sequence_get_qualities(sequence) + pattern->clip_left;
    sequence_qualities_model_process(
        parameters->qualities_format,qualities,
        pattern->key_length,pattern->quality_mask);
  }
}
void pattern_init_encode_rl(
    pattern_t* const pattern,
    const bool run_length_pattern,
    mm_allocator_t* const mm_allocator) {
  if (run_length_pattern) {
    pattern->run_length = true;
    // Allocate
    pattern->rl_key = mm_allocator_calloc(mm_allocator,pattern->key_length,uint8_t,false);
    pattern->rl_runs_acc = mm_allocator_calloc(mm_allocator,pattern->key_length,uint32_t,false);
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
    mm_allocator_t* const mm_allocator) {
  // Allocate pattern memory
  const uint64_t sequence_length = sequence_get_length(sequence);
  pattern->key = mm_allocator_calloc(mm_allocator,sequence_length,uint8_t,false);
  // Set quality search & Build quality model & mask
  *do_quality_search =
      (parameters->qualities_format!=sequence_qualities_ignore) &&
       sequence_has_qualities(sequence);
  if (*do_quality_search) {
    pattern->quality_mask = mm_allocator_calloc(mm_allocator,sequence_length,uint8_t,true);
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
  pattern_init_encode_rl(pattern,run_length_pattern,mm_allocator);
  /*
   * Compute the effective number of differences
   * and compile the BMP-Pattern & k-mer filter
   */
  const uint64_t key_length = pattern->key_length;
  if (key_length == 0) { // Check length
    pattern->max_effective_filtering_error = 0;
    pattern->max_effective_bandwidth = 0;
  } else {
    // Set error
    const uint64_t max_effective_filtering_error =
        parameters->alignment_max_error_nominal +
        pattern->num_wildcards;
    pattern->max_effective_filtering_error = MIN(max_effective_filtering_error,key_length);
    // Set bandwidth and adjust (no increment over BPM_MAX_TILE_LENGTH)
    const uint64_t max_effective_bandwidth =
        parameters->alignment_max_bandwidth_nominal +
        pattern->num_wildcards;
    pattern->max_effective_bandwidth = MIN(max_effective_bandwidth,key_length);
    if (key_length >= BPM_MAX_TILE_LENGTH) {
      const double max_bandwidth_rate = (double)pattern->max_effective_bandwidth/(double)key_length;
      pattern->max_effective_bandwidth = (uint64_t)ceil((double)BPM_MAX_TILE_LENGTH*max_bandwidth_rate);
    }
    // Set extension error
    const uint64_t max_extension_error =
        parameters->alignment_max_extension_error_nominal +
        pattern->num_wildcards;
    pattern->max_extension_error = MIN(max_extension_error,key_length);
    // Compile BPM filter
    pattern_tiled_compile(
        &pattern->pattern_tiled,pattern->key,
        key_length,pattern->max_effective_filtering_error,
        parameters->candidate_verification.kmer_tiles,
        parameters->candidate_verification.kmer_length,
        parameters->gpu_stage_kmer_filter_enabled,mm_allocator);
  }
}
void pattern_destroy(
    pattern_t* const pattern,
    mm_allocator_t* const mm_allocator) {
  // Pattern
  if (pattern->key != NULL) {
    mm_allocator_free(mm_allocator,pattern->key);
  }
  if (pattern->quality_mask != NULL) {
    mm_allocator_free(mm_allocator,pattern->quality_mask);
  }
  if (pattern->rl_key != NULL) {
    mm_allocator_free(mm_allocator,pattern->rl_key);
  }
  if (pattern->rl_runs_acc != NULL) {
    mm_allocator_free(mm_allocator,pattern->rl_runs_acc);
  }
  pattern->key_length = 0;
  // Pattern Tiles
  pattern_tiled_destroy(&pattern->pattern_tiled,mm_allocator);
}
bool pattern_is_null(
    pattern_t* const pattern) {
  return (pattern->key_length == 0);
}
/*
 * Pattern tiling
 */
void pattern_tiling_init(
    pattern_tiling_t* const pattern_tiling,
    const uint64_t pattern_length,
    const uint64_t pattern_tile_length,
    const uint64_t sequence_length,
    const uint64_t max_error) {
  // Calculate default tile dimensions & position
  pattern_tiling->pattern_max_error = max_error;
  pattern_tiling->pattern_tile_tall = pattern_tile_length;
  pattern_tiling->pattern_remaining_length = pattern_length;
  const int64_t pattern_band_width = sequence_length - pattern_length + 2*max_error;
  if (pattern_band_width < 0) {
    // It's supposed to be an invalid case because the filtering-process should have
    // detected (sequence_length < pattern_length) and trim the pattern accordingly
    gem_fatal_error_msg("Pattern tiling went wrong: invalid dimensions");
  }
  pattern_tiling->pattern_band_width = pattern_band_width;
  pattern_tiling->sequence_length = sequence_length;
  pattern_tiling->tile_offset = 0;
  // Calculate current tile dimensions (adjusted to initial conditions)
  pattern_tiling->tile_next_offset_inc = BOUNDED_SUBTRACTION(pattern_tile_length,max_error,0);
  pattern_tiling->tile_tall = MIN(pattern_tile_length,pattern_length);
  pattern_tiling->tile_wide = pattern_tiling->pattern_band_width + pattern_tiling->tile_next_offset_inc;
  if (pattern_tiling->tile_offset+pattern_tiling->tile_wide > pattern_tiling->sequence_length) {
    pattern_tiling->tile_wide = pattern_tiling->sequence_length - pattern_tiling->tile_offset;
  }
}
void pattern_tiling_next(
    pattern_tiling_t* const pattern_tiling) {
  // DEBUG
  //  gem_cond_debug_block(DEBUG_PATTERN_TILE_POSITION) {
  //    fprintf(stderr,">Tile (pos=%"PRIu64",len=%"PRIu64",tall=%"PRIu64") [distance=%"PRIu64",match_col=%"PRIu64"]\n",
  //        pattern_tiling->tile_offset,pattern_tiling->tile_wide,pattern_tiling->tile_tall,
  //        pattern_tiling->tile_distance,pattern_tiling->tile_match_column+pattern_tiling->tile_offset);
  //  }
  // Update tile dimensions
  pattern_tiling->tile_offset += pattern_tiling->tile_next_offset_inc;
  pattern_tiling->pattern_remaining_length -= pattern_tiling->tile_tall;
  // Calculate current tile dimensions
  pattern_tiling->tile_next_offset_inc = pattern_tiling->tile_tall;
  pattern_tiling->tile_tall = MIN(pattern_tiling->tile_tall,pattern_tiling->pattern_remaining_length);
  pattern_tiling->tile_wide = pattern_tiling->pattern_band_width + pattern_tiling->tile_next_offset_inc;
  if (pattern_tiling->tile_offset+pattern_tiling->tile_wide > pattern_tiling->sequence_length) {
    pattern_tiling->tile_wide = pattern_tiling->sequence_length - pattern_tiling->tile_offset;
  }
}
//uint64_t pattern_tiling_bound_matching_path(pattern_tiling_t* const pattern_tiling) {
//  if (pattern_tiling->prev_tile_match_position!=UINT64_MAX) {
//    const int64_t prev_tile_match_position = pattern_tiling->prev_tile_match_position;
//    const int64_t tile_match_position = pattern_tiling->tile_match_column + pattern_tiling->tile_offset;
//    const int64_t tile_match_position_proyection = tile_match_position - pattern_tiling->tile_tall;
//    // Calculate plausible match-band limits
//    const int64_t tile_match_band_begin =
//        BOUNDED_SUBTRACTION(tile_match_position_proyection,pattern_tiling->tile_distance,0);
//    const int64_t tile_match_band_end =
//        BOUNDED_ADDITION(tile_match_position_proyection,pattern_tiling->tile_distance,pattern_tiling->sequence_length-1);
//    // Calculate differences
//    const int64_t band_begin_difference = ABS(prev_tile_match_position - tile_match_band_begin);
//    const int64_t band_end_difference = ABS(prev_tile_match_position - tile_match_band_end);
//    // Keep matching column
//    pattern_tiling->prev_tile_match_position = tile_match_position;
//    // Return bound
//    return MAX(band_begin_difference,band_end_difference);
//  } else {
//    // Keep matching column
//    pattern_tiling->prev_tile_match_position = pattern_tiling->tile_match_column;
//    return 0;
//  }
//}
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
