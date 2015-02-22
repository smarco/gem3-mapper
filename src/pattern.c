/*
 * PROJECT: GEMMapper
 * FILE: pattern.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "pattern.h"

/*
 * Debug
 */
#define DEBUG_PATTERN_TILE_POSITION true

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
//    fprintf(stderr,">Tile (pos=%lu,len=%lu,tall=%lu) [distance=%lu,match_col=%lu]\n",
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
