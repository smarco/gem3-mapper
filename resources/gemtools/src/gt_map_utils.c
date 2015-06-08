/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_utils.h
 * DATE: 03/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_map_utils.h"

GT_INLINE bool gt_map_block_overlap(
    gt_map* const map_block_1,gt_map* const map_block_2,
    uint64_t* const overlap_start,uint64_t* const overlap_end) {
  // Global coordinates
  uint64_t map_1_start = gt_map_get_global_coordinate(map_block_1);
  uint64_t map_2_start = gt_map_get_global_coordinate(map_block_2);
  uint64_t map_1_end = map_1_start + gt_map_get_global_length(map_block_1);
  uint64_t map_2_end = map_2_start + gt_map_get_global_length(map_block_2);
  // Check overlap
  if (map_1_start <= map_2_start) {
    if (map_1_end < map_2_start) return false;
  } else if (map_2_start <= map_1_start) {
    if (map_2_end < map_1_start) return false;
  } else { // map_1_start == map_2_start
    return false;
  }
  // They do overlap
  if (map_1_start <= map_2_start) {
    if (overlap_start) *overlap_start = map_1_start + (map_2_start - map_1_start);
    if (overlap_end) *overlap_end = map_1_end;
  } else {
    if (overlap_start) *overlap_start = map_2_start + (map_1_start - map_2_start);
    if (overlap_end) *overlap_end = map_2_end;
  }
  return true;
}

