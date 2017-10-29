/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_score.c
 * DATE: 03/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_map_score.h"

/*
 * Map Score Accessors
 */
GT_INLINE uint64_t gt_map_get_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->gt_score;
}
GT_INLINE void gt_map_set_score(gt_map* const map,const uint64_t score) {
  GT_MAP_CHECK(map);
  map->gt_score = score;
}
GT_INLINE uint8_t gt_map_get_phred_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->phred_score;
}
GT_INLINE void gt_map_set_phred_score(gt_map* const map,const uint8_t phred_score) {
  GT_MAP_CHECK(map);
  map->phred_score = phred_score;
}
