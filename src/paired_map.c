/*
 * PROJECT: GEMMapper
 * FILE: paired_map.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "paired_map.h"

/*
 * Accessors
 */
GEM_INLINE uint64_t paired_map_compute_distance(
    const match_trace_t* const match_end1,const match_trace_t* const match_end2) {
  return match_end1->distance+match_end2->distance;
}
GEM_INLINE uint64_t paired_map_compute_edit_distance(
    const match_trace_t* const match_end1,const match_trace_t* const match_end2) {
  return match_end1->edit_distance+match_end2->edit_distance;
}
GEM_INLINE uint64_t paired_map_compute_swg_score(
    const match_trace_t* const match_end1,const match_trace_t* const match_end2) {
  return match_end1->swg_score+match_end2->swg_score;
}
GEM_INLINE uint64_t paired_map_get_distance(paired_map_t* const paired_map) {
  return paired_map->distance;
}
GEM_INLINE uint64_t paired_map_get_edit_distance(paired_map_t* const paired_map) {
  return paired_map->edit_distance;
}
GEM_INLINE int32_t paired_map_get_swg_score(paired_map_t* const paired_map) {
  return paired_map->swg_score;
}
/*
 * Utils
 */
GEM_INLINE uint64_t paired_map_get_template_observed_length(
    const uint64_t begin_position_1,const uint64_t end_position_1,
    const uint64_t begin_position_2,const uint64_t end_position_2) {
  const uint64_t most_right = MAX(end_position_1,end_position_2);
  const uint64_t most_left = MIN(begin_position_1,begin_position_2);
  return most_right - most_left;
}
