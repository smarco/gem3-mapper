/*
 * PROJECT: GEMMapper
 * FILE: paired_map.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/paired_map.h"

/*
 * Accessors
 */
uint64_t paired_map_compute_distance(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2) {
  return match_end1->distance+match_end2->distance;
}
uint64_t paired_map_compute_edit_distance(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2) {
  return match_end1->edit_distance+match_end2->edit_distance;
}
uint64_t paired_map_compute_swg_score(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2) {
  return match_end1->swg_score+match_end2->swg_score;
}
uint64_t paired_map_get_distance(paired_map_t* const paired_map) {
  return paired_map->distance;
}
uint64_t paired_map_get_edit_distance(paired_map_t* const paired_map) {
  return paired_map->edit_distance;
}
int32_t paired_map_get_swg_score(paired_map_t* const paired_map) {
  return paired_map->swg_score;
}
