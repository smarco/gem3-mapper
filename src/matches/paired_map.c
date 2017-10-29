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

#include "matches/paired_map.h"

/*
 * Accessors
 */
uint64_t paired_map_compute_event_distance(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2) {
  return match_end1->event_distance+match_end2->event_distance;
}
uint64_t paired_map_compute_edit_distance(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2) {
  return match_end1->edit_distance+match_end2->edit_distance;
}
int32_t paired_map_compute_swg_score(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2) {
  return match_end1->swg_score+match_end2->swg_score;
}
float paired_map_compute_error_quality(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2) {
  return (match_end1->error_quality+match_end2->error_quality)/2.0;
}
uint64_t paired_map_get_event_distance(paired_map_t* const paired_map) {
  return paired_map->event_distance;
}
uint64_t paired_map_get_edit_distance(paired_map_t* const paired_map) {
  return paired_map->edit_distance;
}
int32_t paired_map_get_swg_score(paired_map_t* const paired_map) {
  return paired_map->swg_score;
}
