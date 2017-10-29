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

#ifndef PAIRED_MAP_H_
#define PAIRED_MAP_H_

#include "utils/essentials.h"
#include "archive/search/archive_search_pe_parameters.h"
#include "matches/matches.h"

/*
 * Paired Matches
 */
typedef struct {
  match_trace_t* match_trace_end1;     // Map end1
  match_trace_t* match_trace_end2;     // Map end2
  uint8_t mapq_score;                  // MAPQ score
  uint64_t template_length;            // Template observed length
  double template_length_sigma;        // Number of sigmas deviation from the std-distribution
  uint64_t event_distance;             // Distance of the paired-alignment
  uint64_t edit_distance;              // Edit Distance of the paired-alignment
  int32_t swg_score;                   // Distance of the paired-alignment
  float error_quality;                 // Average errors quality
  pair_relation_t pair_relation;       // Pair relation (concordant/discordant)
  pair_orientation_t pair_orientation; // Pair orientation (FR,RF,FF,RR)
  pair_layout_t pair_layout;           // Pair layout (pair_layout_separate,pair_layout_overlap,pair_layout_contain)
  uint64_t index_position;             // Pair index position (only used for sorting purposes)
} paired_map_t;

/*
 * Accessors
 */
uint64_t paired_map_compute_event_distance(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2);
uint64_t paired_map_compute_edit_distance(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2);
int32_t paired_map_compute_swg_score(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2);
float paired_map_compute_error_quality(
    const match_trace_t* const match_end1,
    const match_trace_t* const match_end2);

uint64_t paired_map_get_event_distance(paired_map_t* const paired_map);
uint64_t paired_map_get_edit_distance(paired_map_t* const paired_map);
int32_t paired_map_get_swg_score(paired_map_t* const paired_map);

#endif /* PAIRED_MAP_H_ */
