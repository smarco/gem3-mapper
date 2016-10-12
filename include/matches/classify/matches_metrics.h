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

#ifndef MATCHES_METRICS_H_
#define MATCHES_METRICS_H_

#include "utils/essentials.h"

/*
 * Constants
 */
#define MAX_TEMPLATE_LENGTH_SIGMAS 10.0

/*
 * Matches-Metrics (Current metrics status -always being updated-)
 */
typedef struct {
  /* Search Magnitudes */
  double proper_length;
  uint64_t read_length;
  int32_t swg_match_score;
  /* Read Mappability */
  uint64_t max_region_length;
  double kmer_frequency;
  /* Matches */
  uint64_t accepted_candidates;      // Total candidates accepted
  bool limited_candidates;           // Number of candidates has been limited due to selection
  uint64_t accepted_matches;         // Total matches accepted
  /* Matches Distance */
  uint64_t min_event_distance;       // Minimum event-distance recorded
  uint64_t min_event_distance_count; // Total matches with minimum event-distance
  uint64_t min_edit_distance;        // Minimum edit-distance recorded
  uint64_t min_edit_distance_count;  // Total matches with minimum edit-distance
  int32_t max_swg_score;             // Maximum smith-waterman-gotoh score recorded
  uint64_t max_swg_score_count;      // Total matches with maximum smith-waterman-gotoh score
  /* Matches template-length */
  double min_template_length_sigma;         // Minimum number of sigma-deviations (from the mean template-length)
  uint64_t min_template_length_sigma_count; // Total matches with minimum number of sigma-deviations
  /* Matches MAPQ */
  uint8_t mapq;                      // MAPQ score of the primary alignment
} matches_metrics_t;

/*
 * Setup
 */
void matches_metrics_init(matches_metrics_t* const metrics);

/*
 * Accessors
 */
uint64_t matches_metrics_get_min_event_distance(matches_metrics_t* const metrics);
uint64_t matches_metrics_get_min_edit_distance(matches_metrics_t* const metrics);
int32_t matches_metrics_get_max_swg_score(matches_metrics_t* const metrics);

void matches_metrics_set_proper_length(
    matches_metrics_t* const metrics,
    const double proper_length);
void matches_metrics_set_read_length(
    matches_metrics_t* const metrics,
    const uint64_t read_length);
void matches_metrics_set_swg_match_score(
    matches_metrics_t* const metrics,
    const int32_t swg_match_score);
void matches_metrics_set_max_region_length(
    matches_metrics_t* const metrics,
    const uint64_t max_region_length);
void matches_metrics_set_kmer_frequency(
    matches_metrics_t* const metrics,
    const double kmer_frequency);
void matches_metrics_add_accepted_candidates(
    matches_metrics_t* const metrics,
    const uint64_t accepted_candidates);
void matches_metrics_set_accepted_candidates(
    matches_metrics_t* const metrics,
    const uint64_t accepted_candidates);
void matches_metrics_set_limited_candidates(
    matches_metrics_t* const metrics,
    const bool limited);
void matches_metrics_set_mapq(
    matches_metrics_t* const metrics,
    const uint8_t mapq);

/*
 * Update
 */
void matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t event_distance,
    const uint64_t edit_distance,
    const int32_t swg_score);
void paired_matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t event_distance,
    const uint64_t edit_distance,
    const int32_t swg_score,
    const double template_length_sigma);

/*
 * Display
 */
void matches_metrics_print(
    FILE* const stream,
    matches_metrics_t* const matches_metrics);

#endif /* MATCHES_METRICS_H_ */
