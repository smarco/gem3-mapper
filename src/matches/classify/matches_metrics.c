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

#include "matches/classify/matches_metrics.h"

/*
 * Setup
 */
void matches_metrics_init(matches_metrics_t* const metrics) {
  // Matches
  metrics->candidates_accepted = 0;
  // Matches Distance
  metrics->min1_event_distance = UINT32_MAX;
  metrics->min2_event_distance = UINT32_MAX;
  metrics->min1_edit_distance = UINT32_MAX;
  metrics->min2_edit_distance = UINT32_MAX;
  metrics->max1_swg_score = INT32_MIN;
  metrics->max2_swg_score = INT32_MIN;
  // MAPQ Score
  metrics->mapq = 0;
}
/* */
void matches_metrics_set_proper_length(
    matches_metrics_t* const metrics,
    const double proper_length) {
  metrics->proper_length = proper_length;
}
void matches_metrics_set_read_length(
    matches_metrics_t* const metrics,
    const uint64_t read_length) {
  metrics->read_length = read_length;
}
void matches_metrics_set_swg_match_score(
    matches_metrics_t* const metrics,
    const int32_t swg_match_score) {
  metrics->swg_match_score = swg_match_score;
}
void matches_metrics_set_region_profile_metrics(
    matches_metrics_t* const metrics,
    const uint64_t avg_region_length,
    const uint64_t max_region_length,
    const double kmer_frequency) {
  metrics->avg_region_length = avg_region_length;
  metrics->max_region_length = max_region_length;
  metrics->kmer_frequency = kmer_frequency;
}
void matches_metrics_add_accepted_candidates(
    matches_metrics_t* const metrics,
    const uint64_t candidates_accepted) {
  metrics->candidates_accepted += candidates_accepted;
}
void matches_metrics_set_mapq(
    matches_metrics_t* const metrics,
    const uint8_t mapq) {
  metrics->mapq = mapq;
}
/*
 * Update
 */
void matches_metrics_update(
    matches_metrics_t* const matches_metrics,
    const uint64_t event_distance,
    const uint64_t edit_distance,
    const int32_t swg_score) {
  // Min counter
  if (event_distance <= matches_metrics->min1_event_distance) {
    matches_metrics->min2_event_distance = matches_metrics->min1_event_distance;
    matches_metrics->min1_event_distance = event_distance;
  } else if (event_distance < matches_metrics->min2_event_distance) {
    matches_metrics->min2_event_distance = event_distance;
  }
  // Min Edit-distance
  if (edit_distance <= matches_metrics->min1_edit_distance) {
    matches_metrics->min2_edit_distance = matches_metrics->min1_edit_distance;
    matches_metrics->min1_edit_distance = edit_distance;
  } else if (edit_distance < matches_metrics->min2_edit_distance) {
    matches_metrics->min2_edit_distance = edit_distance;
  }
  // Max SWG-score
  if (swg_score >= matches_metrics->max1_swg_score) {
    matches_metrics->max2_swg_score = matches_metrics->max1_swg_score;
    matches_metrics->max1_swg_score = swg_score;
  } else if (swg_score > matches_metrics->max2_swg_score) {
    matches_metrics->max2_swg_score = swg_score;
  }
}
/*
 * Display
 */
void matches_metrics_print(
    FILE* const stream,
    matches_metrics_t* const matches_metrics) {
  tab_fprintf(stream,"[GEM]>Metrics\n");
  tab_fprintf(stream,"  => Search.Magnitudes\n");
  tab_fprintf(stream,"    => Read.length        %lu\n",matches_metrics->read_length);
  tab_fprintf(stream,"    => Proper.Length      %2.3f\n",matches_metrics->proper_length);
  tab_fprintf(stream,"    => SWG.Match.Score    %lu\n",matches_metrics->swg_match_score);
  tab_fprintf(stream,"  => Mappability\n");
  tab_fprintf(stream,"    => Max.Region.length  %lu\n",matches_metrics->max_region_length);
  tab_fprintf(stream,"    => Kmer.frequency     %2.3f\n",matches_metrics->kmer_frequency);
  tab_fprintf(stream,"  => Matches\n");
  tab_fprintf(stream,"    => Candidates.accepted       %lu\n",matches_metrics->candidates_accepted);
  tab_fprintf(stream,"    => Min1.counter.value        %ld\n",matches_metrics->min1_event_distance);
  tab_fprintf(stream,"    => Min2.counter.value        %ld\n",matches_metrics->min2_event_distance);
  tab_fprintf(stream,"    => Min1.edit.distance        %ld\n",matches_metrics->min1_edit_distance);
  tab_fprintf(stream,"    => Min2.edit.distance        %ld\n",matches_metrics->min2_edit_distance);
  tab_fprintf(stream,"    => Max1.swg.score            %ld\n",matches_metrics->max1_swg_score);
  tab_fprintf(stream,"    => Max2.swg.score            %ld\n",matches_metrics->max2_swg_score);
  tab_fprintf(stream,"  => MAPQ\n");
  tab_fprintf(stream,"    => MAPQ %d\n",matches_metrics->mapq);
}
