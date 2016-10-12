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

#include "matches/classify/matches_classify.h"
#include "matches/matches_cigar.h"

/*
 * Classify
 */
void matches_classify(matches_t* const matches) {
  // Parameters
  matches_metrics_t* const metrics = &matches->metrics;
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  match_trace_t* const primary_match = match_traces[0];
  match_trace_t* const subdominant_match = match_traces[1];
  // Classify
  if (metrics->accepted_matches == 0) {
    matches->matches_class = matches_class_unmapped;
  } else if (metrics->accepted_matches == 1) {
    matches->matches_class = matches_class_unique;
  } else if (primary_match->swg_score == subdominant_match->swg_score ||
             matche_trace_cigar_cmp(matches->cigar_vector,primary_match,
                                    matches->cigar_vector,subdominant_match)==0) {
    matches->matches_class = matches_class_tie_perfect; // Perfect tie
  } else if (primary_match->event_distance > metrics->min_event_distance ||
             metrics->min_event_distance_count > 1 ||
             primary_match->edit_distance > metrics->min_edit_distance ||
             metrics->min_edit_distance_count > 1) {
    matches->matches_class = matches_class_tie; // General tie
  } else if (primary_match->edit_distance+1 == subdominant_match->edit_distance) {
    matches->matches_class = matches_class_mmap_d1; // MultiMap-D1
  } else {
    matches->matches_class = matches_class_mmap; // MultiMap
  }
}
void paired_matches_classify(paired_matches_t* const paired_matches) {
  // Parameters
  matches_metrics_t* const metrics = &paired_matches->metrics;
  // Unmapped
  if (metrics->accepted_matches == 0) {
    paired_matches->paired_matches_class = paired_matches_class_unmapped;
    return;
  }
  // Unique
  if (metrics->accepted_matches == 1) {
    paired_matches->paired_matches_class = paired_matches_class_unique;
    return;
  }
  // Tie-d0
  paired_map_t* const primary_match = paired_matches_get_maps(paired_matches);
  paired_map_t* const subdominant_match = primary_match + 1;
  if (primary_match->swg_score == subdominant_match->swg_score &&
      primary_match->template_length_sigma == subdominant_match->template_length_sigma) {
    paired_matches->paired_matches_class = paired_matches_class_tie_perfect;
    return;
  }
  // Tie-d1
  if (primary_match->swg_score == subdominant_match->swg_score ||
      primary_match->event_distance > metrics->min_event_distance ||
      metrics->min_event_distance_count > 1 ||
      primary_match->edit_distance > metrics->min_edit_distance ||
      metrics->min_edit_distance_count > 1) {
    paired_matches->paired_matches_class = paired_matches_class_tie;
    return;
  }
  // MultiMap-D1
  if (primary_match->edit_distance+1 == subdominant_match->edit_distance) {
    paired_matches->paired_matches_class = paired_matches_class_mmap_d1;
    return;
  }
  // MultiMap
  paired_matches->paired_matches_class = paired_matches_class_mmap;
}
