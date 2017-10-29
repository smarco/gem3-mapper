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

#include "filtering/candidates/filtering_candidates_classify.h"

/*
 * Filtering candidates classify as subdominant match
 * (Expected best score worse than all the found matches so far)
 */
bool filtering_candidates_classify_subdominant_match_edit(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  alignment_t* const alignment = &filtering_region->alignment;
  const uint64_t max_searched_matches = select_parameters->max_searched_matches;
  // The candidate needs to have a expected max-score than the current max
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  const uint64_t candidate_edit_distance_bound = alignment->distance_min_bound;
  const uint64_t worst_reported_match_edit_distance = match_traces[max_searched_matches-1]->edit_distance; // FIXME Sort by SWG
  return candidate_edit_distance_bound >= worst_reported_match_edit_distance;
}
bool filtering_candidates_classify_subdominant_match_swg(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  swg_penalties_t* const swg_penalties = &search_parameters->swg_penalties;
  alignment_t* const alignment = &filtering_region->alignment;
  const uint64_t max_searched_matches = select_parameters->max_searched_matches;
  // The candidate needs to have a expected max-score than the current max
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  const uint64_t candidate_edit_distance_bound = alignment->distance_min_bound;
  const int32_t candidate_max_score_bound = align_swg_score_compute_max_score_bound(
      swg_penalties,candidate_edit_distance_bound,pattern->key_length);
  const int32_t worst_reported_match_score = match_traces[max_searched_matches-1]->swg_score;
  return candidate_max_score_bound <= worst_reported_match_score;
}
bool filtering_candidates_classify_subdominant_match(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  // Basic cases
  const uint64_t min_reported_strata = select_parameters->min_reported_strata_nominal;
  const uint64_t max_searched_matches = select_parameters->max_searched_matches;
  if (min_reported_strata > 0) return false;
  if (num_matches == 0 || num_matches < max_searched_matches) return false;
  // Score/Distance Bounded Cases (Only pays off to align matches
  //   that can be include within user report limits)
  switch (search_parameters->match_alignment_model) {
    case match_alignment_model_hamming:
    case match_alignment_model_levenshtein:
      return filtering_candidates_classify_subdominant_match_edit(
          filtering_candidates,filtering_region,pattern,matches);
      break;
    case match_alignment_model_gap_affine:
      return filtering_candidates_classify_subdominant_match_swg(
          filtering_candidates,filtering_region,pattern,matches);
      break;
    default:
      return false;
      break;
  }
}
/*
 * Filtering candidates classify as subdominant region
 * compared with a reference score/distance
 */
bool filtering_candidates_classify_subdominant_region_edit(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t filtering_region_rank,
    const uint64_t filtering_region_edit_bound,
    const uint64_t reference_edit_bound) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  // Basic cases
  const uint64_t min_reported_strata = select_parameters->min_reported_strata_nominal;
  const uint64_t max_searched_matches = select_parameters->max_searched_matches;
  if (min_reported_strata > 0) return false;
  if (filtering_region_rank == 0 || filtering_region_rank < max_searched_matches) return false;
  // Score Bounded Case
  return filtering_region_edit_bound >= reference_edit_bound;
}

bool filtering_candidates_classify_subdominant_region_swg(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t filtering_region_rank,
    const int32_t filtering_region_score_bound,
    const int32_t reference_score_bound) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  // Basic cases
  const uint64_t min_reported_strata = select_parameters->min_reported_strata_nominal;
  const uint64_t max_searched_matches = select_parameters->max_searched_matches;
  if (min_reported_strata > 0) return false;
  if (filtering_region_rank == 0 || filtering_region_rank < max_searched_matches) return false;
  // Score Bounded Case
  return filtering_region_score_bound <= reference_score_bound;
}



