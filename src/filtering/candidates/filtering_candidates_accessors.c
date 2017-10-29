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
 * DESCRIPTION:
 */

#include "align/alignment.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/region/filtering_region.h"
#include "filtering/region/filtering_region_verify.h"
#include "filtering/region/filtering_region_align.h"
#include "matches/classify/matches_classify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Counting
 */
uint64_t filtering_candidates_count_regions_by_status(
    const filtering_candidates_t* const filtering_candidates,
    const filtering_region_status_t filtering_region_status) {
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  filtering_region_t** const filtering_region = filtering_candidates_get_regions(filtering_candidates);
  uint64_t count = 0, n;
  for (n=0;n<num_filtering_regions;++n) {
    if (filtering_region[n]->status == filtering_region_status) ++count;
  }
  return count;
}
/*
 * Compose filtering-region (from a group of candidate-positions)
 */
void filtering_candidates_compose_filtering_region_from_positions_exact(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const uint64_t candidate_idx,
    const uint64_t align_offset,
    const bool run_length_text) {
  // Parameters
  filtering_position_t* const candidate_position = filtering_candidates_get_positions(filtering_candidates)[candidate_idx];
  // Prepare alignment
  alignment_t* const alignment = &filtering_region->alignment;
  alignment->distance_min_bound = 0;
  const uint64_t key_length = (pattern->run_length) ? pattern->rl_key_length : pattern->key_length;
  alignment->num_tiles = 1;
  alignment->alignment_tiles = filtering_candidates_allocate_alignment_tiles(filtering_candidates,1);
  alignment->alignment_tiles->distance = 0;
  alignment->alignment_tiles->text_begin_offset = align_offset;
  alignment->alignment_tiles->text_end_offset = align_offset + key_length;
  // Compose scaffolding
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_ALIGNMENT_REGIONS_TOTAL,1);
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  match_scaffold_init(match_scaffold);
  match_scaffold->alignment_regions =
      match_scaffold_allocate_alignment_region(
          match_scaffold,1,filtering_candidates->mm_allocator);
  match_scaffold->num_alignment_regions = 1;
  match_scaffold->alignment_regions_rl = run_length_text;
  match_scaffold->scaffolding_coverage = 0;
  // Compose alignment-regions
  match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions;
  const uint64_t region_length =
      candidate_position->source_region_end -
      candidate_position->source_region_begin;
  // Read coordinates
  const uint64_t region_key_begin = candidate_position->source_region_begin;
  const uint64_t region_key_end = candidate_position->source_region_end;
  // Text coordinates (relative to the effective begin position)
  const uint64_t region_text_begin =
      candidate_position->region_text_position -
      filtering_region->text_begin_position;
  const uint64_t region_text_end = region_text_begin + region_length;
  // Init alignment-region
  match_alignment_region_init(
      match_alignment_region,match_alignment_region_exact,0,0,0,
      region_key_begin,region_key_end,region_text_begin,region_text_end);
  match_scaffold->scaffolding_coverage += region_length;
  // PROF
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_ALIGNMENT_COVERAGE,
      (100*match_scaffold->scaffolding_coverage)/
      ((pattern->run_length) ? pattern->rl_key_length : pattern->key_length));
}
void filtering_candidates_compose_filtering_region_from_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const uint64_t first_candidate_idx,
    const uint64_t last_candidate_idx,
    const bool compose_alignment_regions,
    const bool run_length_text) {
  // Parameters
  filtering_position_t** const candidate_positions = filtering_candidates_get_positions(filtering_candidates);
  // Prepare Region-Alignment
  alignment_t* const alignment = &filtering_region->alignment;
  alignment->distance_min_bound = ALIGN_DISTANCE_UNKNOWN;
  alignment->alignment_tiles = NULL;
  // Compose alignment-regions
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_ALIGNMENT_REGIONS_TOTAL,last_candidate_idx-first_candidate_idx+1);
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  match_scaffold_init(match_scaffold);
  if (compose_alignment_regions) {
    const uint64_t num_alignment_regions = last_candidate_idx-first_candidate_idx+1;
    match_scaffold->alignment_regions =
        match_scaffold_allocate_alignment_region(
            match_scaffold,num_alignment_regions,filtering_candidates->mm_allocator);
    match_scaffold->num_alignment_regions = num_alignment_regions;
    match_scaffold->alignment_regions_rl = run_length_text;
    match_scaffold->scaffolding_coverage = 0;
    uint64_t i;
    for (i=0;i<num_alignment_regions;++i) {
      match_alignment_region_t* const match_alignment_region = match_scaffold->alignment_regions + i;
      filtering_position_t* const candidate_position = candidate_positions[first_candidate_idx + i];
      // Region error
      const uint64_t region_length =
          candidate_position->source_region_end -
          candidate_position->source_region_begin;
      // Alignment-region type
      const match_alignment_region_type region_type = (candidate_position->source_region_error==0) ?
          match_alignment_region_exact : match_alignment_region_approximate;
      const uint64_t region_error = candidate_position->source_region_error;
      // Read coordinates
      const uint64_t region_key_begin = candidate_position->source_region_begin;
      const uint64_t region_key_end = candidate_position->source_region_end;
      // Text coordinates (relative to the effective begin position)
      const uint64_t region_text_begin =
          candidate_position->region_text_position -
          filtering_region->text_begin_position;
      const uint64_t region_text_end = region_text_begin + region_length;
      // Init alignment-region
      match_alignment_region_init(
          match_alignment_region,region_type,region_error,0,0,
          region_key_begin,region_key_end,region_text_begin,region_text_end);
      match_scaffold->scaffolding_coverage += region_length;
    }
    // PROF
    PROF_ADD_COUNTER(GP_CANDIDATE_REGION_ALIGNMENT_COVERAGE,
        (100*match_scaffold->scaffolding_coverage)/
        ((pattern->run_length) ? pattern->rl_key_length : pattern->key_length));
  }
}
/*
 * Adding Positions (Candidate Positions)
 */
void filtering_candidates_add_positions_from_interval(
    filtering_candidates_t* const filtering_candidates,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern,
    const uint64_t interval_lo,
    const uint64_t interval_hi,
    const uint64_t region_begin_pos,
    const uint64_t region_end_pos,
    const uint64_t region_errors) {
  // Check total candidates
  const uint64_t total_candidates = interval_hi-interval_lo;
  if (gem_expect_false(total_candidates==0)) return;
  // Check for exact matches (whole read)
  const uint64_t key_length = pattern->key_length;
  const bool exact_match = region_errors==0 && region_begin_pos==0 && region_end_pos==key_length;
  // Select matches
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  uint64_t interval_top;
  if (exact_match &&
      select_parameters->min_reported_strata_nominal==0 &&
      total_candidates > select_parameters->max_searched_matches) {
    interval_top = interval_lo + select_parameters->max_searched_matches;
  } else {
    interval_top = interval_hi;
  }
  // Store candidate positions
  uint64_t index_position;
  for (index_position=interval_lo;index_position<interval_top;++index_position) {
    // Allocate
    filtering_position_t* const filtering_position = filtering_candidates_allocate_position(filtering_candidates);
    // Configure
    filtering_position->source_region_begin = region_begin_pos;
    filtering_position->source_region_end = region_end_pos;
    filtering_position->source_region_error = region_errors;
    filtering_position->region_index_position = index_position;
    filtering_position->align_distance = exact_match ? 0 : ALIGN_DISTANCE_UNKNOWN;
  }
}
/*
 * Adding Region (filtering regions)
 */
void filtering_candidates_add_region_from_group_positions(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const uint64_t first_candidate_idx,
    const uint64_t last_candidate_idx,
    const uint64_t align_distance,
    const uint64_t align_offset,
    const bool compose_alignment_regions,
    const bool run_length_text) {
  // Parameters
  filtering_position_t** const candidate_positions = filtering_candidates_get_positions(filtering_candidates);
  // Allow new filtering-region & store it
  filtering_region_t* const filtering_region = filtering_candidates_allocate_region(filtering_candidates);
  // State
  filtering_region->status = filtering_region_unverified; // Newly created region (unverified)
  // Location
  filtering_position_t* const first_candidate = candidate_positions[first_candidate_idx];
  filtering_position_t* const last_candidate = candidate_positions[last_candidate_idx];
  filtering_region->text_trace.text = NULL; // Not retrieved
  filtering_region->text_begin_position = first_candidate->text_begin_position;
  filtering_region->text_end_position = last_candidate->text_end_position;
  // Source-region offsets
  filtering_region->text_source_region_offset = first_candidate->source_region_text_offset;
  filtering_region->key_source_region_offset = first_candidate->source_region_begin;
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_LENGTH,
      filtering_region->text_end_position-filtering_region->text_begin_position);
  // Compute key trims (if we know the real text dimensions)
  if (!run_length_text) {
    filtering_region_compute_key_trims(filtering_region,pattern);
  } else {
    filtering_region->key_trim_left = 0;
    filtering_region->key_trim_right = 0;
    filtering_region->key_trimmed = false;
  }
  // Prepare region-alignment & compose regions-matching
  if (align_distance==0) {
    filtering_candidates_compose_filtering_region_from_positions_exact(
        filtering_candidates,filtering_region,pattern,
        first_candidate_idx,align_offset,run_length_text);
  } else {
    filtering_candidates_compose_filtering_region_from_positions(
        filtering_candidates,filtering_region,pattern,
        first_candidate_idx,last_candidate_idx,
        compose_alignment_regions,run_length_text);
  }
}
void filtering_candidates_add_region_verified(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const uint64_t text_begin_offset,
    const uint64_t text_end_offset,
    const uint64_t begin_position,
    const uint64_t end_position,
    const uint64_t align_distance) {
  // Allocate new filtering-region
  filtering_region_t* const filtering_region =
      filtering_candidates_allocate_region(filtering_candidates);
  // State
  filtering_region->status = filtering_region_accepted;
  filtering_region->text_trace.text = NULL; // Not retrieved
  // Location
  filtering_region->text_source_region_offset = 0;
  filtering_region->key_source_region_offset = 0;
  filtering_region->text_begin_position = begin_position;
  filtering_region->text_end_position = end_position;
  filtering_region->key_trim_left = 0;
  filtering_region->key_trim_right = 0;
  filtering_region->key_trimmed = false;
  // Scaffolding
  match_scaffold_init(&filtering_region->match_scaffold);
  // Regions-Alignment
  alignment_t* const alignment = &filtering_region->alignment;
  alignment->distance_min_bound = align_distance;
  alignment->num_tiles = 1;
  alignment->alignment_tiles = filtering_candidates_allocate_alignment_tiles(filtering_candidates,1);
  alignment->alignment_tiles->text_begin_offset = text_begin_offset;
  alignment->alignment_tiles->text_end_offset = text_end_offset;
  alignment->alignment_tiles->distance = align_distance;
}
