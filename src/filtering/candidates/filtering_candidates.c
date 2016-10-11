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
 *   Filtering candidates module provides functions to store and handle
 *   all positions/regions during a search based on filtering, that is,
 *   generation of candidates & verification of candidates
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
 * Constants
 */
#define REGIONS_BUFFER_INIT                      100
#define CANDIDATE_POSITIONS_INIT                 1000

/*
 * Setup
 */
void filtering_candidates_init(filtering_candidates_t* const filtering_candidates) {
  // Candidates
  filtering_candidates->filtering_positions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_position_t*);
  filtering_candidates->filtering_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t*);
  filtering_candidates->discarded_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t*);
  // Cache
  filtering_region_cache_init(&filtering_candidates->filtering_region_cache);
  // Text-Collection
  text_collection_init(&filtering_candidates->text_collection);
}
void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates) {
  // Candidates
  vector_clear(filtering_candidates->filtering_positions);
  vector_clear(filtering_candidates->filtering_regions);
  vector_clear(filtering_candidates->discarded_regions);
  // Text-Collection
  text_collection_clear(&filtering_candidates->text_collection);
}
void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates) {
  // Candidates
  vector_delete(filtering_candidates->filtering_positions);
  vector_delete(filtering_candidates->filtering_regions);
  vector_delete(filtering_candidates->discarded_regions);
  // Text-Collection
  text_collection_destroy(&filtering_candidates->text_collection);
}
/*
 * Handlers Injection (Support Data Structures)
 */
void filtering_candidates_inject_handlers(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    filtering_candidates_mm_t* const filtering_candidates_mm,
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm) {
  filtering_candidates->archive = archive;
  filtering_candidates->search_parameters = search_parameters;
  filtering_candidates->mm = filtering_candidates_mm;
  filtering_candidates->buffered_mm = filtering_candidates_buffered_mm;
  if (filtering_candidates_buffered_mm != NULL) {
    filtering_candidates->mm->mm_alignment_regions = filtering_candidates_buffered_mm->mm_regions_attr;
    filtering_candidates->mm->mm_alignment_tiles = filtering_candidates_buffered_mm->mm_regions_attr;
  }
  text_collection_inject_mm(&filtering_candidates->text_collection,filtering_candidates->mm->mm_text);
}
/*
 * Allocators
 */
filtering_position_t* filtering_candidates_allocate_position(
    filtering_candidates_t* const filtering_candidates) {
  filtering_position_t* const position = mm_stack_alloc(filtering_candidates->mm->mm_positions,filtering_position_t);
  vector_insert(filtering_candidates->filtering_positions,position,filtering_position_t*);
  return position;
}
filtering_region_t* filtering_candidates_allocate_region(
    filtering_candidates_t* const filtering_candidates) {
  filtering_region_t* const region = mm_stack_alloc(filtering_candidates->mm->mm_regions,filtering_region_t);
  vector_insert(filtering_candidates->filtering_regions,region,filtering_region_t*);
  return region;
}
filtering_region_t* filtering_candidates_allocate_discarded_region(
    filtering_candidates_t* const filtering_candidates) {
  filtering_region_t* const region = mm_stack_alloc(filtering_candidates->mm->mm_regions,filtering_region_t);
  vector_insert(filtering_candidates->discarded_regions,region,filtering_region_t*);
  return region;
}
alignment_tile_t* filtering_candidates_allocate_alignment_tiles(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t num_alignment_tiles) {
  return mm_stack_calloc(filtering_candidates->mm->mm_alignment_tiles,
      num_alignment_tiles,alignment_tile_t,false);
}
match_alignment_region_t* filtering_candidates_allocate_alignment_regions(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t num_alignment_regions) {
  return mm_stack_calloc(filtering_candidates->mm->mm_alignment_regions,
      num_alignment_regions,match_alignment_region_t,false);
}
/*
 * Prepare Alignment
 */
void filtering_candidates_init_alignment_tiles(
    filtering_candidates_t* const filtering_candidates,
    alignment_t* const alignment,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    const uint64_t text_begin_offset,
    const uint64_t text_end_offset,
    const uint64_t max_error) {
  // Allocate alignment
  const uint64_t num_tiles = bpm_pattern_tiles->num_pattern_tiles;
  alignment->num_tiles = num_tiles;
  alignment->distance_min_bound = bpm_pattern->pattern_length;
  if (alignment->alignment_tiles==NULL) {
    alignment->alignment_tiles =
        filtering_candidates_allocate_alignment_tiles(filtering_candidates,num_tiles);
  }
  // Init all tiles
  const uint64_t text_length = text_end_offset-text_begin_offset;
  alignment_tile_t* const alignment_tiles = alignment->alignment_tiles;
  if (num_tiles==1) {
    alignment_tiles->distance = ALIGN_DISTANCE_INF;
    alignment_tiles->text_end_offset = text_end_offset;
    alignment_tiles->text_begin_offset = text_begin_offset;
  } else {
    // Calculate tile dimensions
    const uint64_t effective_max_error = MIN(max_error,bpm_pattern->pattern_length);
    pattern_tiled_t pattern_tiled;
    pattern_tiled_init(&pattern_tiled,bpm_pattern->pattern_length,
        bpm_pattern_tiles->tile_length,text_length,effective_max_error);
    uint64_t tile_pos;
    for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
      // Init Tile
      alignment_tiles[tile_pos].distance = ALIGN_DISTANCE_INF;
      alignment_tiles[tile_pos].text_end_offset = text_begin_offset+pattern_tiled.tile_offset+pattern_tiled.tile_wide;
      alignment_tiles[tile_pos].text_begin_offset = text_begin_offset+pattern_tiled.tile_offset;
      // Calculate next tile
      pattern_tiled_calculate_next(&pattern_tiled);
    }
  }
}
void filtering_candidates_init_alignment(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    const bool force_reset) {
  // Check alignment
  alignment_t* const alignment = &filtering_region->alignment;
  if (alignment->alignment_tiles!=NULL && !force_reset) return; // Already initialized
  // Prepare alignment
  const uint64_t text_length = filtering_region->text_end_position-filtering_region->text_begin_position;
  filtering_candidates_init_alignment_tiles(filtering_candidates,alignment,
      bpm_pattern,bpm_pattern_tiles,0,text_length,filtering_region->max_error);
}
/*
 * Accessors
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
 * Filtering Positions
 */
uint64_t filtering_candidates_get_num_positions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->filtering_positions);
}
void filtering_candidates_set_num_positions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_positions) {
  vector_set_used(filtering_candidates->filtering_positions,num_positions);
}
filtering_position_t** filtering_candidates_get_positions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_mem(filtering_candidates->filtering_positions,filtering_position_t*);
}
void filtering_candidates_clear_positions(
    const filtering_candidates_t* const filtering_candidates) {
  vector_clear(filtering_candidates->filtering_positions);
}
/*
 * Filtering Regions
 */
uint64_t filtering_candidates_get_num_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->filtering_regions);
}
void filtering_candidates_set_num_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_regions) {
  vector_set_used(filtering_candidates->filtering_regions,num_regions);
}
filtering_region_t** filtering_candidates_get_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t*);
}
void filtering_candidates_clear_regions(
    const filtering_candidates_t* const filtering_candidates) {
  vector_clear(filtering_candidates->filtering_regions);
}
/*
 * Discarded Regions
 */
uint64_t filtering_candidates_get_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->discarded_regions);
}
void filtering_candidates_set_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_discarded_regions) {
  vector_set_used(filtering_candidates->discarded_regions,num_discarded_regions);
}
void filtering_candidates_add_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_discarded_regions) {
  vector_add_used(filtering_candidates->discarded_regions,num_discarded_regions);
}
filtering_region_t** filtering_candidates_get_discarded_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_mem(filtering_candidates->discarded_regions,filtering_region_t*);
}
filtering_region_t** filtering_candidates_reserve_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_regions) {
  vector_reserve_additional(filtering_candidates->discarded_regions,num_regions);
  return vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t*);
}
void filtering_candidates_clear_discarded_regions(
    const filtering_candidates_t* const filtering_candidates) {
  vector_clear(filtering_candidates->discarded_regions);
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
  match_scaffold->alignment_regions = filtering_candidates_allocate_alignment_regions(filtering_candidates,1);
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
  alignment->distance_min_bound = ALIGN_DISTANCE_INF;
  alignment->alignment_tiles = NULL;
  // Compose alignment-regions
  PROF_ADD_COUNTER(GP_CANDIDATE_REGION_ALIGNMENT_REGIONS_TOTAL,last_candidate_idx-first_candidate_idx+1);
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  match_scaffold_init(match_scaffold);
  if (compose_alignment_regions) {
    const uint64_t num_alignment_regions = last_candidate_idx-first_candidate_idx+1;
    match_scaffold->alignment_regions =
        filtering_candidates_allocate_alignment_regions(filtering_candidates,num_alignment_regions);
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
    const uint64_t region_errors,
    bool* const candidates_limited) {
  // Check total candidates
  const uint64_t total_candidates = interval_hi-interval_lo;
  if (gem_expect_false(total_candidates==0)) return;
  // Check for exact matches (whole read)
  const uint64_t key_length = pattern->key_length;
  const bool exact_match = region_errors==0 && region_begin_pos==0 && region_end_pos==key_length;
  // Select matches
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_align;
  uint64_t interval_top;
  if (exact_match &&
      select_parameters->min_reported_strata_nominal==0 &&
      total_candidates > select_parameters->max_reported_matches) {
    interval_top = interval_lo + select_parameters->max_reported_matches;
    *candidates_limited = true;
  } else {
    interval_top = interval_hi;
    *candidates_limited = false;
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
    filtering_position->align_distance = exact_match ? 0 : ALIGN_DISTANCE_INF;
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
  filtering_region->text_trace_offset = UINT64_MAX; // Unassigned
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
  // Set max-error & max-bandwidth
  filtering_region->max_error = pattern->max_effective_filtering_error;
  filtering_region->max_bandwidth = pattern->max_effective_bandwidth;
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
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    const uint64_t text_trace_offset,
    const uint64_t begin_position,
    const uint64_t end_position,
    const uint64_t align_distance,
    const uint64_t max_effective_bandwidth,
    const uint64_t text_begin_offset,
    const uint64_t text_end_offset) {
  // Allocate new filtering-region
  filtering_region_t* const filtering_region = filtering_candidates_allocate_region(filtering_candidates);
  // State
  filtering_region->status = filtering_region_accepted;
  // Text-trace
  filtering_region->text_trace_offset = text_trace_offset;
  // Location
  filtering_region->text_source_region_offset = 0;
  filtering_region->key_source_region_offset = 0;
  filtering_region->text_begin_position = begin_position;
  filtering_region->text_end_position = end_position;
  filtering_region->key_trim_left = 0;
  filtering_region->key_trim_right = 0;
  filtering_region->key_trimmed = false;
  // Trimmed Pattern
  filtering_region->max_error = align_distance;
  filtering_region->max_bandwidth = max_effective_bandwidth;
  filtering_region->bpm_pattern_trimmed = NULL;
  filtering_region->bpm_pattern_trimmed_tiles = NULL;
  // Scaffolding
  match_scaffold_init(&filtering_region->match_scaffold);
  // Regions-Alignment
  alignment_t* const alignment = &filtering_region->alignment;
  alignment->alignment_tiles = NULL;
  filtering_candidates_init_alignment_tiles(
      filtering_candidates,alignment,bpm_pattern,
      bpm_pattern_tiles,text_begin_offset,text_end_offset,align_distance);
}
/*
 * Sorting
 */
int filtering_region_cmp_align_distance(filtering_region_t** a,filtering_region_t** b) {
  return (*a)->alignment.distance_min_bound - (*b)->alignment.distance_min_bound;
}
int filtering_region_cmp_scaffold_coverage(filtering_region_t** a,filtering_region_t** b) {
  return (*b)->match_scaffold.scaffolding_coverage - (*a)->match_scaffold.scaffolding_coverage;
}
int64_t filtering_position_cmp_position(filtering_position_t** a,filtering_position_t** b) {
  return (int64_t)(*a)->text_end_position - (int64_t)(*b)->text_end_position;
}
#define VECTOR_SORT_NAME                 align_distance
#define VECTOR_SORT_TYPE                 filtering_region_t*
#define VECTOR_SORT_CMP(a,b)             filtering_region_cmp_align_distance(a,b)
#include "utils/vector_sort.h"
void filtering_candidates_sort_regions_by_align_distance(filtering_candidates_t* const filtering_candidates) {
  vector_sort_align_distance(filtering_candidates->filtering_regions);
}
#define VECTOR_SORT_NAME                 scaffold_coverage
#define VECTOR_SORT_TYPE                 filtering_region_t*
#define VECTOR_SORT_CMP(a,b)             filtering_region_cmp_scaffold_coverage(a,b)
#include "utils/vector_sort.h"
void filtering_candidates_sort_regions_by_scaffold_coverage(filtering_candidates_t* const filtering_candidates) {
  vector_sort_scaffold_coverage(filtering_candidates->filtering_regions);
}
#define VECTOR_SORT_NAME                 filtering_positions
#define VECTOR_SORT_TYPE                 filtering_position_t*
#define VECTOR_SORT_CMP(a,b)             filtering_position_cmp_position(a,b)
#include "utils/vector_sort.h"
void filtering_candidates_sort_positions(filtering_candidates_t* const filtering_candidates) {
  vector_sort_filtering_positions(filtering_candidates->filtering_positions);
}
/*
 * Display
 */
void filtering_candidates_print_regions_by_status(
    FILE* const stream,
    vector_t* const filtering_regions,
    const filtering_region_status_t status,
    const text_collection_t* const text_collection,
    const bool print_alignment_regions) {
  uint64_t i, total_printed = 0;
  const uint64_t num_regions = vector_get_used(filtering_regions);
  filtering_region_t** const fregion = vector_get_mem(filtering_regions,filtering_region_t*);
  // Count
  for (i=0;i<num_regions;++i) {
    if (fregion[i]->status!=status) continue;
    ++total_printed;
  }
  if (total_printed == 0) return;
  tab_fprintf(stream,"  => Regions.%s  (%"PRIu64")\n",filtering_region_status_label[status],total_printed);
  // Print
  tab_global_inc();
  for (i=0;i<num_regions;++i) {
    if (fregion[i]->status!=status) continue;
    filtering_region_print(stream,fregion[i],text_collection,false,print_alignment_regions,true);
  }
  tab_global_dec();
}
void filtering_candidates_print_regions(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates,
    const bool print_alignment_regions) {
  tab_fprintf(stream,"[GEM]>Filtering.Regions\n");
  text_collection_t* const text_collection = &filtering_candidates->text_collection;
  vector_t* const filtering_regions = filtering_candidates->filtering_regions;
  vector_t* const discarded_regions = filtering_candidates->discarded_regions;
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_pending,text_collection,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_unverified,text_collection,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,discarded_regions,filtering_region_verified_discarded,text_collection,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_accepted,text_collection,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,discarded_regions,filtering_region_accepted_subdominant,text_collection,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_aligned,text_collection,print_alignment_regions);
  const uint64_t total_regions =
      filtering_candidates_get_num_regions(filtering_candidates) +
      vector_get_used(filtering_candidates->discarded_regions);
  if (total_regions > 0) tab_fprintf(stream,"  => Total.Regions %"PRIu64"\n",total_regions);
}
void filtering_candidates_print_positions(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates) {
  const uint64_t num_candidate_positions = filtering_candidates_get_num_positions(filtering_candidates);
  filtering_position_t** const candidate_positions = filtering_candidates_get_positions(filtering_candidates);
  uint64_t i;
  tab_fprintf(stream,"[GEM]>Filtering.Positions\n");
  for (i=0;i<num_candidate_positions;++i) {
    tab_fprintf(stream,"  => Key[%lu,%lu) ~> Text[%lu,%lu) Distance=%lu\n",
        candidate_positions[i]->source_region_begin,
        candidate_positions[i]->source_region_end,
        candidate_positions[i]->text_begin_position,
        candidate_positions[i]->text_end_position,
        candidate_positions[i]->align_distance);
  }
}
