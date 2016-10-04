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

#ifndef FILTERING_CANDIDATES_H_
#define FILTERING_CANDIDATES_H_

#include "utils/essentials.h"
#include "archive/archive.h"
#include "archive/search/archive_search_se_parameters.h"
#include "filtering/candidates/filtering_candidates_mm.h"
#include "filtering/region/filtering_region.h"
#include "filtering/region/filtering_region_cache.h"
#include "matches/align/match_alignment_region.h"

/*
 * Candidate Position
 */
typedef struct {
  // Locator Interval
  locator_interval_t* locator_interval;
  // Decode data
  uint64_t decode_distance;
  uint64_t decode_sampled_pos;
  // Region location
  uint64_t region_index_position;      // Begin region position (index-space)
  uint64_t region_text_position;       // Begin region position (text-space)
  uint64_t align_distance;             // Align-distance if already know (NSed positions)
  // Source Region
  uint64_t source_region_begin;        // Source-region Begin
  uint64_t source_region_end;          // Source-region End
  uint64_t source_region_error;        // Source-region Error
  uint64_t source_region_text_offset;  // Source-region Text-Offset
  // Position
  uint64_t sequence_id;                // Id of the sequence the position belongs to
  uint64_t text_begin_position;        // Region effective begin position (adjusted to error boundaries)
  uint64_t text_end_position;          // Region effective end position (adjusted to error boundaries)
} filtering_position_t;
/*
 * Filtering Candidates
 */
typedef struct {
  /* Index Structures & Parameters */
  archive_t* archive;                              // Archive
  search_parameters_t* search_parameters;          // Search Parameters
  /* Candidates */
  vector_t* filtering_positions;                   // Candidate positions (filtering_position_t*)
  vector_t* filtering_regions;                     // Candidate regions (filtering_region_t*)
  vector_t* discarded_regions;                     // Discarded regions (filtering_region_t*)
  /* Cache */
  filtering_region_cache_t filtering_region_cache; // Filtering-Region Cache
  /* Text-Collection */
  text_collection_t text_collection;               // Stores text-traces
  /* MM */
  filtering_candidates_mm_t* mm;                   // Filtering-Candidates MM
  filtering_candidates_buffered_mm_t* buffered_mm; // Filtering-Candidates MM Buffered
} filtering_candidates_t;

/*
 * Setup
 */
void filtering_candidates_init(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates);

/*
 * Handlers Injection (Support Data Structures)
 */
void filtering_candidates_inject_handlers(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    filtering_candidates_mm_t* const filtering_candidates_mm,
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm);

/*
 * Allocators
 */
filtering_position_t* filtering_candidates_allocate_position(
    filtering_candidates_t* const filtering_candidates);
filtering_region_t* filtering_candidates_allocate_region(
    filtering_candidates_t* const filtering_candidates);
filtering_region_t* filtering_candidates_allocate_discarded_region(
    filtering_candidates_t* const filtering_candidates);
alignment_tile_t* filtering_candidates_allocate_alignment_tiles(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t num_alignment_tiles);
match_alignment_region_t* filtering_candidates_allocate_alignment_regions(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t num_alignment_regions);

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
    const uint64_t max_error);
void filtering_candidates_init_alignment(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    const bool force_reset);

/*
 * Accessors
 */
uint64_t filtering_candidates_count_regions_by_status(
    const filtering_candidates_t* const filtering_candidates,
    const filtering_region_status_t filtering_region_status);

/*
 * Filtering Positions
 */
uint64_t filtering_candidates_get_num_positions(
    const filtering_candidates_t* const filtering_candidates);
void filtering_candidates_set_num_positions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_positions);

filtering_position_t** filtering_candidates_get_positions(
    const filtering_candidates_t* const filtering_candidates);
void filtering_candidates_clear_positions(
    const filtering_candidates_t* const filtering_candidates);

/*
 * Filtering Regions
 */
uint64_t filtering_candidates_get_num_regions(
    const filtering_candidates_t* const filtering_candidates);
void filtering_candidates_set_num_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_regions);

filtering_region_t** filtering_candidates_get_regions(
    const filtering_candidates_t* const filtering_candidates);
void filtering_candidates_clear_regions(
    const filtering_candidates_t* const filtering_candidates);

/*
 * Discarded Regions
 */
uint64_t filtering_candidates_get_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates);
void filtering_candidates_set_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_discarded_regions);
void filtering_candidates_add_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_discarded_regions);

filtering_region_t** filtering_candidates_get_discarded_regions(
    const filtering_candidates_t* const filtering_candidates);
filtering_region_t** filtering_candidates_reserve_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_regions);
void filtering_candidates_clear_discarded_regions(
    const filtering_candidates_t* const filtering_candidates);

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
    bool* const candidates_limited);

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
    const bool run_length_text);
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
    const uint64_t text_end_offset);

/*
 * Sorting
 */
void filtering_candidates_sort_regions_by_align_distance(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_sort_regions_by_scaffold_coverage(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_sort_positions(filtering_candidates_t* const filtering_candidates);

/*
 * Display
 */
void filtering_candidates_print_regions(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates,
    const bool print_alignment_regions);
void filtering_candidates_print_positions(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates);

#endif /* FILTERING_CANDIDATES_H_ */
