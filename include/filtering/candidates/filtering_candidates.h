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
#include "filtering/region/filtering_region.h"
#include "filtering/region/filtering_region_cache.h"
#include "matches/align/match_alignment_region.h"

#define GEM_HIST_CAND_ALIGNED	64

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
  uint64_t region_index_position;        // Begin region position (index-space)
  uint64_t region_text_position;         // Begin region position (text-space)
  uint64_t align_distance;               // Align-distance if already know (NSed positions)
  // Source Region
  uint64_t source_region_begin;          // Source-region Begin
  uint64_t source_region_end;            // Source-region End
  uint64_t source_region_error;          // Source-region Error
  uint64_t source_region_text_offset;    // Source-region Text-Offset
  // Position
  uint64_t sequence_id;                  // Id of the sequence the position belongs to
  uint64_t text_begin_position;          // Region effective begin position (adjusted to error boundaries)
  uint64_t text_end_position;            // Region effective end position (adjusted to error boundaries)
} filtering_position_t;
/*
 * Filtering Candidates
 */
typedef struct {
  /* Index Structures & Parameters */
  archive_t* archive;                                              // Archive
  search_parameters_t* search_parameters;                          // Search Parameters
  /* Candidates */
  vector_t* filtering_positions;                                   // Candidate positions (filtering_position_t*)
  vector_t* filtering_regions;                                     // Candidate regions (filtering_region_t*)
  vector_t* discarded_regions;                                     // Discarded regions (filtering_region_t*)
  /* Cache */
  filtering_region_cache_t filtering_region_cache;                 // Filtering-Region Cache
  /* Stats */
  uint64_t total_candidates_analized;                              // Total number of candidates analized
  uint64_t total_candidates_canonical_gpu;                         // Minimum number of candidates to realign on GPU
  uint64_t total_candidates_realigned_swg;                         // Total number of candidates realigned (gap-affine)
  gem_counter_t candidates_aligned_histo[GEM_HIST_CAND_ALIGNED];   // Tracks candidates aligned
  /* MM */
  mm_allocator_t* mm_allocator;                    	               // MM-Allocator
} filtering_candidates_t;

/*
 * Setup
 */
void filtering_candidates_init(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_init_alignment(
    filtering_candidates_t* const filtering_candidates,
    alignment_t* const alignment,
    pattern_t* const pattern,
    const uint64_t text_length,
    const uint64_t max_error);
void filtering_candidates_clear(
    filtering_candidates_t* const filtering_candidates,
    const bool free_memory);
void filtering_candidates_destroy(
    filtering_candidates_t* const filtering_candidates,
    const bool free_memory);

/*
 * Handlers Injection (Support Data Structures)
 */
void filtering_candidates_inject_handlers(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    mm_allocator_t* const mm_allocator);

/*
 * Allocators
 */
filtering_position_t* filtering_candidates_allocate_position(
    filtering_candidates_t* const filtering_candidates);
void filtering_candidates_free_position(
    const filtering_candidates_t* const filtering_candidates,
    filtering_position_t* const filtering_position);

filtering_region_t* filtering_candidates_allocate_region(
    filtering_candidates_t* const filtering_candidates);
void filtering_candidates_free_region(
    const filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region);

alignment_tile_t* filtering_candidates_allocate_alignment_tiles(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t num_alignment_tiles);
void filtering_candidates_free_alignment_tiles(
    const filtering_candidates_t* const filtering_candidates,
    alignment_tile_t* const alignment_tile);

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
    const filtering_candidates_t* const filtering_candidates,
    const bool free_positions);

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
    const filtering_candidates_t* const filtering_candidates,
    const bool free_regions);

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
    const filtering_candidates_t* const filtering_candidates,
    const bool free_regions);

/*
 * Sorting
 */
void filtering_candidates_sort_positions(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_sort_regions_by_align_distance(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_sort_regions_by_scaffold_coverage(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_sort_discarded_by_scaffold_coverage(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_sort_discarded_by_rank(filtering_candidates_t* const filtering_candidates);

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
