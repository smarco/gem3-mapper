/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_CANDIDATES_H_
#define FILTERING_CANDIDATES_H_

#include "utils/essentials.h"
#include "data_structures/interval_set.h"
#include "archive/archive.h"
#include "archive/archive_search_parameters.h"
#include "filtering/filtering_region.h"
#include "filtering/filtering_region_cache.h"

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
typedef struct {
  // Source Region
  uint64_t source_region_begin;        // Source-region Begin
  uint64_t source_region_end;          // Source-region End
} filtering_position_buffered_t;
/*
 * Filtering Candidates
 */
typedef struct {
  /* Index Structures & Parameters */
  archive_t* archive;                              // Archive
  search_parameters_t* search_parameters;          // Search Parameters
  /* Candidates */
  vector_t* filtering_positions;                   // Candidate positions (filtering_position_t)
  vector_t* filtering_regions;                     // Candidate regions (filtering_region_t)
  vector_t* discarded_regions;                     // Discarded regions (filtering_region_t)
  vector_t* verified_regions;                      // Verified regions (verified_region_t)
  /* Cache */
  filtering_region_cache_t filtering_region_cache; // Filtering-Region Cache
  /* Text-Collection */
  text_collection_t* text_collection;              // Stores text-traces
  /* MM */
  mm_stack_t* mm_stack;                            // MM-Stack
} filtering_candidates_t;

/*
 * Setup
 */
void filtering_candidates_init(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_clear(filtering_candidates_t* const filtering_candidates);
void filtering_candidates_destroy(filtering_candidates_t* const filtering_candidates);

/*
 * Memory Injection (Support Data Structures)
 */
void filtering_candidates_inject_search(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,
    search_parameters_t* const search_parameters);
void filtering_candidates_inject_mm_stack(
    filtering_candidates_t* const filtering_candidates,
    mm_stack_t* const mm_stack);
void filtering_candidates_inject_text_collection(
    filtering_candidates_t* const filtering_candidates,
    text_collection_t* const text_collection);

/*
 * Accessors
 */
uint64_t filtering_candidates_get_num_candidate_positions(
    const filtering_candidates_t* const filtering_candidates);
uint64_t filtering_candidates_get_num_candidate_regions(
    const filtering_candidates_t* const filtering_candidates);
uint64_t filtering_candidates_count_candidate_regions(
    filtering_candidates_t* const filtering_candidates_end,
    const filtering_region_status_t filtering_region_status);

/*
 * Adding candidate positions
 */
void filtering_candidates_add_read_interval(
    filtering_candidates_t* const filtering_candidates,
    search_parameters_t* const search_parameters,
    const uint64_t interval_lo,
    const uint64_t interval_hi,
    const uint64_t key_length,
    const uint64_t align_distance);

void filtering_candidates_add_region_interval(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t interval_lo,
    const uint64_t interval_hi,
    const uint64_t region_begin_pos,
    const uint64_t region_end_pos,
    const uint64_t region_errors);
void filtering_candidates_add_region_interval_set(
    filtering_candidates_t* const filtering_candidates,
    interval_set_t* const interval_set,
    const uint64_t region_begin_pos,
    const uint64_t region_end_pos);
void filtering_candidates_add_region_interval_set_thresholded(
    filtering_candidates_t* const filtering_candidates,
    interval_set_t* const interval_set,
    const uint64_t region_begin_pos,
    const uint64_t region_end_pos,
    const uint64_t max_error);

/*
 * Sorting
 */
void filtering_positions_sort_positions(vector_t* const filtering_positions);
void filtering_regions_sort_align_distance(vector_t* const filtering_regions);
void filtering_regions_sort_scaffold_coverage(vector_t* const filtering_regions);
void verified_regions_sort_positions(vector_t* const verified_regions);

/*
 * Display
 */
void filtering_candidates_print_regions(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates,
    const bool print_matching_regions);

#endif /* FILTERING_CANDIDATES_H_ */
