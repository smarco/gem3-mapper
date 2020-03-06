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

#ifndef MATCHES_H_
#define MATCHES_H_

#include "utils/essentials.h"
#include "archive/search/archive_select_parameters.h"
#include "archive/locator.h"
#include "matches/classify/matches_metrics.h"
#include "matches/matches_counters.h"
#include "matches/match_trace.h"

/*
 * Matches Classes
 */
typedef enum {
  matches_class_unmapped    = 0,
  matches_class_mmap        = 1,
  matches_class_unique      = 2,
} matches_class_type;
typedef struct {
  matches_class_type matches_class; // Matches class
  int64_t delta_group;              // Delta-strata group
  int64_t wdelta_group;             // Actual delta-strata group (considering maximum complete strata)
} matches_classification_t;

/*
 * Match Blocks
 * (non-overlapping regions of read with local alignments - a global alignment will create a single block stretching the length of the read)
 */

typedef struct {
	uint64_t lo;              // Start position of block in read
	uint64_t hi;						  // End position of block in read
	uint64_t clip_right;      // clip_left == lo
	uint64_t num_matches;     // Number of matches to this block
	bool overlaps;
	bool printed;
} match_block_t;

/*
 * Matches
 */
typedef struct {
  /* Matches Classification */
  matches_classification_t classification; //  Matches Classification
  /* Matches Counters */
  matches_counters_t* counters;            // Global counters
  uint64_t max_complete_stratum;           // Maximum complete stratum
  uint64_t key_length;										 // Length of pattern;
  bool limited_exact_matches;              // Limited exact matches
  bool matches_extended;                   // Matches added from PE-extension of the other end
  bool local_search;                       // Are we in global or local search phase?
  /* Matches */
  vector_t* match_traces;                  // Matches (match_trace_t*)
  vector_t* match_traces_local;            // Local Matches (match_trace_t*)
  vector_t* match_traces_extended;         // Extended Matches (match_trace_t*)
  ihash_t* match_traces_begin;             // Begin position (of the aligned match) in the text-space
  ihash_t* match_traces_end;               // End position (of the aligned match) in the text-space
  /* Blocks */
  vector_t* match_blocks;                  // Almost non-overlapping blocks of the read with local alignments
  uint64_t match_total_covered_bases;      // Number of bases covered by blocks
  /* Identified restrictions sites */
  vector_t* match_potential_split_sites;   //Iindexes into read of potential split sites (uint64_t)
  bool match_replaced;                     // A match has been replaced (can affect paired-end metrics)
  /* CIGAR */
  vector_t* cigar_vector;                  // CIGAR operations storage (cigar_element_t)
  /* Metrics */
  matches_metrics_t metrics;               // Metrics
  /* MM */
  mm_slab_t* mm_slab;                      // MM-Slab
  mm_allocator_t* mm_allocator;            // MM-Allocator
} matches_t;

/*
 * Setup
 */
matches_t* matches_new(void);
void matches_clear(matches_t* const matches);
void matches_delete(matches_t* const matches);

/*
 * Accessors
 */
bool matches_is_mapped(const matches_t* const matches);
void matches_recompute_metrics(matches_t* const matches);
uint64_t matches_get_first_stratum_matches(matches_t* const matches);
uint64_t matches_get_subdominant_stratum_matches(matches_t* const matches);
uint8_t matches_get_primary_mapq(matches_t* const matches);

void matches_update_mcs(
    matches_t* const matches,
    const uint64_t current_mcs);
void matches_update_limited_exact_matches(
    matches_t* const matches,
    const uint64_t num_exact_matches_limited);

/*
 * Matches Accessors
 */
void matches_clear_match_traces(const matches_t* const matches);
uint64_t matches_get_num_match_traces(const matches_t* const matches);
uint64_t matches_get_num_match_traces_extended(const matches_t* const matches);
match_trace_t* matches_get_primary_match(const matches_t* const matches);
match_trace_t** matches_get_match_traces(const matches_t* const matches);

/*
 * Sort
 */
void matches_traces_sort_by_genomic_position(
    match_trace_t** const match_traces,
    const uint64_t num_match_traces);
void matches_traces_sort_by_primary_and_read_position(
    match_trace_t** const match_traces,
    const uint64_t num_match_traces);

/*
 * Adding Match-traces
 */
match_trace_t* matches_add_match_trace(
		search_parameters_t* const search,
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace,
    bool* const match_replaced);
void matches_add_match_trace_extended(
		search_parameters_t* const search,
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace);
void matches_add_match_trace_local_pending(
    matches_t* const matches,
    match_trace_t* const match_trace);

void matches_local_pending_add_to_regular_matches(
		search_parameters_t* const search,
    matches_t* const matches,
    const locator_t* const locator);
void matches_local_pending_add_to_extended_matches(
    matches_t* const matches,
    const locator_t* const locator);
uint64_t matches_num_match_blocks(
		const matches_t * const matches);

/*
 * Display
 */
void matches_print(
    FILE* const stream,
    matches_t* const matches);

#endif /* MATCHES_H_ */
