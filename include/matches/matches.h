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
#include "text/text_collection.h"
#include "align/align_swg_score.h"
#include "archive/locator.h"
#include "archive/search/archive_select_parameters.h"
#include "matches/matches_counters.h"
#include "matches/align/match_alignment.h"
#include "matches/matches_cigar.h"
#include "matches/classify/matches_metrics.h"

/*
 * Matches Type
 */
typedef enum {
  match_type_regular,           // Regular Match
  match_type_local,             // Local Match (product of local alignment)
  match_type_extended           // Extended Match (product of extending PE)
} match_type;

/*
 * Matches Classes
 */
typedef enum {
  matches_class_unmapped    = 0,
  matches_class_tie_perfect = 1,
  matches_class_tie         = 2,
  matches_class_mmap_d1     = 3,
  matches_class_mmap        = 4,
  matches_class_unique      = 5,
} matches_class_t;
extern const char* matches_class_label[6];

/*
 * Match (Trace-Match)
 */
typedef struct {
  /* Type */
  match_type type;                   // Match type
  /* Location */
  char* sequence_name;               // Sequence name (After decoding.Eg Chr1)
  strand_t strand;                   // Mapping Strand
  bs_strand_t bs_strand;             // Bisulfite Strand
  uint64_t text_position;            // Position of the match in the text. Local text (Eg wrt Chr1)
  /* Reference-Text */
  uint64_t text_trace_offset;        // Trace-offset in the text-collection
  uint8_t* text;                     // Pointer to the matching-text
  uint64_t text_length;              // Length of the matching-text
  /* Distance/Score */
  uint64_t event_distance;           // Distance
  uint64_t edit_distance;            // Edit-Distance
  int32_t swg_score;                 // SWG Distance/Score
  uint8_t mapq_score;                // MAPQ Score
  /* Alignment */
  match_alignment_t match_alignment; // Match Alignment (CIGAR + ...)
  void* match_scaffold;              // Supporting Scaffolding
} match_trace_t;

/*
 * Matches
 */
typedef struct {
  /* State */
  matches_class_t matches_class;
  uint64_t max_complete_stratum;
  /* Matches Counters */
  matches_counters_t* counters;        // Global counters
  /* Matches */
  vector_t* match_traces;              // Matches (match_trace_t*)
  vector_t* match_traces_local;        // Local Matches (match_trace_t)
  ihash_t* match_traces_begin;         // Begin position (of the aligned match) in the text-space
  ihash_t* match_traces_end;           // End position (of the aligned match) in the text-space
  /* CIGAR */
  vector_t* cigar_vector;              // CIGAR operations storage (cigar_element_t)
  /* Metrics */
  matches_metrics_t metrics;           // Metrics
  /* MM */
  mm_slab_t* mm_slab;                  // MM-Slab
  mm_stack_t* mm_stack;                // MM-Stack
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

void matches_update_mcs(matches_t* const matches,const uint64_t current_mcs);

/*
 * Matches Accessors
 */
void matches_clear_match_traces(const matches_t* const matches);
uint64_t matches_get_num_match_traces(const matches_t* const matches);
match_trace_t* matches_get_primary_match(const matches_t* const matches);
match_trace_t** matches_get_match_traces(const matches_t* const matches);

/*
 * Match-Trace
 */
cigar_element_t* match_trace_get_cigar_buffer(
    const matches_t* const matches,
    const match_trace_t* const match_trace);
uint64_t match_trace_get_cigar_length(const match_trace_t* const match_trace);
uint64_t match_trace_get_event_distance(const match_trace_t* const match_trace);
int64_t match_trace_get_effective_length(
    matches_t* const matches,
    const uint64_t read_length,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);
int matche_trace_cigar_cmp(
    vector_t* const cigar_vector_match0,
    match_trace_t* const match0,
    vector_t* const cigar_vector_match1,
    match_trace_t* const match1);

/*
 * Add Match-Trace
 */
match_trace_t* matches_add_match_trace(
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace,
    bool* const match_replaced);

/*
 * Sort
 */
void matches_traces_sort_by_genomic_position(
    match_trace_t** const match_traces,
    const uint64_t num_match_traces);

/*
 * Local Matches
 */
void matches_add_local_match_pending(
    matches_t* const matches,
    match_trace_t* const match_trace);
void matches_add_pending_local_matches(
    matches_t* const matches,
    const locator_t* const locator);

/*
 * Filters
 */
void matches_filter_by_mapq(
    matches_t* const matches,
    const uint8_t mapq_threshold);

/*
 * Display
 */
void matches_print(
    FILE* const stream,
    matches_t* const matches);

#endif /* MATCHES_H_ */
