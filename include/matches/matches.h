/*
 * PROJECT: GEMMapper
 * FILE: matches.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure to store alignment matches {sequence,position,strand,...}
 */

#ifndef MATCHES_H_
#define MATCHES_H_

#include "utils/essentials.h"
#include "data_structures/interval_set.h"
#include "data_structures/text_collection.h"
#include "archive/locator.h"
#include "archive/archive_select_parameters.h"
#include "matches/matches_counters.h"
#include "matches/match_alignment.h"
#include "matches/matches_metrics.h"

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
  matches_class_unmapped = 0,
  matches_class_tie_d0 = 1,
  matches_class_tie_d1 = 2,
  matches_class_mmap = 3,
  matches_class_unique = 4,
} matches_class_t;
extern const char* matches_class_label[5];

/*
 * Match (Trace-Match)
 */
typedef struct {
  /* Type */
  match_type type;                   // Match type
  uint64_t* match_trace_offset;      // Match-Trace offset in the the matches vector
  /* Location */
  char* sequence_name;               // Sequence name (After decoding.Eg Chr1)
  strand_t strand;                   // Mapping Strand
  bs_strand_t bs_strand;             // Bisulfite Strand
  uint64_t text_position;            // Position of the match in the text. Local text (Eg wrt Chr1)
  bool emulated_rc_search;           // Match resulting from a RC-emulated search (using the forward-strand)
  /* Reference-Text */
  uint64_t text_trace_offset;        // Trace-offset in the text-collection
  uint8_t* text;                     // Pointer to the matching-text
  uint64_t text_length;              // Length of the matching-text
  /* Distance/Score */
  uint64_t distance;                 // Distance
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
  /* Text Collection Buffer */
  text_collection_t* text_collection;  // Stores text-traces (candidates/matches/regions/...)
  /* Matches Counters */
  matches_counters_t* counters;        // Global counters
  /* Position Matches */
  vector_t* position_matches;          // Matches (match_trace_t)
  ihash_t* begin_pos_matches;          // Begin position (of the aligned match) in the text-space
  ihash_t* end_pos_matches;            // End position (of the aligned match) in the text-space
  /* Local Matches */
  vector_t* local_matches;             // Local Matches (match_trace_t)
  /* CIGAR */
  vector_t* cigar_vector;              // CIGAR operations storage (cigar_element_t)
  /* Metrics */
  matches_metrics_t metrics;           // Metrics
} matches_t;

/*
 * Setup
 */
matches_t* matches_new();
void matches_configure(matches_t* const matches,text_collection_t* const text_collection);
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

/*
 * Index
 */
void matches_index_rebuild(matches_t* const matches,mm_stack_t* const mm_stack);
void matches_index_clear(matches_t* const matches);

/*
 * Matches Accessors
 */
match_trace_t* matches_get_match_trace_buffer(const matches_t* const matches);
match_trace_t* matches_get_match_trace(const matches_t* const matches,const uint64_t offset);
uint64_t matches_get_num_match_traces(const matches_t* const matches);
void matches_get_clear_match_traces(const matches_t* const matches);

cigar_element_t* match_trace_get_cigar_buffer(const matches_t* const matches,const match_trace_t* const match_trace);
uint64_t match_trace_get_cigar_length(const match_trace_t* const match_trace);
uint64_t match_trace_get_event_distance(const match_trace_t* const match_trace);
int64_t match_trace_get_effective_length(
    matches_t* const matches,
    const uint64_t read_length,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length);

/*
 * Matches Rank Consistency
 */
match_trace_t* matches_get_ranked_match_trace(
    matches_t* const matches,
    select_parameters_t* const select_parameters);

/*
 * Add Matches
 */
bool matches_add_match_trace(
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace,
    mm_stack_t* const mm_stack);
void matches_add_match_trace__preserve_rank(
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace,
    select_parameters_t* const select_parameters,
    const alignment_model_t alignment_model,
    match_trace_t** const match_trace_added,
    bool* const match_added,
    bool* const match_replaced,
    mm_stack_t* const mm_stack);

/*
 * Local Matches
 */
void matches_add_local_match_pending(
    matches_t* const matches,
    match_trace_t* const match_trace);
void matches_add_pending_local_matches(
    matches_t* const matches,
    const locator_t* const locator,
    mm_stack_t* const mm_stack);

/*
 * Matches hints
 */
void matches_hint_allocate_match_trace(
    matches_t* const matches,
    const uint64_t num_matches_trace_to_add);

/*
 * Sorting Matches
 */
void matches_sort_by_distance(matches_t* const matches);
void matches_sort_by_swg_score(matches_t* const matches);
void matches_sort_by_mapq_score(matches_t* const matches);
void matches_sort_by_sequence_name__position(matches_t* const matches);

/*
 * Filters
 */
void matches_filter_by_mapq(matches_t* const matches,const uint8_t mapq_threshold,mm_stack_t* const mm_stack);

/*
 * Display
 */
void matches_print(FILE* const stream,matches_t* const matches);

#endif /* MATCHES_H_ */
