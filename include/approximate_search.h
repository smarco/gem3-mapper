/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_H_
#define APPROXIMATE_SEARCH_H_

#include "search_parameters.h"
#include "archive.h"
#include "pattern.h"
#include "region_profile.h"
#include "interval_set.h"
#include "filtering_candidates.h"
#include "mapper_profile.h"

/*
 * Debug
 */
#define DEBUG_SEARCH_STATE          GEM_DEEP_DEBUG

/*
 * Search States
 */
typedef enum {
  // Null
  asearch_processing_state_begin = 0,                 // Begin
  // Region Profile
  asearch_processing_state_region_partitioned = 1,    // Region Partition Performed
  asearch_processing_state_region_profiled = 2,       // Region Profile Performed
  asearch_processing_state_no_regions = 3,            // Region Profile Performed (No regions were found)
  asearch_processing_state_exact_matches = 4,         // Region Profile Performed (Exact Matches; One maximum region was found)
  // Verify Candidates
  asearch_processing_state_candidates_processed = 5,  // Candidates processed
  asearch_processing_state_candidates_verified = 6,   // Candidates verified
} asearch_processing_state_t;
extern const char* asearch_processing_state_label[7];
/*
 * Search Stages
 */
typedef enum {
  asearch_stage_begin = 0,                       // Beginning of the search
  asearch_stage_read_recovery = 1,               // Read recovery (Reads containing multiple uncalled bases)
  asearch_stage_filtering_adaptive = 2,          // Adaptive filtering search
  asearch_stage_filtering_boost = 4,             // Boost filtering
  asearch_stage_inexact_filtering = 5,           // Inexact Filtering (using approximate-candidate generation)
  asearch_stage_unbounded_alignment = 6,         // Unbounded filtering Search (unbounded alignment)
  asearch_stage_neighborhood = 7,                // Neighborhood search
  asearch_stage_end = 8,                         // End of the search
} asearch_stage_t;
extern const char* asearch_stage_label[9];

/*
 * Approximate Search
 */
typedef struct {
  /* Index Structures, Pattern & Parameters */
  archive_t* archive;                                     // Archive
  pattern_t pattern;                                      // Search Pattern
  as_parameters_t* as_parameters; // Search Parameters (Evaluated to read-length)
  /* Search State */
  bool emulated_rc_search;                                // Currently searching on the RC (emulated on the forward strand)
  bool do_quality_search;                                 // Quality search
  asearch_stage_t search_stage;                           // Current Search Stage
  asearch_processing_state_t processing_state;            // Current Processing State
  bool stop_before_neighborhood_search;                   // Stop before Neighborhood Search
  uint64_t max_complete_error;
  uint64_t max_complete_stratum;                          // Maximum complete stratum reached by the search
  uint64_t max_matches_reached;                           // Quick abandon due to maximum matches found
  uint64_t lo_exact_matches;                              // Interval Lo (Exact matching)
  uint64_t hi_exact_matches;                              // Interval Hi (Exact matching)
  /* Search Structures */
  region_profile_t region_profile;                        // Region Profile
  filtering_candidates_t* filtering_candidates;           // Filtering Candidates
  /* Buffer Offsets */
  uint64_t gpu_buffer_fmi_search_offset;
  uint64_t gpu_buffer_fmi_search_total;
  uint64_t gpu_buffer_fmi_decode_offset;
  uint64_t gpu_buffer_fmi_decode_total;
  uint64_t gpu_buffer_align_offset;
  uint64_t gpu_buffer_align_total;
  /* Search Auxiliary Structures (external) */
  text_collection_t* text_collection;                   // Stores text-traces
  interval_set_t* interval_set;                         // Interval Set
  /* MM */
  mm_stack_t* mm_stack;                                 // MM-Stack
} approximate_search_t;

/*
 * Setup
 */
void approximate_search_init(
    approximate_search_t* const search,archive_t* const archive,
    as_parameters_t* const as_parameters,const bool emulated_rc_search);
void approximate_search_configure(
    approximate_search_t* const search,filtering_candidates_t* const filtering_candidates,
    text_collection_t* text_collection,interval_set_t* const interval_set,mm_stack_t* const mm_stack);
void approximate_search_reset(approximate_search_t* const search);
void approximate_search_destroy(approximate_search_t* const search);

/*
 * Accessors
 */
uint64_t approximate_search_get_num_filtering_candidates(const approximate_search_t* const search);
uint64_t approximate_search_get_num_exact_filtering_candidates(const approximate_search_t* const search);
void approximate_search_update_mcs(approximate_search_t* const search,const uint64_t max_complete_stratum);

/*
 * Modifiers
 */
void approximate_search_hold_verification_candidates(approximate_search_t* const search);
void approximate_search_release_verification_candidates(approximate_search_t* const search);

/*
 * Aproximate String Search
 */
void approximate_search(approximate_search_t* const search,matches_t* const matches);

/*
 * Display
 */
void approximate_search_print(FILE* const stream,approximate_search_t* const search);

#endif /* APPROXIMATE_SEARCH_H_ */
