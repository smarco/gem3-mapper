/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_H_
#define APPROXIMATE_SEARCH_H_

#include "approximate_search/approximate_search_metrics.h"
#include "archive/archive_search_parameters.h"
#include "archive/archive.h"
#include "data_structures/pattern.h"
#include "data_structures/interval_set.h"
#include "filtering/region_profile.h"
#include "filtering/filtering_candidates.h"

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
  // Verify Candidates
  asearch_processing_state_candidates_processed = 4,  // Candidates processed
  asearch_processing_state_candidates_verified = 5,   // Candidates verified
} asearch_processing_state_t;
extern const char* asearch_processing_state_label[6];
/*
 * Search Stages
 */
typedef enum {
  asearch_stage_begin = 0,                            // Beginning of the search
  asearch_stage_read_recovery = 1,                    // Read recovery (Reads containing multiple uncalled bases)
  asearch_stage_filtering_adaptive = 2,               // Adaptive filtering search
  asearch_stage_inexact_filtering = 3,                // Inexact Filtering (using approximate-candidate generation)
  asearch_stage_local_alignment = 4,                  // Local-Alignment filtering Search
  asearch_stage_neighborhood = 5,                     // Neighborhood search
  asearch_stage_end = 6,                              // End of the search
} asearch_stage_t;
extern const char* asearch_stage_label[7];

/*
 * Approximate Search
 */
typedef struct {
  /* Index Structures, Pattern & Parameters */
  archive_t* archive;                                      // Archive
  pattern_t pattern;                                       // Search Pattern
  search_parameters_t* search_parameters;                  // Search Parameters
  /* Search State */
  bool emulated_rc_search;                                 // Currently searching on the RC (emulated on the forward strand)
  bool do_quality_search;                                  // Quality search
  asearch_stage_t search_stage;                            // Current Search Stage
  asearch_processing_state_t processing_state;             // Current Processing State
  bool stop_before_neighborhood_search;                    // Stop before Neighborhood Search
  uint64_t max_complete_error;
  uint64_t max_complete_stratum;                           // Maximum complete stratum reached by the search
  /* Search Structures */
  region_profile_t region_profile;                         // Region Profile
  filtering_candidates_t* filtering_candidates;            // Filtering Candidates
  /* Buffered Search Data */
  uint64_t gpu_buffer_fmi_search_offset;
  uint64_t gpu_buffer_fmi_search_total;
  uint64_t gpu_buffer_fmi_decode_offset;
  uint64_t gpu_buffer_fmi_decode_total;
  uint64_t gpu_buffer_align_offset;
  filtering_position_buffered_t* gpu_filtering_positions;  // Filtering-Positions Info
  uint64_t gpu_num_filtering_positions;                    // Total Buffered Filtering-Positions
  filtering_region_buffered_t* gpu_filtering_regions;      // Filtering-Regions Info
  uint64_t gpu_num_filtering_regions;                      // Total Buffered Filtering-Regions
  /* Search Auxiliary Structures (external) */
  text_collection_t* text_collection;                      // Stores text-traces
  interval_set_t* interval_set;                            // Interval Set
  /* Search Metrics */
  approximate_search_metrics_t metrics;                    // Search Metrics
  /* MM */
  mm_stack_t* mm_stack;                                    // MM-Stack
} approximate_search_t;

/*
 * Setup
 */
void approximate_search_init(
    approximate_search_t* const search,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool emulated_rc_search);
void approximate_search_reset(approximate_search_t* const search);
void approximate_search_destroy(approximate_search_t* const search);

/*
 * Memory Injection (Support Data Structures)
 */
void approximate_search_inject_mm_stack(
    approximate_search_t* const search,
    mm_stack_t* const mm_stack);
void approximate_search_inject_interval_set(
    approximate_search_t* const search,
    interval_set_t* const interval_set);
void approximate_search_inject_text_collection(
    approximate_search_t* const search,
    text_collection_t* const text_collection);
void approximate_search_inject_filtering_candidates(
    approximate_search_t* const search,
    filtering_candidates_t* const filtering_candidates,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack);

/*
 * Accessors
 */
void approximate_search_update_mcs(
    approximate_search_t* const search,
    const uint64_t max_complete_stratum);

uint64_t approximate_search_get_num_regions_profile(const approximate_search_t* const search);
uint64_t approximate_search_get_num_decode_candidates(const approximate_search_t* const search);
uint64_t approximate_search_get_num_verify_candidates(const approximate_search_t* const search);

/*
 * Aproximate String Search
 */
void approximate_search(approximate_search_t* const search,matches_t* const matches);

/*
 * Display
 */
void approximate_search_print(FILE* const stream,approximate_search_t* const search);

#endif /* APPROXIMATE_SEARCH_H_ */
