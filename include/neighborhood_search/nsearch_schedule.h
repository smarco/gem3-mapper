/*
 * PROJECT: GEMMapper
 * FILE: nsearch_schedule.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef NSEARCH_SCHEDULE_H_
#define NSEARCH_SCHEDULE_H_

#include "utils/essentials.h"
#include "system/profiler_timer.h"
#include "data_structures/interval_set.h"
#include "fm_index/fm_index.h"
#include "filtering/region_profile.h"
#include "neighborhood_search/dp_matrix.h"
#include "neighborhood_search/nsearch_operation.h"
#include "approximate_search/approximate_search.h"
#include "matches/matches.h"

/*
 * Enumaration Mode
 */
//#define NSEARCH_ENUMERATE

/*
 * Profile
 */
#ifdef GEM_DEBUG
  #define NSEARCH_PROF_ADD_NODE(nsearch_schedule)      ++((nsearch_schedule)->profile.ns_nodes)
  #define NSEARCH_PROF_CLOSE_NODE(nsearch_schedule)    ++((nsearch_schedule)->profile.ns_nodes_closed)
  #define NSEARCH_PROF_ADD_SOLUTION(nsearch_schedule)  ++((nsearch_schedule)->profile.ns_nodes_success)
  #define NSEARCH_PROF_ACCOUNT_DEPTH(nsearch_schedule,depth) \
    COUNTER_ADD(&(nsearch_schedule)->profile.ns_nodes_closed_depth,depth);
#else
  #define NSEARCH_PROF_ADD_NODE(nsearch_schedule)
  #define NSEARCH_PROF_CLOSE_NODE(nsearch_schedule)
  #define NSEARCH_PROF_ADD_SOLUTION(nsearch_schedule)
  #define NSEARCH_PROF_ACCOUNT_DEPTH(nsearch_schedule,depth)
#endif

/*
 * Neighborhood Search Metric
 */
typedef enum {
  nsearch_model_hamming,
  nsearch_model_levenshtein
} nsearch_model_t;

/*
 * Neighborhood Search Schedule
 */
typedef struct {
  uint64_t ns_nodes;
  uint64_t ns_nodes_mtable;
  uint64_t ns_nodes_success;
  uint64_t ns_nodes_closed;
  uint64_t ns_nodes_fail_optimize;
  gem_counter_t ns_nodes_closed_depth;
  gem_timer_t ns_timer;
} nsearch_schedule_profile_t;
typedef struct {
  // Search Structures
  uint64_t search_id;                         // Search ID
  approximate_search_t* search;               // ASM-Search
  matches_t* matches;                         // Matches
  // Search Parameters
  nsearch_model_t nsearch_model;              // Search error model
  uint8_t* key;
  uint64_t key_length;
  uint64_t max_error;
  // Scheduler Operations
  nsearch_operation_t* pending_searches;      // Pending search operations
  uint64_t num_pending_searches;              // Total pending operations
  // Misc
  nsearch_schedule_profile_t profile;         // Profiler
  mm_stack_t* mm_stack;                       // MM
} nsearch_schedule_t;

/*
 * Setup
 */
void nsearch_schedule_init(
    nsearch_schedule_t* const nsearch_schedule,
    const nsearch_model_t nsearch_model,
    approximate_search_t* const search,
    matches_t* const matches);

/*
 * Schedule the search
 */
void nsearch_schedule_search(nsearch_schedule_t* const nsearch_schedule);
void nsearch_schedule_search_preconditioned(nsearch_schedule_t* const nsearch_schedule);

/*
 * Display
 */
void nsearch_schedule_print(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule);
void nsearch_schedule_print_pretty(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule);
void nsearch_schedule_print_profile(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule);

#endif /* NSEARCH_SCHEDULE_H_ */
