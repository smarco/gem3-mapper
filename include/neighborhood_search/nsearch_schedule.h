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
#include "neighborhood_search/nsearch_levenshtein_state.h"
#include "neighborhood_search/dp_matrix.h"

/*
 * Neighborhood Search Metric
 */
typedef enum {
  nsearch_model_hamming,
  nsearch_model_levenshtein
} nsearch_model_t;

/*
 * Neighborhood Search Operation
 */
typedef enum { direction_forward, direction_backward } search_direction_t;
typedef struct {
  // Pattern chunk
  uint64_t begin;
  uint64_t end;
  // Search direction
  search_direction_t search_direction;
  // Error
  uint64_t max_local_error;
  uint64_t min_local_error;
  uint64_t max_global_error;
  uint64_t min_global_error;
  uint64_t max_text_length;
  // Search State
  nsearch_levenshtein_state_t nsearch_state;
} nsearch_operation_t;

/*
 * Neighborhood Search Schedule
 */
typedef struct {
  // FM-Index
  fm_index_t* fm_index;
  // Search Parameters
  nsearch_model_t nsearch_model;
  uint8_t* key;
  uint64_t key_length;
  region_profile_t* region_profile;
  uint64_t max_error;
  uint64_t max_text_length;
  // Scheduler Progress
  nsearch_operation_t* pending_searches; // Pending search operations
  uint64_t num_pending_searches;         // Total pending operations
  uint64_t search_id;                    // Search ID
  // Output results
  interval_set_t* intervals_result;
  // Profiler/Stats
  uint64_t ns_nodes;
  uint64_t ns_nodes_mtable;
  uint64_t ns_nodes_success;
  uint64_t ns_nodes_closed;
  uint64_t ns_nodes_fail_optimize;
  gem_counter_t ns_nodes_closed_depth;
  gem_timer_t ns_timer;
  // MM
  mm_stack_t* mm_stack;
  // Debug
  char* search_string;
} nsearch_schedule_t;

/*
 * Setup
 */
void nsearch_schedule_init(
    nsearch_schedule_t* const nsearch_schedule,const nsearch_model_t nsearch_model,
    fm_index_t* const fm_index,region_profile_t* const region_profile,
    uint8_t* const key,const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack);

/*
 * Schedule the search
 */
void nsearch_schedule_search(nsearch_schedule_t* const nsearch_schedule);
void nsearch_schedule_search_preconditioned(nsearch_schedule_t* const nsearch_schedule);

/*
 * Display
 */
void nsearch_schedule_print(FILE* const stream,nsearch_schedule_t* const nsearch_schedule);
void nsearch_schedule_print_pretty(FILE* const stream,nsearch_schedule_t* const nsearch_schedule);
void nsearch_schedule_print_profile(FILE* const stream,nsearch_schedule_t* const nsearch_schedule);

void nsearch_schedule_print_search_string(FILE* const stream,nsearch_schedule_t* const nsearch_schedule);

#endif /* NSEARCH_SCHEDULE_H_ */
