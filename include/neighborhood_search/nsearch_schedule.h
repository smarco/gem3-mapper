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

#ifndef NSEARCH_SCHEDULE_H_
#define NSEARCH_SCHEDULE_H_

#include "utils/essentials.h"
#include "archive/archive.h"
#include "archive/search/archive_search_se_parameters.h"
#include "align/pattern/pattern.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/region_profile/region_profile.h"
#include "matches/matches.h"
#include "neighborhood_search/nsearch_operation.h"

/*
 * Enumeration Mode
 */
//#define NSEARCH_ENUMERATE

/*
 * Profile
 */
#ifdef GEM_PROFILE
  #define NSEARCH_PROF_NODE(nsearch_schedule,matches_found) \
    ++((nsearch_schedule)->profile.ns_nodes); \
    if (matches_found==0) { \
      ++((nsearch_schedule)->profile.ns_nodes_fail); \
    } else { \
      ++((nsearch_schedule)->profile.ns_nodes_success); \
    }
#else
  #define NSEARCH_PROF_NODE(nsearch_schedule,matches_found)
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
  /* Nodes */
  uint64_t ns_nodes;
  uint64_t ns_nodes_success;
  uint64_t ns_nodes_fail;
} nsearch_schedule_profile_t;
typedef struct {
  // Index Structures & Pattern
  uint64_t search_id;                           // Search ID
  archive_t* archive;                           // Archive
  pattern_t* pattern;                           // Search Pattern
  region_profile_t* region_profile;             // Region Profile
  filtering_candidates_t* filtering_candidates; // Filtering Candidates
  matches_t* matches;                           // Matches
  // Search Parameters
  search_parameters_t* search_parameters;       // Search Parameters
  nsearch_model_t nsearch_model;                // Search error model
  uint64_t max_error;                           // Max-error search allowed
  // Scheduler Operations
  nsearch_operation_t* pending_searches;        // Pending search operations
  uint64_t num_pending_searches;                // Total pending operations
  bool quick_abandon;                           // Quick abandon search
  // Misc
  nsearch_schedule_profile_t profile;           // Profiler
  mm_allocator_t* mm_allocator;                 // MM-Allocator
} nsearch_schedule_t;

/*
 * Setup
 */
void nsearch_schedule_init(
    nsearch_schedule_t* const nsearch_schedule,
    const nsearch_model_t nsearch_model,
    const uint64_t max_complete_error,
    archive_t* const archive,
    pattern_t* const pattern,
    region_profile_t* const region_profile,
    search_parameters_t* const search_parameters,
    filtering_candidates_t* const filtering_candidates,
    matches_t* const matches);
void nsearch_schedule_inject_mm(
    nsearch_schedule_t* const nsearch_schedule,
    mm_allocator_t* const mm_allocator);

/*
 * Schedule the search
 */
void nsearch_schedule_search(nsearch_schedule_t* const nsearch_schedule);
void nsearch_schedule_search_preconditioned(nsearch_schedule_t* const nsearch_schedule);

/*
 * Utils
 */
uint64_t nsearch_schedule_compute_min_error(
    nsearch_schedule_t* const nsearch_schedule);

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
