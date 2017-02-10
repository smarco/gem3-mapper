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
 *   Approximate-String-Matching (ASM) main module.
 *   Dispatch the search depending on the search-approach selected
 *   and provides data structures for the search
 */

#ifndef APPROXIMATE_SEARCH_H_
#define APPROXIMATE_SEARCH_H_

#include "archive/search/archive_search_se_parameters.h"
#include "archive/archive.h"
#include "align/pattern/pattern.h"
#include "filtering/region_profile/region_profile.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/candidates/filtering_candidates_buffered.h"
#include "neighborhood_search/nsearch_schedule.h"

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
  asearch_stage_filtering_adaptive = 1,               // Adaptive filtering search
  asearch_stage_filtering_adaptive_finished = 2,      // Adaptive filtering search
  asearch_stage_neighborhood = 3,                     // Neighborhood search
  asearch_stage_local_alignment = 4,                  // Local-Alignment filtering Search
  asearch_stage_end = 5,                              // End of the search
} asearch_stage_t;
extern const char* asearch_stage_label[6];

/*
 * Approximate Search
 */
typedef struct {
  /* Index Structures, Pattern & Parameters */
  archive_t* archive;                                            // Archive
  pattern_t pattern;                                             // Search Pattern
  search_parameters_t* search_parameters;                        // Search Parameters
  /* Search State */
  bool do_quality_search;                                        // Quality search
  asearch_stage_t search_stage;                                  // Current Search Stage
  asearch_processing_state_t processing_state;                   // Current Processing State
  uint64_t current_max_complete_error;                           // Current max-error allowed (adjusted by matches found)
  uint64_t current_max_complete_stratum;                         // Current mcs reached by the search (adjusted by matches found)
  /* Filtering Structures */
  region_profile_t region_profile;                               // Region Profile
  filtering_candidates_t* filtering_candidates;                  // Filtering Candidates
  filtering_candidates_buffered_t filtering_candidates_buffered; // Filtering Candidates Buffered
  /* GPU buffer(s) offsets */
  uint64_t gpu_buffer_fmi_search_offset;                         // FMI-Buffer search offset
  uint64_t gpu_buffer_fmi_search_total;
  uint64_t gpu_buffer_fmi_decode_offset;                         // FMI-Buffer decode offset
  uint64_t gpu_buffer_fmi_decode_total;
  uint64_t gpu_buffer_align_offset;                              // FMI-Buffer align offset
  /* Neighborhood Search Structures */
  nsearch_schedule_t* nsearch_schedule;                          // Nsearch Scheduler
} approximate_search_t;

/*
 * Setup
 */
void approximate_search_init(
    approximate_search_t* const search,
    archive_t* const archive,
    search_parameters_t* const search_parameters);
void approximate_search_reset(approximate_search_t* const search);

void approximate_search_inject_handlers(
    approximate_search_t* const search,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_mm_t* const filtering_candidates_mm,
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm,
    nsearch_schedule_t* const nsearch_schedule,
    mm_stack_t* const mm_region_profile,
    mm_stack_t* const mm_nsearch);

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
