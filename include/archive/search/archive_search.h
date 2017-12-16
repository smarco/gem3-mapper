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
 *   Archive-Search main module & support data structures
 */

#ifndef ARCHIVE_SEARCH_H_
#define ARCHIVE_SEARCH_H_

#include "utils/essentials.h"
#include "text/sequence.h"
#include "archive/archive.h"
#include "approximate_search/approximate_search.h"
#include "matches/classify/matches_classify.h"
#include "mapper/mapper_stats.h"

/*
 * Archive Search State
 */
typedef enum {
  archive_search_pe_state_begin = 0,                     // Beginning of the search
  archive_search_pe_state_search_end1 = 1,               // Generate candidates for end1
  archive_search_pe_state_search_end2 = 2,               // Generate candidates for end2
  archive_search_pe_state_find_pairs = 3,                // Cross-link matches from both ends
  archive_search_pe_state_end = 4                        // End of the current workflow
} archive_search_pe_state_t;
extern const char* archive_search_pe_state_label[7];

/*
 * Archive Search
 */
typedef struct {
  /* Archive */
  archive_t* archive;                        // Archive
  /* PE State */
  archive_search_pe_state_t pe_search_state; // Search State (only end/1 updated)
  bool searched;                             // Archive performed
  /* Parameters */
  search_parameters_t search_parameters;     // Search parameters
  sequence_t* sequence;                      // Input Sequence
  bool buffered_search;                      // Buffered Search
  /* Approximate Search */
  approximate_search_t approximate_search;   // Approximate Search State
  /* Stats */
  mapper_stats_t* mapper_stats;              // Mapping statistics
  /* MM */
  mm_allocator_t* mm_allocator;              // MM-Allocator
} archive_search_t;

/*
 * Setup
 */
void archive_search_se_new(
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    archive_search_t** const archive_search);
void archive_search_pe_new(
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);
void archive_search_destroy(
    archive_search_t* const archive_search);
void archive_search_delete(
    archive_search_t* const archive_search);

/*
 * Prepare Search
 */
void archive_search_prepare_sequence(
    archive_search_t* const archive_search,
    sequence_t* const sequence);
void archive_search_inject_handlers(
    archive_search_t* const archive_search,
    archive_t* const archive,
    filtering_candidates_t* const filtering_candidates,
    nsearch_schedule_t* const nsearch_schedule,
    mapper_stats_t* const mapper_stats,
    mm_allocator_t* const mm_allocator);

/*
 * Accessors
 */
bool archive_search_finished(const archive_search_t* const archive_search);

uint64_t archive_search_get_num_regions_profile(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_decode_candidates(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_kmer_filter_candidates(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_bpm_distance_candidates(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_bpm_align_candidates(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_bpm_align_canonical_candidates(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_bpm_align_candidate_tiles_length(const archive_search_t* const archive_search);

/*
 * Errors
 */
#define GEM_ERROR_ARCHIVE_SEARCH_INDEX_COMPLEMENT_REQUIRED "Archive Search. Explicit indexed complement required"

#endif /* ARCHIVE_SEARCH_H_ */
